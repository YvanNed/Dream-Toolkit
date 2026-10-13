"""
Sleep macrostructure: the metric registry, the computations and the figures of tool 5
(5_compute_sleep_macrostructure_voila.ipynb), shared with tool 5bis (registry, reference ranges) and tool 0
(step 4). Tool-0 redesign, tools/plan_tool0_redesign.md.

MOVED VERBATIM from tool 5: its constants and helpers (cell 1), computations (cell 3) and figures / HTML
(cell 4), plus `compute_recording`, the computation part of its `process_participant` lifted unchanged. The
TXT -> CSV -> XML scored-event chain is edf_io_lib's (it was tool 5's copy, identical). Every metric is
declared once, in METRICS: add a metric there, never ad hoc. The policies (OTHER_STAGE_POLICY,
EVENT_COUNT_WINDOW) are set here, once for every tool.
"""
import os
import re
import io
import json
import difflib
import datetime
import unicodedata
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import mne

from edf_io_lib import (read_edf_start_datetime, event_companion_paths, _events_df_from_txt, _events_df_from_csv,
                        _events_df_from_xml, _events_df_from_manual, augment_remap_for_manual, load_events,
                        load_event_remap)

try:
    import openpyxl                     # noqa: F401  (pandas needs it to read/write .xlsx)
    HAS_OPENPYXL = True
except Exception:
    HAS_OPENPYXL = False

# --- sleep stages ------------------------------------------------------------
AASM_STAGES = ['W', 'N1', 'N2', 'N3', 'R']
SLEEP_STAGES = ['N1', 'N2', 'N3', 'R']      # everything counted as sleep

# MT (movement time), custom stages and unknown labels are "other": they follow YASA's convention -
# counted in time in bed (and inside the sleep period when they fall there) but excluded from both
# total sleep time and WASO. Set to 'wake' to count them as wake instead (WASO when inside the
# sleep period); the choice is recorded in the outputs.
OTHER_STAGE_POLICY = 'exclude'
# OTHER_STAGE_POLICY = 'wake'

# Window over which the scored events are counted for the per-hour indices: 'tib' = every event
# between lights-off and lights-on (denominator TST, the clinical convention); 'sleep_epochs' =
# only events whose onset falls in a sleep epoch.
EVENT_COUNT_WINDOW = 'tib'
# EVENT_COUNT_WINDOW = 'sleep_epochs'

DEFAULT_EPOCH_S = 30

# Pipeline cost units for the per-participant progress bar
COST_HEADER, COST_HYPNO, COST_LIGHTS, COST_EVENTS, COST_SPO2, COST_WRITE = 2, 3, 20, 10, 20, 5

# --- metric registry ----------------------------------------------------------
# (key, label, unit, group, one-line definition) - the single source for the long TSV, the
# glossary, the workbook and the report. Order = column order of the global table.
METRICS = [
    ('tib_min', 'Time in bed (TIB)', 'min', 'macrostructure',
     'Lights-off to lights-on: number of epochs in bed x epoch length.'),
    ('sol_min', 'Sleep onset latency (SOL)', 'min', 'macrostructure',
     'Lights-off to sleep onset (rule chosen in Section 2; default = first non-W epoch).'),
    ('spt_min', 'Sleep period time (SPT)', 'min', 'macrostructure',
     'Sleep onset to the last sleep epoch, both included.'),
    ('tst_min', 'Total sleep time (TST)', 'min', 'macrostructure',
     'Number of N1 + N2 + N3 + R epochs within the sleep period x epoch length.'),
    ('waso_min', 'Wake after sleep onset (WASO)', 'min', 'macrostructure',
     'Number of W epochs within the sleep period x epoch length.'),
    ('other_min', 'Other epochs within SPT', 'min', 'macrostructure',
     'MT / custom / unknown epochs within the sleep period (neither sleep nor wake; see policy).'),
    ('se_pct', 'Sleep efficiency (SE)', '%', 'macrostructure', '100 x TST / TIB.'),
    ('sme_pct', 'Sleep maintenance efficiency (SME)', '%', 'macrostructure', '100 x TST / SPT.'),
    ('lat_n1_min', 'Latency to N1', 'min', 'macrostructure',
     'Sleep onset to the first N1 epoch (0 under the default onset rule when onset is N1).'),
    ('lat_n2_min', 'Latency to N2', 'min', 'macrostructure', 'Sleep onset to the first N2 epoch.'),
    ('lat_n3_min', 'Latency to N3', 'min', 'macrostructure', 'Sleep onset to the first N3 epoch.'),
    ('lat_rem_min', 'REM latency', 'min', 'macrostructure', 'Sleep onset to the first R epoch.'),
    ('n1_min', 'N1 duration', 'min', 'macrostructure', 'Number of N1 epochs within SPT x epoch length.'),
    ('n2_min', 'N2 duration', 'min', 'macrostructure', 'Number of N2 epochs within SPT x epoch length.'),
    ('n3_min', 'N3 duration', 'min', 'macrostructure', 'Number of N3 epochs within SPT x epoch length.'),
    ('rem_min', 'REM duration', 'min', 'macrostructure', 'Number of R epochs within SPT x epoch length.'),
    ('n1_pct', 'N1 (% of TST)', '%', 'macrostructure', '100 x N1 duration / TST.'),
    ('n2_pct', 'N2 (% of TST)', '%', 'macrostructure', '100 x N2 duration / TST.'),
    ('n3_pct', 'N3 (% of TST)', '%', 'macrostructure', '100 x N3 duration / TST.'),
    ('rem_pct', 'REM (% of TST)', '%', 'macrostructure', '100 x REM duration / TST.'),
    ('n_apnea_obstructive', 'Obstructive apneas', 'count', 'respiratory',
     'Scored obstructive apneas within the counting window.'),
    ('n_apnea_central', 'Central apneas', 'count', 'respiratory',
     'Scored central apneas within the counting window.'),
    ('n_apnea_mixed', 'Mixed apneas', 'count', 'respiratory',
     'Scored mixed apneas within the counting window.'),
    ('n_hypopnea', 'Hypopneas', 'count', 'respiratory',
     'Scored hypopneas (any canonical label starting with "hypopnea") within the counting window.'),
    ('oai', 'Obstructive apnea index (OAI)', 'events/h', 'respiratory', 'Obstructive apneas / TST (h).'),
    ('cai', 'Central apnea index (CAI)', 'events/h', 'respiratory', 'Central apneas / TST (h).'),
    ('mai', 'Mixed apnea index (MAI)', 'events/h', 'respiratory', 'Mixed apneas / TST (h).'),
    ('hi', 'Hypopnea index (HI)', 'events/h', 'respiratory', 'Hypopneas / TST (h).'),
    ('ahi', 'Apnea-hypopnea index (AHI)', 'events/h', 'respiratory',
     '(all apneas + hypopneas) / TST (h) = OAI + CAI + MAI + HI.'),
    ('n_plm', 'Periodic limb movements', 'count', 'plm',
     'Scored PLM events (canonical "plm") within the counting window.'),
    ('plm_index', 'PLM index', 'events/h', 'plm', 'PLM events / TST (h).'),
    ('n_arousal', 'Arousals', 'count', 'arousal',
     'Scored arousals (any canonical label starting with "arousal") within the counting window.'),
    ('arousal_index', 'Arousal index', 'events/h', 'arousal', 'Arousals / TST (h).'),
    ('n_desaturation', 'SpO2 desaturations', 'count', 'odi',
     'Scored desaturations within the counting window (depth >= threshold when the depth is known).'),
    ('odi', 'Oxygen desaturation index (ODI)', 'events/h', 'odi', 'Desaturations / TST (h).'),
    ('t90_min', 'T90', 'min', 'spo2', 'Time with SpO2 < threshold (default 90 %) between lights-off and lights-on.'),
    ('t90_pct_tst', 'T90 (% of TST)', '%', 'spo2', '100 x T90 / TST.'),
    ('t88_sleep_min', 'T88 during sleep', 'min', 'spo2',
     'Time with SpO2 <= 88 % during the sleep epochs (N1/N2/N3/R) between lights-off and lights-on '
     '(ICSD-3-TR sleep-related hypoxemia, criterion A: >= 5 min).'),
    ('spo2_mean_pct', 'Mean SpO2', '%', 'spo2', 'Mean of the valid SpO2 samples in bed.'),
    ('spo2_nadir_pct', 'SpO2 nadir', '%', 'spo2', 'Minimum of the valid SpO2 samples in bed.'),
    ('spo2_artifact_pct', 'SpO2 artefact', '%', 'spo2',
     'Percentage of SpO2 samples in bed discarded as artefact (<= 0 or below the floor).'),
]
METRIC_KEYS = [m[0] for m in METRICS]
METRIC_INFO = {m[0]: {'label': m[1], 'unit': m[2], 'group': m[3], 'definition': m[4]} for m in METRICS}

# Provenance written with every participant (text columns of the global table)
PROVENANCE_KEYS = [
    'lights_source', 'lights_off_clock', 'lights_on_clock', 'lights_off_epoch', 'lights_on_epoch',
    'sleep_onset_rule', 'sleep_onset_epoch', 'other_stage_policy', 'event_count_window',
    'epoch_length_s', 'n_epochs_hypno', 'n_epochs_in_bed', 'n_mt_custom', 'event_source',
    'n_events_total', 'n_unmapped_events', 'n_events_outside_window', 'has_spo2', 'odi_threshold_pct',
    'odi_rule', 'spo2_channel', 'light_channel', 'recording_start', 'edf_duration_min', 'n_warnings']

# --- indicative normal ranges (healthy adults) --------------------------------
# INDICATIVE ONLY: ranges for healthy adults recorded in a sleep laboratory; they depend on age
# (N3, REM and SE fall, WASO and N1 rise with age), on the first-night effect and on the scoring
# rules. Where each value comes from: REFERENCE_SOURCES / REFERENCES below (also shown in the
# collapsible section at the end of the notebook and in the database report). Edit freely; they are
# only used to shade the figures of the database report and to annotate the glossary. (low, high, note)
REFERENCE_RANGES = {
    'sol_min':        (0, 30,   '~15-20 min at every adult age, < 10 min longer at 80 than at 20 y (Ohayon 2004); '
                                '> 30 min suggests difficulty initiating sleep'),
    'tst_min':        (360, 480, '~7.5 h at 20-25 y, ~6.7 h at 50 y, ~6 h after 75 y in the laboratory '
                                 '(-10 min per decade, Ohayon 2004)'),
    'waso_min':       (0, 30,   'young adults; +10 min per decade from 30 y (~40 min at 55 y, ~55-75 min at '
                                '65-75 y in healthy subjects, Ohayon 2004) - > 30 min is not abnormal after 50 y'),
    'se_pct':         (85, 100, '~95 % at 20 y, ~90 % at 45, ~85 % at 65, ~80 % after 75 y (-3 % per decade '
                                'after 40 y, Ohayon 2004); lower on a first laboratory night'),
    'sme_pct':        (90, 100, 'sleep continuity once asleep (convention)'),
    'lat_rem_min':    (50, 110, '~90 min at 25 y, ~70 at 65, ~55 at 90 y (Ohayon 2004); <= 15 min = sleep-onset '
                                'REM period (ICSD-3-TR, narcolepsy); < 50 min sometimes associated with depression; '
                                'prolonged by REM-suppressing medication'),
    'n1_pct':         (2, 8,    '~6 % in young adults, ~7 % and up to ~10 % in older adults (Ohayon 2004); '
                                'high N1 = fragmented sleep'),
    'n2_pct':         (45, 58,  '~50 % at 25 y, ~55 % at 65, ~58 % at 90 y (Ohayon 2004)'),
    'n3_pct':         (13, 23,  'young adults (~18-20 %); ~12 % at 45 y, ~5-10 % after 65 y (-2 % per decade, '
                                'Ohayon 2004); depends on the amplitude criterion of the scorer'),
    'rem_pct':        (18, 25,  '~21-22 % in young adults, ~18-19 % after 70 y (Ohayon 2004)'),
    'ahi':            (0, 5,    'AASM severity: < 5 normal, 5-15 mild, 15-30 moderate, > 30 severe'),
    'plm_index':      (0, 15,   'PLMD threshold (ICSD-3-TR): > 15/h in adults (> 5/h in children)'),
    'arousal_index':  (0, 15,   'convention: ~10/h in young adults, higher with age; > 15-20/h = fragmented sleep'),
    'odi':            (0, 5,    'convention: same severity bands as the AHI; state the threshold (3 % or 4 %)'),
    't90_pct_tst':    (0, 10,   'research convention (CT90), not an ICSD-3-TR criterion'),
    't88_sleep_min':  (0, 5,    'ICSD-3-TR sleep-related hypoxemia (G47.36), criterion A: >= 5 min with SpO2 <= 88 % '
                                'during sleep (adults; <= 90 % in children); criterion B (not explained by '
                                'hypoventilation / OSA) needs a clinical evaluation'),
    'spo2_mean_pct':  (94, 100, 'convention'),
    'spo2_nadir_pct': (88, 100, 'convention - not an ICSD-3-TR criterion: an isolated nadir <= 88 % is not '
                                'hypoxemia (see T88)'),
}

# Where each indicative value comes from: (metric, indicative value, source, remark).
# Ohayon 2004 values are read on the fitted age curves of its Figure 1 (approximate).
REFERENCE_SOURCES = [
    ('SOL', '< 30 min (typically 15-20)', 'Ohayon 2004',
     'fitted curve ~16 min at 20 y, ~19 min at 90 y; < 10 min increase between 20 and 80 y'),
    ('TST', '360-480 min', 'Ohayon 2004', '~450 min at 20-25 y, ~400 at 50 y, ~360 after 75 y; -10 min per decade'),
    ('WASO', '< 30 min', 'Ohayon 2004',
     'valid up to ~45 y; +10 min per decade from 30 y (~40 min at 55 y, ~55-75 min at 65-75 y)'),
    ('SE', '> 85 %', 'Ohayon 2004',
     '~95 % at 20 y, ~90 % at 45, ~85 % at 65, ~80 % after 75 y; the only parameter still decreasing after 60 y'),
    ('SME', '> 90 %', 'convention', 'no normative reference'),
    ('REM latency', '50-110 min', 'Ohayon 2004; ICSD-3-TR',
     '~90 min at 25 y, ~70 at 65, ~55 at 90 y; <= 15 min = sleep-onset REM period (narcolepsy, ICSD-3-TR)'),
    ('N1 %', '2-8 %', 'Ohayon 2004; Carskadon & Dement',
     '~6 % in young adults, ~7 % (up to ~10 %) in older adults; textbooks give 2-5 %'),
    ('N2 %', '45-58 %', 'Ohayon 2004', '~50 % at 25 y, ~55 % at 65, ~58 % at 90 y'),
    ('N3 %', '13-23 %', 'Ohayon 2004',
     'young adults (~18-20 %); ~12 % at 45 y, ~5-10 % after 65 y; Ohayon uses R&K S3 + S4, not AASM N3'),
    ('REM %', '18-25 %', 'Ohayon 2004', '~21-22 % in young adults, ~18-19 % after 70 y'),
    ('AHI (severity bands)', '< 5/h; 5-15 / 15-30 / > 30', 'AASM 1999', 'normal / mild / moderate / severe'),
    ('PLM index', '< 15/h (adults)', 'ICSD-3-TR', 'PLMD criterion: > 15/h in adults, > 5/h in children'),
    ('Arousal index', '< 15/h', 'convention', 'no normative reference; rises with age'),
    ('ODI', '< 5/h', 'convention',
     'mirrors the AHI bands; state the desaturation threshold (3 % or 4 %, AASM scoring manual)'),
    ('T90', '< 10 % of TST', 'convention', 'research convention (CT90); not an ICSD-3-TR criterion'),
    ('T88 during sleep', '< 5 min', 'ICSD-3-TR',
     'sleep-related hypoxemia (G47.36), criterion A: SpO2 <= 88 % (adults; <= 90 % in children) for >= 5 min '
     'during sleep; criterion B (not fully explained by hypoventilation, OSA or another breathing disorder) '
     'needs a clinical evaluation'),
    ('Mean SpO2', '>= 94 %', 'convention', 'no normative reference'),
    ('SpO2 nadir', '>= 88 %', 'convention', 'not an ICSD-3-TR criterion: an isolated nadir is not hypoxemia'),
]
REFERENCES = [   # (short name, citation, DOI or '' - the DOI is appended as 'doi:<link>')
    ('Ohayon 2004', 'Ohayon MM, Carskadon MA, Guilleminault C, Vitiello MV. Meta-analysis of quantitative sleep '
                    'parameters from childhood to old age in healthy individuals: developing normative sleep '
                    'values across the human lifespan. Sleep. 2004;27(7):1255-73.',
     '10.1093/sleep/27.7.1255'),
    ('ICSD-3-TR', 'American Academy of Sleep Medicine. International classification of sleep disorders. 3rd ed., '
                  'text revision (ICSD-3-TR). Darien (IL): American Academy of Sleep Medicine; 2023.', ''),
    ('AASM 1999', 'American Academy of Sleep Medicine Task Force. Sleep-related breathing disorders in adults: '
                  'recommendations for syndrome definition and measurement techniques in clinical research. '
                  'Sleep. 1999;22(5):667-89.', '10.1093/sleep/22.5.667'),
    ('AASM scoring manual', 'Troester MM, Quan SF, Berry RB, et al. The AASM manual for the scoring of sleep and '
                            'associated events: rules, terminology and technical specifications. Version 3. '
                            'Darien (IL): American Academy of Sleep Medicine; 2023.', ''),
    ('Carskadon & Dement', 'Carskadon MA, Dement WC. Normal human sleep: an overview. In: Kryger MH, Roth T, '
                           'Dement WC, editors. Principles and practice of sleep medicine. 6th ed. '
                           'Philadelphia: Elsevier; 2017. p. 15-24.', ''),
]
AHI_SEVERITY_BANDS = [(0, 5, 'normal'), (5, 15, 'mild'), (15, 30, 'moderate'), (30, None, 'severe')]


def norm(p):
    """Normalised path/id string, used on BOTH sides of every comparison (Windows safety)."""
    return os.path.normcase(str(p))


def strip_accents(text):
    """'Luminosité' -> 'Luminosite' (accent-insensitive column matching)."""
    return ''.join(c for c in unicodedata.normalize('NFKD', str(text)) if not unicodedata.combining(c))


def parse_custom_field(text):
    """Parse the comma-separated 'Custom stages' field into a clean, de-duplicated list."""
    seen = []
    for tok in str(text).split(','):
        tok = tok.strip()
        if tok and tok not in seen:
            seen.append(tok)
    return seen


def load_custom_stages(folder):
    """Read config_param/custom_stages.json under the study folder ([] when absent)."""
    if not folder:
        return []
    path = Path(folder) / 'config_param' / 'custom_stages.json'
    if not path.exists():
        return []
    try:
        with open(path, encoding='utf-8') as f:
            data = json.load(f)
        stages = data.get('custom_stages', []) if isinstance(data, dict) else []
        return [str(s) for s in stages]
    except Exception:
        return []


def load_subject_info(path):
    """Read the optional participant-info table (.csv/.tsv/.xlsx). Returns (DataFrame, error)."""
    if not path:
        return None, ''
    try:
        if norm(path).endswith(('.xlsx', '.xls')):
            if not HAS_OPENPYXL:
                return None, 'openpyxl is not installed - save the table as .csv or .tsv'
            df = pd.read_excel(path)
        else:
            sep = '\t' if norm(path).endswith('.tsv') else None
            df = pd.read_csv(path, sep=sep, engine='python')
        return df, ''
    except Exception as exc:
        return None, str(exc)


# =============================================================================
# Core computation - EDF header, hypnogram, lights, macrostructure, events, SpO2, checks
# (pure functions: no widget is touched here)
# =============================================================================

# ---------------------------------------------------------------------------
# EDF header (header-only binary read - the signal is never loaded here)
# ---------------------------------------------------------------------------
def read_edf_header_info(edf_path):
    """Read the fixed EDF header and the per-channel labels without touching the signal.

    Same byte offsets as tool 1's custom parser: 168 start date 'dd.mm.yy', 176 start time
    'hh.mm.ss', 184 header size, 236 number of data records, 244 record duration (s), 252 number
    of channels, then 16 bytes x n labels ... 8 bytes x n samples-per-record. The recording
    duration is `n_records x record_duration`; when the record count is -1 (streaming export) it
    is recovered from the file size. edfio is NOT used on purpose: it fails on files whose EDF+
    annotation record is empty (tools/test_data/73.edf). Raises on any failure (fatal)."""
    edf_path = Path(edf_path)
    with open(edf_path, 'rb') as f:
        head = f.read(256)
        if len(head) < 256:
            raise ValueError('file shorter than the 256-byte EDF header')
        date_str = head[168:176].decode('ascii', 'replace').strip()
        time_str = head[176:184].decode('ascii', 'replace').strip()
        header_bytes = int(head[184:192].decode('ascii', 'replace').strip())
        n_records = int(head[236:244].decode('ascii', 'replace').strip())
        rec_dur = float(head[244:252].decode('ascii', 'replace').strip())
        n_ch = int(head[252:256].decode('ascii', 'replace').strip())
        block = f.read(n_ch * 256)
    dd, mm, yy = (int(x) for x in date_str.split('.'))
    hh, mi, ss = (int(x) for x in time_str.split('.'))
    year = 2000 + yy if yy <= 84 else 1900 + yy          # EDF 2-digit-year rule
    start_dt = datetime.datetime(year, mm, dd, hh, mi, ss)

    def field(offset, width):
        return [block[offset + i * width: offset + (i + 1) * width].decode('latin-1').strip()
                for i in range(n_ch)]
    labels = field(0, 16)
    dims = field(16 * n_ch + 80 * n_ch, 8)
    spr_off = 16 * n_ch + 80 * n_ch + 8 * n_ch + 8 * n_ch * 4 + 80 * n_ch
    samples_per_record = [int(float(x)) if x else 0 for x in field(spr_off, 8)]
    if n_records < 0:                                     # unknown length: infer from the size
        bytes_per_record = 2 * sum(samples_per_record)
        n_records = (edf_path.stat().st_size - header_bytes) // max(bytes_per_record, 1)
    duration_s = n_records * rec_dur
    sfreq = {lab: (spr / rec_dur if rec_dur > 0 else np.nan)
             for lab, spr in zip(labels, samples_per_record)}
    return {'start_dt': start_dt, 'duration_s': float(duration_s), 'n_data_records': int(n_records),
            'record_duration_s': rec_dur, 'n_channels': n_ch, 'ch_names': labels,
            'samples_per_record': samples_per_record, 'sfreq_per_ch': sfreq, 'dims': dims}


def find_channel(ch_names, pattern):
    """First channel label matching a case-insensitive regex, or None."""
    rx = re.compile(pattern, re.IGNORECASE)
    for ch in ch_names:
        if rx.search(ch):
            return ch
    return None


LIGHT_CH_PATTERN = r'^\s*lights?\s*$|lumin'
SPO2_CH_PATTERN = r'spo2|sao2'


# ---------------------------------------------------------------------------
# Hypnogram
# ---------------------------------------------------------------------------
def load_hypnogram(path):
    """One label per line (tool 3 output) -> array of str."""
    return np.loadtxt(str(path), dtype=str).astype('<U10')


def reconcile_hypno_length(stages, duration_s, epoch_s):
    """Compare the hypnogram with the EDF duration. One extra scored epoch (the recording ends
    mid-epoch) is dropped; one missing epoch is kept as is; anything else is a mismatch.
    Returns (stages, check) with check['status'] in {'ok', 'warning', 'fail'}."""
    expected = int(np.floor(duration_s / epoch_s))
    n = len(stages)
    if n == expected:
        return stages, {'check': 'hypno_vs_edf_length', 'status': 'ok',
                        'detail': f'{n} epochs in both'}
    if n - expected == 1:
        return stages[:expected], {'check': 'hypno_vs_edf_length', 'status': 'warning',
                                   'detail': f'hypnogram has 1 epoch more than the EDF ({n} vs '
                                             f'{expected}) - last epoch dropped'}
    if expected - n == 1:
        return stages, {'check': 'hypno_vs_edf_length', 'status': 'warning',
                        'detail': f'hypnogram has 1 epoch less than the EDF ({n} vs {expected}) - kept'}
    return stages, {'check': 'hypno_vs_edf_length', 'status': 'fail',
                    'detail': f'hypnogram has {n} epochs, the EDF {expected} - mismatch > 1 epoch'}


def classify_stages(stages, custom_stages):
    """Boolean masks over the hypnogram: sleep (N1/N2/N3/R), wake (W) and "other" (MT, custom
    stages, unknown labels). Under OTHER_STAGE_POLICY = 'wake' the other epochs join the wake
    mask instead (YASA-style exclusion is the default)."""
    stages = np.asarray(stages, dtype=str)
    sleep = np.isin(stages, SLEEP_STAGES)
    wake = stages == 'W'
    other = ~(sleep | wake)
    known = set(AASM_STAGES) | {'MT'} | set(custom_stages)
    unknown = sorted(str(x) for x in set(stages[other]) - known)
    if OTHER_STAGE_POLICY == 'wake':
        wake = wake | other
        other = np.zeros_like(other)
    return {'sleep_mask': sleep, 'wake_mask': wake, 'other_mask': other, 'unknown_labels': unknown}


# ---------------------------------------------------------------------------
# Lights-off / lights-on
# ---------------------------------------------------------------------------
def load_lights_from_summary_txt(path, col='Luminosité'):
    """Read the per-epoch light state from a Profusion *_Summary_Export.txt (comma-separated,
    French header, UTF-16 with BOM or UTF-8). Column `col` (accent/case-insensitive match):
    1 = lights ON, 0 = lights OFF. Returns a boolean array (True = ON). Raises when the column
    is missing (the caller falls through to the next source)."""
    raw = open(path, 'rb').read()
    if raw[:2] in (b'\xff\xfe', b'\xfe\xff'):
        text = raw.decode('utf-16')
    else:
        text = raw.decode('utf-8-sig', errors='replace')
    df = pd.read_csv(io.StringIO(text), sep=',')
    target = strip_accents(col).strip().lower()
    match = [c for c in df.columns if strip_accents(c).strip().lower() == target]
    if not match:
        raise ValueError(f'column "{col}" not found (columns: {", ".join(map(str, df.columns))})')
    values = pd.to_numeric(df[match[0]], errors='coerce').fillna(0).to_numpy()
    return values == 1


def detect_lights_from_channel(edf_path, ch_name, epoch_s):
    """Per-epoch light state from an EDF step channel (Compumedics 'Light': 1 = ON, 0 = OFF).
    Only that channel is read (`include=` at read time). Returns a boolean array (True = ON)."""
    raw = mne.io.read_raw_edf(str(edf_path), include=[ch_name], preload=True, verbose=False)
    x = raw.get_data()[0]
    sf = raw.info['sfreq']
    n_per_epoch = int(round(epoch_s * sf))
    n = len(x) // n_per_epoch
    if n == 0:
        raise ValueError('Light channel shorter than one epoch')
    med = np.median(x[:n * n_per_epoch].reshape(n, n_per_epoch), axis=1)
    return med > 0.5


def align_lights_to_hypno(on, n_hyp):
    """Bring a per-epoch light array to the hypnogram length: truncate a longer one, pad a
    one-epoch-short one with its last value, reject a larger mismatch (returns None)."""
    on = np.asarray(on, dtype=bool)
    n = len(on)
    if n == n_hyp:
        return on, None
    if n > n_hyp:
        return on[:n_hyp], {'check': 'lights_txt_vs_hypno_length', 'status': 'warning',
                            'detail': f'light file has {n} epochs, hypnogram {n_hyp} - truncated'}
    if n_hyp - n == 1:
        return np.concatenate([on, on[-1:]]), {
            'check': 'lights_txt_vs_hypno_length', 'status': 'warning',
            'detail': f'light file has {n} epochs, hypnogram {n_hyp} - last value repeated'}
    return None, {'check': 'lights_txt_vs_hypno_length', 'status': 'fail',
                  'detail': f'light file has {n} epochs, hypnogram {n_hyp} - source rejected'}


def lights_from_on_flags(on, sleep_mask):
    """Lights-off / lights-on epoch indices from a per-epoch ON/OFF array.

    lights-off = the first epoch AFTER the last ON epoch that precedes the first sleep epoch
    (a bare "first OFF epoch" rule fails when the recording starts with the light already off
    and the technician switches it on briefly before the night). lights-on = the first ON epoch
    after the last sleep epoch (exclusive index). Without an ON epoch on a side the recording
    start / end is used, with a warning. ON epochs inside the sleep period are reported, not used.
    Returns (off_idx, on_idx, checks)."""
    on = np.asarray(on, dtype=bool)
    sleep = np.asarray(sleep_mask, dtype=bool)
    n = len(on)
    checks = []
    if not sleep.any():
        before = np.where(on)[0]
        off_idx = int(before[-1] + 1) if len(before) and before[-1] + 1 < n else 0
        checks.append({'check': 'no_lights_on', 'status': 'warning',
                       'detail': 'no sleep epoch - lights-on set to the recording end'})
        return off_idx, n, checks
    first_sleep = int(np.argmax(sleep))
    last_sleep = int(n - 1 - np.argmax(sleep[::-1]))
    before = np.where(on[:first_sleep])[0]
    if len(before):
        off_idx = int(before[-1] + 1)
    else:
        off_idx = 0
        checks.append({'check': 'no_lights_off', 'status': 'warning',
                       'detail': 'no lights-ON epoch before the first sleep epoch - lights-off '
                                 'set to the recording start'})
    after = np.where(on[last_sleep + 1:])[0]
    if len(after):
        on_idx = int(last_sleep + 1 + after[0])
    else:
        on_idx = n
        checks.append({'check': 'no_lights_on', 'status': 'warning',
                       'detail': 'no lights-ON epoch after the last sleep epoch - lights-on set '
                                 'to the recording end'})
    n_inside = int(on[first_sleep:last_sleep + 1].sum())
    if n_inside:
        checks.append({'check': 'lights_on_inside_sleep', 'status': 'warning',
                       'detail': f'{n_inside} lights-ON epoch(s) between the first and last sleep '
                                 f'epoch - ignored (bathroom break, technician?)'})
    return off_idx, on_idx, checks


def _parse_lights_cell(value, start_dt, epoch_s):
    """One participant-table cell -> epoch index. Accepts 'HH:MM', 'HH:MM:SS' (clock time on the
    recording date, rolled to the next day when it precedes the start by more than 1 h) or a
    bare integer (epoch index). Returns None for a blank/unparseable cell."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    text = str(value).strip()
    if not text or text.lower() in ('nan', 'none', 'nat'):
        return None
    if re.fullmatch(r'\d+', text):
        return int(text)
    m = re.fullmatch(r'(\d{1,2})[:h.](\d{2})(?::(\d{2}))?', text)
    if m is None:
        try:                                           # a pandas/Excel time object
            t = pd.to_datetime(text).time()
        except Exception:
            return None
    else:
        t = datetime.time(int(m.group(1)), int(m.group(2)), int(m.group(3) or 0))
    dt = datetime.datetime.combine(start_dt.date(), t)
    if (start_dt - dt).total_seconds() > 3600:         # after midnight -> next calendar day
        dt += datetime.timedelta(days=1)
    return int(round((dt - start_dt).total_seconds() / epoch_s))


def lights_from_participant_table(row, start_dt, epoch_s, n_hyp, off_col, on_col):
    """Lights-off / lights-on from two columns of the participant table (clock time or epoch
    index). A blank cell leaves that side to the recording bound. Returns (off, on, checks) or
    (None, None, checks) when neither cell is usable."""
    checks = []
    off = _parse_lights_cell(row.get(off_col) if row is not None else None, start_dt, epoch_s)
    on = _parse_lights_cell(row.get(on_col) if row is not None else None, start_dt, epoch_s)
    if off is None and on is None:
        return None, None, checks
    if off is None:
        off = 0
        checks.append({'check': 'no_lights_off', 'status': 'warning',
                       'detail': f'participant table: no "{off_col}" value - recording start used'})
    if on is None:
        on = n_hyp
        checks.append({'check': 'no_lights_on', 'status': 'warning',
                       'detail': f'participant table: no "{on_col}" value - recording end used'})
    if off < 0 or off > n_hyp or on < 0 or on > n_hyp:
        checks.append({'check': 'lights_order', 'status': 'warning',
                       'detail': f'participant table lights outside the recording (off={off}, '
                                 f'on={on}, n={n_hyp}) - clamped'})
        off, on = int(np.clip(off, 0, n_hyp)), int(np.clip(on, 0, n_hyp))
    return off, on, checks


def resolve_lights(cfg, p, header, n_hyp, sleep_mask, subj_row):
    """Try the lights sources in the configured order and return
    {'off_idx', 'on_idx', 'source', 'light_channel', 'checks'}.

    Sources: 'summary_txt' (Profusion export next to the EDF), 'light_channel' (EDF step
    channel), 'participant_table' (two columns), 'recording_bounds' (fallback, always warned).
    'auto' walks them in that order; a forced source that is unavailable falls back to the
    recording bounds with a warning naming what was missing."""
    order = (['summary_txt', 'light_channel', 'participant_table'] if cfg['lights_source'] == 'auto'
             else [cfg['lights_source']])
    checks = []
    light_ch = p.get('light_ch')
    for src in order:
        try:
            if src == 'summary_txt':
                path = p.get('summary_txt')
                if path is None:
                    checks.append({'check': 'lights_source', 'status': 'info',
                                   'detail': 'no summary txt next to the EDF'})
                    continue
                on = load_lights_from_summary_txt(path, cfg['lights_col'])
                on, chk = align_lights_to_hypno(on, n_hyp)
                if chk is not None:
                    checks.append(chk)
                if on is None:
                    continue
                off_idx, on_idx, c2 = lights_from_on_flags(on, sleep_mask)
                return {'off_idx': off_idx, 'on_idx': on_idx, 'source': 'summary_txt',
                        'light_channel': '', 'checks': checks + c2}
            if src == 'light_channel':
                if not light_ch:
                    checks.append({'check': 'lights_source', 'status': 'info',
                                   'detail': 'no Light channel in the EDF'})
                    continue
                on = detect_lights_from_channel(p['edf'], light_ch, cfg['epoch_s'])
                on, chk = align_lights_to_hypno(on, n_hyp)
                if chk is not None:
                    checks.append(chk)
                if on is None:
                    continue
                off_idx, on_idx, c2 = lights_from_on_flags(on, sleep_mask)
                return {'off_idx': off_idx, 'on_idx': on_idx, 'source': 'light_channel',
                        'light_channel': light_ch, 'checks': checks + c2}
            if src == 'participant_table':
                off_idx, on_idx, c2 = lights_from_participant_table(
                    subj_row, header['start_dt'], cfg['epoch_s'], n_hyp,
                    cfg['lights_off_col'], cfg['lights_on_col'])
                if off_idx is None:
                    checks.append({'check': 'lights_source', 'status': 'info',
                                   'detail': 'no lights columns / values in the participant table'})
                    continue
                return {'off_idx': off_idx, 'on_idx': on_idx, 'source': 'participant_table',
                        'light_channel': '', 'checks': checks + c2}
            if src == 'recording_bounds':
                break
        except Exception as exc:
            checks.append({'check': 'lights_source', 'status': 'warning',
                           'detail': f'{src} unusable ({exc}) - next source tried'})
    tried = ', '.join(order) if cfg['lights_source'] == 'auto' else cfg['lights_source']
    if cfg['lights_source'] != 'recording_bounds':
        checks.append({'check': 'lights_fallback', 'status': 'warning',
                       'detail': f'no usable lights source ({tried}) - the recording start and end '
                                 f'are used as lights-off / lights-on'})
    else:
        checks.append({'check': 'lights_fallback', 'status': 'info',
                       'detail': 'recording start / end used as lights-off / lights-on (as requested)'})
    return {'off_idx': 0, 'on_idx': n_hyp, 'source': 'recording_bounds', 'light_channel': '',
            'checks': checks}


# ---------------------------------------------------------------------------
# Macrostructure
# ---------------------------------------------------------------------------
ONSET_RULES = ['first non-W epoch', 'first N2 epoch', 'first of 3 consecutive sleep epochs']


def find_sleep_onset(seg, sleep_seg, rule):
    """Index (within the in-bed segment) of the sleep onset under `rule`, or None."""
    seg = np.asarray(seg, dtype=str)
    sleep_seg = np.asarray(sleep_seg, dtype=bool)
    if rule == 'first N2 epoch':
        hits = np.where(seg == 'N2')[0]
        return int(hits[0]) if len(hits) else None
    if rule == 'first of 3 consecutive sleep epochs':
        for i in range(len(sleep_seg) - 2):
            if sleep_seg[i] and sleep_seg[i + 1] and sleep_seg[i + 2]:
                return i
        return None
    hits = np.where(sleep_seg)[0]                       # default: first non-W (= first sleep) epoch
    return int(hits[0]) if len(hits) else None


def compute_macrostructure(stages, cls, off_idx, on_idx, epoch_s, onset_rule):
    """Every hypnogram-based metric on the epochs [off_idx, on_idx).

    Indices are 0-based, lights-on exclusive. Stage durations are counted WITHIN the sleep
    period (identical to the whole-night count under the default onset rule; under the N2 /
    3-consecutive rules the sleep epochs before onset are left out and reported).
    No sleep at all -> TIB only, everything else NaN.
    Returns (metrics dict, checks list, extras dict)."""
    ep = epoch_s / 60.0
    stages = np.asarray(stages, dtype=str)
    seg = stages[off_idx:on_idx]
    sleep = cls['sleep_mask'][off_idx:on_idx]
    wake = cls['wake_mask'][off_idx:on_idx]
    other = cls['other_mask'][off_idx:on_idx]
    n = len(seg)
    metrics = {k: np.nan for k, _, _, g, _ in METRICS if g == 'macrostructure'}
    metrics['tib_min'] = n * ep
    checks = []
    extras = {'onset_idx_abs': None, 'last_sleep_idx_abs': None, 'n_epochs_in_bed': n,
              'n_other_in_bed': int(other.sum()), 'n_sleep_before_onset': 0}
    onset = find_sleep_onset(seg, sleep, onset_rule)
    if onset is None or not sleep.any():
        checks.append({'check': 'no_sleep', 'status': 'fail',
                       'detail': 'no sleep onset found between lights-off and lights-on - every '
                                 'sleep metric is NaN'})
        return metrics, checks, extras
    last = int(n - 1 - np.argmax(sleep[::-1]))
    spt_seg, spt_sleep, spt_wake, spt_other = (seg[onset:last + 1], sleep[onset:last + 1],
                                               wake[onset:last + 1], other[onset:last + 1])
    metrics['sol_min'] = onset * ep
    metrics['spt_min'] = (last - onset + 1) * ep
    metrics['tst_min'] = int(spt_sleep.sum()) * ep
    metrics['waso_min'] = int(spt_wake.sum()) * ep
    metrics['other_min'] = int(spt_other.sum()) * ep
    metrics['se_pct'] = 100.0 * metrics['tst_min'] / metrics['tib_min'] if metrics['tib_min'] > 0 else np.nan
    metrics['sme_pct'] = 100.0 * metrics['tst_min'] / metrics['spt_min'] if metrics['spt_min'] > 0 else np.nan
    for stage, key in [('N1', 'n1'), ('N2', 'n2'), ('N3', 'n3'), ('R', 'rem')]:
        hits = np.where(spt_seg == stage)[0]
        metrics[f'lat_{key}_min'] = hits[0] * ep if len(hits) else np.nan
        metrics[f'{key}_min'] = len(hits) * ep
        metrics[f'{key}_pct'] = (100.0 * metrics[f'{key}_min'] / metrics['tst_min']
                                 if metrics['tst_min'] > 0 else np.nan)
    extras['onset_idx_abs'] = off_idx + onset
    extras['last_sleep_idx_abs'] = off_idx + last
    extras['n_sleep_before_onset'] = int(sleep[:onset].sum())
    if extras['n_sleep_before_onset']:
        checks.append({'check': 'sleep_before_onset', 'status': 'info',
                       'detail': f"{extras['n_sleep_before_onset']} sleep epoch(s) before the "
                                 f'"{onset_rule}" onset are not counted in TST / stage durations'})
    return metrics, checks, extras


# ---------------------------------------------------------------------------
# Event indices (AHI family, PLM, arousals) and ODI
# ---------------------------------------------------------------------------
def canonical_event_class(canonical):
    """Map a canonical label (tool 4 vocabulary) to the class counted here, or None.
    `hypopnea*` and `arousal*` are matched by prefix: users type variants such as
    `hypopnea_obstructive` / `arousal_respiratory` and they all belong to the same index."""
    if canonical is None:
        return None
    c = str(canonical).strip().lower()
    if c in ('apnea_obstructive', 'apnea_central', 'apnea_mixed', 'plm', 'spo2_desaturation'):
        return c
    if c.startswith('hypopnea'):
        return 'hypopnea'
    if c.startswith('arousal'):
        return 'arousal'
    return None


def select_events_in_window(events_df, cls, off_idx, on_idx, epoch_s, mode='tib'):
    """Keep the events whose onset falls in the counting window.
    'tib': onset in [lights-off, lights-on); 'sleep_epochs': onset in a sleep epoch.
    Returns (df, n_outside_window, n_bad_onset)."""
    if events_df is None or not len(events_df):
        return events_df, 0, 0
    start = pd.to_numeric(events_df['Start'], errors='coerce').to_numpy(dtype=float)
    bad = ~np.isfinite(start)
    ep_idx = np.floor(np.where(bad, -1, start) / epoch_s).astype(int)
    inside = (~bad) & (ep_idx >= off_idx) & (ep_idx < on_idx)
    if mode == 'sleep_epochs':
        sleep = np.asarray(cls['sleep_mask'], dtype=bool)
        ok_idx = (ep_idx >= 0) & (ep_idx < len(sleep))
        in_sleep = np.zeros_like(inside)
        in_sleep[ok_idx] = sleep[ep_idx[ok_idx]]
        inside = inside & in_sleep
    n_outside = int(((~inside) & (~bad)).sum())
    return events_df[inside].copy(), n_outside, int(bad.sum())


def compute_event_indices(events_df, remap, tst_min):
    """Counts and per-hour indices of the scored events (already windowed).
    Returns (metrics dict, checks list, n_unmapped). A TST of 0 / NaN gives NaN indices."""
    keys = [k for k, _, _, g, _ in METRICS if g in ('respiratory', 'plm', 'arousal')]
    metrics = {k: np.nan for k in keys}
    checks = []
    if events_df is None:
        return metrics, checks, 0
    counts = {'apnea_obstructive': 0, 'apnea_central': 0, 'apnea_mixed': 0, 'hypopnea': 0,
              'plm': 0, 'arousal': 0}
    unmapped = {}
    for name in events_df['Name'].astype(str):
        canonical = (remap or {}).get(name)
        if canonical is None:
            if remap is None or name not in remap:      # absent from the remap (never seen by tool 4)
                unmapped[name] = unmapped.get(name, 0) + 1
            continue                                    # explicitly ignored (null) or unmapped
        klass = canonical_event_class(canonical)
        if klass in counts:
            counts[klass] += 1
    metrics.update({'n_apnea_obstructive': counts['apnea_obstructive'],
                    'n_apnea_central': counts['apnea_central'], 'n_apnea_mixed': counts['apnea_mixed'],
                    'n_hypopnea': counts['hypopnea'], 'n_plm': counts['plm'],
                    'n_arousal': counts['arousal']})
    tst_h = tst_min / 60.0 if (tst_min is not None and np.isfinite(tst_min) and tst_min > 0) else np.nan
    metrics['oai'] = counts['apnea_obstructive'] / tst_h
    metrics['cai'] = counts['apnea_central'] / tst_h
    metrics['mai'] = counts['apnea_mixed'] / tst_h
    metrics['hi'] = counts['hypopnea'] / tst_h
    metrics['ahi'] = metrics['oai'] + metrics['cai'] + metrics['mai'] + metrics['hi']
    metrics['plm_index'] = counts['plm'] / tst_h
    metrics['arousal_index'] = counts['arousal'] / tst_h
    n_unmapped = int(sum(unmapped.values()))
    if unmapped:
        listing = ', '.join(f'{k} x{v}' for k, v in sorted(unmapped.items(), key=lambda kv: -kv[1]))
        checks.append({'check': 'unmapped_events', 'status': 'warning',
                       'detail': f'{n_unmapped} event(s) with a label absent from event_remap.json '
                                 f'(not counted; map them with tool 4): {listing}'})
    return metrics, checks, n_unmapped


def compute_odi(events_df, remap, tst_min, threshold_pct, xml_path, cls, off_idx, on_idx, epoch_s,
                window_mode):
    """Oxygen desaturation index from the scored desaturations.
    When the events carry a `Desaturation` depth (XML source, or the XML companion re-read for
    that purpose) only the events with depth >= threshold are counted; otherwise every scored
    desaturation is counted and the rule says so. Returns (metrics dict, checks, rule)."""
    metrics = {'n_desaturation': np.nan, 'odi': np.nan}
    checks = []
    df = events_df
    rule = 'all_events'
    has_depth = df is not None and 'Desaturation' in df.columns and df['Desaturation'].notna().any()
    if not has_depth and xml_path is not None:
        try:
            xdf = _events_df_from_xml(xml_path)
            if xdf is not None and xdf['Desaturation'].notna().any():
                xdf, _, _ = select_events_in_window(xdf, cls, off_idx, on_idx, epoch_s, window_mode)
                df, has_depth = xdf, True
                checks.append({'check': 'odi_rule', 'status': 'info',
                               'detail': 'desaturation depths read from the .edf.XML companion'})
        except Exception as exc:
            checks.append({'check': 'odi_rule', 'status': 'warning',
                           'detail': f'could not read the XML desaturation depths ({exc})'})
    if df is None:
        return metrics, checks, rule
    is_desat = np.array([canonical_event_class((remap or {}).get(str(n))) == 'spo2_desaturation'
                         for n in df['Name']], dtype=bool)
    n_all = int(is_desat.sum())
    if has_depth:
        depth = pd.to_numeric(df['Desaturation'], errors='coerce').to_numpy(dtype=float)
        deep = is_desat & np.isfinite(depth) & (depth >= threshold_pct)
        # depth unknown for a desaturation event -> counted (benefit of the doubt), reported
        no_depth = is_desat & ~np.isfinite(depth)
        n = int(deep.sum() + no_depth.sum())
        rule = f'depth>={threshold_pct:g}%'
        if no_depth.any():
            checks.append({'check': 'odi_rule', 'status': 'warning',
                           'detail': f'{int(no_depth.sum())} desaturation(s) without a depth value '
                                     f'were counted regardless of the threshold'})
    else:
        n = n_all
        if n_all:
            checks.append({'check': 'odi_rule', 'status': 'warning',
                           'detail': 'no desaturation depth available (txt/csv events, no XML) - '
                                     'every scored desaturation counted, threshold not applied'})
    metrics['n_desaturation'] = n
    tst_h = tst_min / 60.0 if (tst_min is not None and np.isfinite(tst_min) and tst_min > 0) else np.nan
    metrics['odi'] = n / tst_h
    return metrics, checks, rule


# ---------------------------------------------------------------------------
# SpO2 signal -> T90
# ---------------------------------------------------------------------------
def load_spo2_signal(edf_path, ch_name):
    """Read ONE channel (`include=` at read time) -> (values in %, sfreq). A channel stored as a
    fraction (max <= 1) is rescaled to %, which the caller reports."""
    raw = mne.io.read_raw_edf(str(edf_path), include=[ch_name], preload=True, verbose=False)
    x = raw.get_data()[0].astype(float)
    sf = float(raw.info['sfreq'])
    rescaled = False
    finite = x[np.isfinite(x)]
    if len(finite) and np.nanmax(finite) <= 1.0:
        x = x * 100.0
        rescaled = True
    return x, sf, rescaled


def compute_t90(x, sfreq, off_idx, on_idx, epoch_s, tst_min, thr=90.0, floor=50.0, sleep_mask=None):
    """Time below `thr` % between lights-off and lights-on, on the valid samples only
    (a sample is artefact when <= 0 or below `floor` - sensor off, motion). Also T88 (time
    <= 88 % on the sleep epochs of `sleep_mask`), the mean, the nadir and the artefact
    percentage. Returns (metrics dict, checks)."""
    metrics = {'t90_min': np.nan, 't90_pct_tst': np.nan, 't88_sleep_min': np.nan, 'spo2_mean_pct': np.nan,
               'spo2_nadir_pct': np.nan, 'spo2_artifact_pct': np.nan}
    checks = []
    i0 = int(round(off_idx * epoch_s * sfreq))
    i1 = int(round(on_idx * epoch_s * sfreq))
    if i1 > len(x):
        checks.append({'check': 'spo2_length', 'status': 'warning',
                       'detail': f'SpO2 signal ends {(i1 - len(x)) / sfreq:.0f} s before lights-on - '
                                 f'computed on the available samples'})
        i1 = len(x)
    seg = x[i0:i1]
    if len(seg) == 0:
        checks.append({'check': 'spo2_length', 'status': 'fail', 'detail': 'no SpO2 sample in bed'})
        return metrics, checks
    valid = np.isfinite(seg) & (seg > 0) & (seg >= floor)
    n_valid = int(valid.sum())
    metrics['spo2_artifact_pct'] = 100.0 * (1 - n_valid / len(seg))
    if n_valid == 0:
        checks.append({'check': 'spo2_artifact', 'status': 'fail',
                       'detail': 'every SpO2 sample in bed is artefact (<= 0 or below the floor)'})
        return metrics, checks
    below = valid & (seg < thr)
    metrics['t90_min'] = float(below.sum()) / sfreq / 60.0
    metrics['t90_pct_tst'] = (100.0 * metrics['t90_min'] / tst_min
                              if (tst_min is not None and np.isfinite(tst_min) and tst_min > 0) else np.nan)
    # ICSD-3-TR sleep-related hypoxemia, criterion A (adults): SpO2 <= 88 % for >= 5 min DURING
    # SLEEP. The threshold is fixed (not the T90 widget) and only the sleep epochs count, so each
    # sample is mapped to its epoch. No sleep -> NaN, like every other sleep metric.
    if sleep_mask is not None and tst_min is not None and np.isfinite(tst_min) and tst_min > 0:
        epoch_of_sample = (np.arange(i0, i0 + len(seg)) / (epoch_s * sfreq)).astype(int)
        epoch_of_sample = np.clip(epoch_of_sample, 0, len(sleep_mask) - 1)
        asleep = np.asarray(sleep_mask, dtype=bool)[epoch_of_sample]
        low_in_sleep = valid & asleep & (seg <= 88.0)
        metrics['t88_sleep_min'] = float(low_in_sleep.sum()) / sfreq / 60.0
        if metrics['t88_sleep_min'] >= 5.0:
            checks.append({'check': 'icsd_hypoxemia', 'status': 'info',
                           'detail': f"ICSD-3-TR sleep-related hypoxemia criterion A met "
                                     f"({metrics['t88_sleep_min']:.1f} min with SpO2 <= 88 % during sleep) - "
                                     f'criterion B (not explained by hypoventilation / OSA) needs a '
                                     f'clinical evaluation'})
    metrics['spo2_mean_pct'] = float(np.mean(seg[valid]))
    metrics['spo2_nadir_pct'] = float(np.min(seg[valid]))
    if metrics['spo2_artifact_pct'] > 10:
        checks.append({'check': 'spo2_artifact', 'status': 'warning',
                       'detail': f"{metrics['spo2_artifact_pct']:.1f} % of the SpO2 samples in bed "
                                 f'discarded as artefact - T90 may be underestimated'})
    return metrics, checks


# ---------------------------------------------------------------------------
# Consistency checks and output tables
# ---------------------------------------------------------------------------
def run_checks(metrics, extras, lights, n_hyp, cls, event_info):
    """Arithmetic and plausibility checks on one participant's metrics.
    Returns a list of {'check', 'status', 'detail'} (status: ok / info / warning / fail)."""
    out = []

    def add(check, ok, detail, bad='fail'):
        out.append({'check': check, 'status': 'ok' if ok else bad, 'detail': detail})

    m = metrics
    has_sleep = np.isfinite(m.get('tst_min', np.nan))
    if has_sleep:
        total = m['tst_min'] + m['waso_min'] + m['other_min']
        add('spt_consistency', abs(m['spt_min'] - total) < 1e-6,
            f"SPT {m['spt_min']:.1f} = TST {m['tst_min']:.1f} + WASO {m['waso_min']:.1f}"
            + (f" + other {m['other_min']:.1f}" if m['other_min'] > 0 else ''))
        add('se_range', 0 <= m['se_pct'] <= 100 + 1e-9, f"SE = {m['se_pct']:.1f} %")
        add('tst_le_spt_le_tib', m['tst_min'] <= m['spt_min'] + 1e-9 <= m['tib_min'] + 1e-9,
            f"TST {m['tst_min']:.1f} <= SPT {m['spt_min']:.1f} <= TIB {m['tib_min']:.1f}")
        add('sol_plus_spt_le_tib', m['sol_min'] + m['spt_min'] <= m['tib_min'] + 1e-9,
            f"SOL {m['sol_min']:.1f} + SPT {m['spt_min']:.1f} <= TIB {m['tib_min']:.1f}")
        pct_sum = m['n1_pct'] + m['n2_pct'] + m['n3_pct'] + m['rem_pct']
        add('stage_pct_sum', abs(pct_sum - 100) < 0.1, f'N1+N2+N3+REM = {pct_sum:.2f} %')
        lats = {k: m[k] for k in ('lat_n1_min', 'lat_n2_min', 'lat_n3_min', 'lat_rem_min')}
        bad = [k for k, v in lats.items() if np.isfinite(v) and not (0 <= v <= m['spt_min'])]
        add('latencies_range', not bad,
            'every latency within [0, SPT]' if not bad else f'out of range: {", ".join(bad)}')
        missing = [k for k, v in lats.items() if not np.isfinite(v)]
        if missing:
            out.append({'check': 'stage_absent', 'status': 'info',
                        'detail': 'stage never reached (latency NaN): '
                                  + ', '.join(k.replace('lat_', '').replace('_min', '').upper()
                                              for k in missing)})
    add('lights_order', 0 <= lights['off_idx'] < lights['on_idx'] <= n_hyp,
        f"lights-off epoch {lights['off_idx']} < lights-on epoch {lights['on_idx']} <= {n_hyp}")
    n_other = int(cls['other_mask'].sum())
    if n_other or cls['unknown_labels']:
        labels = ', '.join(cls['unknown_labels']) if cls['unknown_labels'] else 'MT / custom'
        policy = ('counted in TIB and SPT, excluded from TST and WASO (YASA convention)'
                  if OTHER_STAGE_POLICY == 'exclude' else 'counted as wake (WASO inside SPT)')
        out.append({'check': 'mt_custom_present', 'status': 'warning',
                    'detail': f'{n_other} epoch(s) neither sleep nor wake ({labels}): {policy}. '
                              f'Re-map them with tool 3 if they should be scored as a stage.'})
    if event_info.get('source') is None and event_info.get('enabled'):
        out.append({'check': 'events_source', 'status': 'info',
                    'detail': 'no scored-event companion found - event indices are NaN'})
    if event_info.get('n_bad_onset'):
        out.append({'check': 'event_onsets_in_recording', 'status': 'warning',
                    'detail': f"{event_info['n_bad_onset']} event(s) with a missing / unreadable onset "
                              f'were skipped'})
    if event_info.get('n_outside'):
        out.append({'check': 'events_outside_window', 'status': 'info',
                    'detail': f"{event_info['n_outside']} event(s) outside the counting window "
                              f'({EVENT_COUNT_WINDOW}) were not counted'})
    return out


def metrics_long_df(file_id, metrics, provenance):
    """Long table: one row per registered metric (NaN when not computed) + provenance rows."""
    rows = []
    for key, label, unit, group, definition in METRICS:
        val = metrics.get(key, np.nan)
        rows.append({'file_id': file_id, 'group': group, 'metric': key,
                     'value': (round(float(val), 4) if val is not None and np.isfinite(float(val))
                               else np.nan),
                     'unit': unit, 'definition_short': definition})
    for key in PROVENANCE_KEYS:
        rows.append({'file_id': file_id, 'group': 'provenance', 'metric': key,
                     'value': provenance.get(key, ''), 'unit': '', 'definition_short': ''})
    return pd.DataFrame(rows, columns=['file_id', 'group', 'metric', 'value', 'unit', 'definition_short'])


def checks_df(file_id, checks):
    return pd.DataFrame([{'file_id': file_id, **c} for c in checks],
                        columns=['file_id', 'check', 'status', 'detail'])


def long_to_wide(long_df):
    """Pivot the concatenated long tables back to one row per participant: metric columns as
    float (METRICS order), provenance columns as text, in a fixed order."""
    if long_df is None or not len(long_df):
        return pd.DataFrame(columns=['file_id'] + METRIC_KEYS + PROVENANCE_KEYS)
    rows = []
    for fid, g in long_df.groupby('file_id', sort=True):
        vals = dict(zip(g['metric'], g['value']))
        row = {'file_id': fid}
        for k in METRIC_KEYS:
            row[k] = pd.to_numeric(pd.Series([vals.get(k)]), errors='coerce').iloc[0]
        for k in PROVENANCE_KEYS:
            v = vals.get(k, '')
            row[k] = '' if (v is None or (isinstance(v, float) and np.isnan(v))) else v
        rows.append(row)
    return pd.DataFrame(rows, columns=['file_id'] + METRIC_KEYS + PROVENANCE_KEYS)


# ---------------------------------------------------------------------------
# Comparison with an external reference table (Section 5)
# ---------------------------------------------------------------------------
# Aliases used to guess which reference column holds which metric (normalised: lower-case,
# accents stripped, non-alphanumerics removed). Extend freely.
METRIC_ALIASES = {
    'tib_min': ['tib', 'timeinbed', 'tempsaulit', 'tempsaulitmin', 'tpsaulit'],
    'sol_min': ['sol', 'sleeponsetlatency', 'sleeplatency', 'latenceendormissement', 'latencedendormissement', 'latencesommeil'],
    'spt_min': ['spt', 'sleepperiodtime', 'periodedesommeil', 'tempsperiodesommeil'],
    'tst_min': ['tst', 'totalsleeptime', 'tempstotaldesommeil', 'tempsdesommeiltotal', 'tts'],
    'waso_min': ['waso', 'wakeaftersleeponset', 'eveilintrasommeil', 'eveilapresendormissement'],
    'se_pct': ['se', 'sepct', 'sleepefficiency', 'sleepefficiencypct', 'efficacitedusommeil',
               'efficacitesommeil', 'efficacitesommeilpct', 'effsommeil', 'effsommeilpct', 'efficacite'],
    'sme_pct': ['sme', 'sleepmaintenanceefficiency'],
    'lat_n1_min': ['n1latency', 'latencyn1', 'latencen1'],
    'lat_n2_min': ['n2latency', 'latencyn2', 'latencen2'],
    'lat_n3_min': ['n3latency', 'latencyn3', 'latencen3', 'swslatency', 'latencesws'],
    'lat_rem_min': ['remlatency', 'latencyrem', 'latencerem', 'latencesp', 'latenceparadoxal'],
    'n1_min': ['n1min', 'n1duration', 'n1'],
    'n2_min': ['n2min', 'n2duration', 'n2'],
    'n3_min': ['n3min', 'n3duration', 'n3', 'sws', 'swsmin'],
    'rem_min': ['remmin', 'remduration', 'rem', 'sp', 'spmin'],
    'n1_pct': ['n1pct', 'n1percent', 'n1tst', 'pctn1', 'n1ptst'],
    'n2_pct': ['n2pct', 'n2percent', 'n2tst', 'pctn2', 'n2ptst'],
    'n3_pct': ['n3pct', 'n3percent', 'n3tst', 'pctn3', 'n3ptst', 'swspct'],
    'rem_pct': ['rempct', 'rempercent', 'remtst', 'pctrem', 'remptst', 'sppct'],
    'ahi': ['ahi', 'iah', 'apneahypopneaindex', 'indexapneeshypopnees'],
    'oai': ['oai', 'obstructiveapneaindex', 'iao'],
    'cai': ['cai', 'centralapneaindex', 'iac'],
    'mai': ['mai', 'mixedapneaindex', 'iam'],
    'hi': ['hi', 'hypopneaindex', 'indexhypopnees', 'ih'],
    'n_apnea_obstructive': ['obstructiveapneas', 'apneesobstructives', 'nobstructiveapnea'],
    'n_apnea_central': ['centralapneas', 'apneescentrales', 'ncentralapnea'],
    'n_apnea_mixed': ['mixedapneas', 'apneesmixtes', 'nmixedapnea'],
    'n_hypopnea': ['hypopneas', 'hypopnees', 'nhypopnea'],
    'plm_index': ['plmindex', 'plmi', 'indexplm', 'mpjindex', 'indexmpj', 'plm'],
    'n_plm': ['nplm', 'plmcount', 'nombreplm'],
    'arousal_index': ['arousalindex', 'ai', 'indexmicroeveils', 'microeveilsindex', 'indexeveils'],
    'n_arousal': ['arousals', 'narousal', 'microeveils', 'nombremicroeveils'],
    'odi': ['odi', 'oxygendesaturationindex', 'indexdesaturation', 'ido', 'odi3', 'odi4'],
    'n_desaturation': ['desaturations', 'ndesaturation', 'nombredesaturations'],
    't90_min': ['t90', 't90min', 'tempsspo290', 'timebelow90'],
    't90_pct_tst': ['t90pct', 't90percent', 'ct90'],
    't88_sleep_min': ['t88', 't88min', 't88sleep', 'timebelow88'],
    'spo2_mean_pct': ['meanspo2', 'spo2mean', 'spo2moyenne', 'saturationmoyenne'],
    'spo2_nadir_pct': ['nadirspo2', 'spo2nadir', 'minspo2', 'spo2min', 'saturationminimale'],
}


def _norm_name(text):
    """'Eff. sommeil %' -> 'effsommeilpct' (a % sign is kept as 'pct' so 'N3 %' cannot match 'N3')."""
    return re.sub(r'[^a-z0-9]', '', strip_accents(str(text).replace('%', ' pct ')).lower())


def guess_metric_mapping(ref_columns):
    """{metric key: reference column or '(none)'}: exact alias match first, then a fuzzy match on
    the normalised names (cutoff 0.8 so 'n1' never steals 'n1pct'); each column used once."""
    cols = [str(c) for c in ref_columns]
    normed = {_norm_name(c): c for c in cols}
    mapping, used = {}, set()
    for key in METRIC_KEYS:                                # exact alias pass
        aliases = [_norm_name(key)] + METRIC_ALIASES.get(key, [])
        hit = next((normed[a] for a in aliases if a in normed and normed[a] not in used), None)
        if hit is not None:
            mapping[key], used = hit, used | {hit}
    for key in METRIC_KEYS:                                # fuzzy pass for what is left
        if key in mapping:
            continue
        aliases = [_norm_name(key)] + METRIC_ALIASES.get(key, [])
        best = None
        for a in aliases:
            cands = difflib.get_close_matches(a, [n for n, c in normed.items() if c not in used],
                                              n=1, cutoff=0.8)
            if cands:
                best = normed[cands[0]]
                break
        mapping[key] = best if best is not None else '(none)'
        if best is not None:
            used.add(best)
    return mapping


def default_tolerance(key):
    """Absolute tolerance by unit: 1 min, 1 %-point, 1 event/h, exact for counts."""
    unit = METRIC_INFO[key]['unit']
    return 0.0 if unit == 'count' else 1.0


def compare_with_reference(wide_df, ref_df, id_col, mapping, tolerances):
    """One row per (participant, mapped metric): our value, the reference value, the difference
    and a status - match / mismatch (|diff| > tol) / missing_ours / missing_reference /
    not_in_reference_table (participant absent from the reference)."""
    ours = wide_df.copy()
    ours['_key'] = ours['file_id'].astype(str).map(norm)
    ref = ref_df.copy()
    ref['_key'] = ref[id_col].astype(str).str.strip().map(norm)
    ref = ref.drop_duplicates('_key')
    ref_keys = set(ref['_key'])
    ref_idx = ref.set_index('_key')
    rows = []
    for _, r in ours.iterrows():
        in_ref = r['_key'] in ref_keys
        for key, col in mapping.items():
            if col in (None, '(none)'):
                continue
            our_v = pd.to_numeric(pd.Series([r.get(key)]), errors='coerce').iloc[0]
            ref_v = (pd.to_numeric(pd.Series([ref_idx.loc[r['_key'], col]]), errors='coerce').iloc[0]
                     if in_ref else np.nan)
            tol = float(tolerances.get(key, default_tolerance(key)))
            if not in_ref:
                status = 'not_in_reference_table'
            elif not np.isfinite(our_v) and not np.isfinite(ref_v):
                status = 'missing_both'
            elif not np.isfinite(our_v):
                status = 'missing_ours'
            elif not np.isfinite(ref_v):
                status = 'missing_reference'
            else:
                status = 'match' if abs(our_v - ref_v) <= tol + 1e-9 else 'mismatch'
            rows.append({'file_id': r['file_id'], 'metric': key, 'ours': our_v, 'reference': ref_v,
                         'ref_column': col, 'diff': (our_v - ref_v) if np.isfinite(our_v) and np.isfinite(ref_v) else np.nan,
                         'tol': tol, 'status': status})
    return pd.DataFrame(rows, columns=['file_id', 'metric', 'ours', 'reference', 'ref_column',
                                       'diff', 'tol', 'status'])


# =============================================================================
# Figures and HTML blocks of the database report
# =============================================================================
STATUS_COLORS = {'ok': '#2e7d32', 'info': '#1565c0', 'warning': '#ef6c00', 'fail': '#c62828'}
SOURCE_COLORS = {'summary_txt': '#1b9e77', 'light_channel': '#7570b3', 'participant_table': '#d95f02',
                 'recording_bounds': '#999999'}


def params_html(cfg):
    """Two-column 'parameter -> value' table."""
    rows = ''.join(f'<tr><td style="padding:2px 12px 2px 0"><b>{k}</b></td>'
                   f'<td>{v}</td></tr>' for k, v in cfg.items())
    return f'<table style="font-size:90%">{rows}</table>'


def save_report_html(path, title, items):
    """Write an mne.Report from an ORDERED list of ('fig'|'html', title, payload) items
    (`None` payloads are skipped)."""
    report = mne.Report(title=title, verbose=False)
    for kind, item_title, payload in items:
        if payload is None:
            continue
        if kind == 'fig':
            report.add_figure(fig=payload, title=item_title, image_format='PNG')
        else:
            report.add_html(html=payload, title=item_title)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    report.save(str(path), overwrite=True, open_browser=False, verbose=False)
    return path


def close_figs(items):
    for kind, _, payload in items:
        if kind == 'fig' and payload is not None:
            plt.close(payload)


def _shade_reference(ax, key, horizontal=False):
    """Shade the indicative normal range of a metric (and the AHI severity bands)."""
    if key == 'ahi' or key == 'odi':
        colors = ['#e8f5e9', '#fff8e1', '#ffe0b2', '#ffcdd2']
        if not horizontal:
            ax.set_ylim(min(ax.get_ylim()[0], 0), max(ax.get_ylim()[1], 6))
        top = ax.get_ylim()[1] if not horizontal else ax.get_xlim()[1]
        for (lo, hi, name), col in zip(AHI_SEVERITY_BANDS, colors):
            hi_v = hi if hi is not None else max(top, lo + 5)
            if horizontal:
                ax.axvspan(lo, hi_v, color=col, zorder=0)
                ax.text((lo + hi_v) / 2, ax.get_ylim()[1], name, ha='center', va='top', fontsize=7, color='#555')
            else:
                ax.axhspan(lo, hi_v, color=col, zorder=0)
                ax.text(ax.get_xlim()[1], (lo + hi_v) / 2, name, ha='right', va='center', fontsize=7, color='#555')
        return
    rng = REFERENCE_RANGES.get(key)
    if rng is None:
        return
    lo, hi, _ = rng
    if horizontal:
        ax.axvspan(lo, hi, color='#e8f5e9', zorder=0)
    else:
        ax.axhspan(lo, hi, color='#e8f5e9', zorder=0)
        y0, y1 = ax.get_ylim()                 # keep the band visible even when every point is outside
        ax.set_ylim(min(y0, lo), max(y1, hi))


def plot_metric_distributions(wide, groups, title):
    """One panel per metric of `groups` with at least one value: every participant as a dot
    (jittered) over a boxplot, the indicative normal range shaded green. None when empty."""
    keys = [k for k in METRIC_KEYS if METRIC_INFO[k]['group'] in groups
            and k in wide.columns and pd.to_numeric(wide[k], errors='coerce').notna().any()]
    if not keys:
        return None
    n_cols = min(4, len(keys))
    n_rows = int(np.ceil(len(keys) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.4 * n_cols, 3.0 * n_rows), squeeze=False)
    rng = np.random.default_rng(0)
    for ax, key in zip(axes.ravel(), keys):
        vals = pd.to_numeric(wide[key], errors='coerce').dropna().to_numpy(dtype=float)
        ax.boxplot(vals, positions=[0], widths=0.5, showfliers=False,
                   medianprops={'color': '#c62828'}, zorder=2)
        ax.scatter(rng.uniform(-0.15, 0.15, len(vals)), vals, s=14, color='#455a64', alpha=0.75, zorder=3)
        ax.set_xlim(-0.6, 0.6)
        ax.set_xticks([])
        ax.set_title(METRIC_INFO[key]['label'], fontsize=9)
        ax.set_ylabel(METRIC_INFO[key]['unit'], fontsize=8)
        ax.tick_params(labelsize=8)
        _shade_reference(ax, key)
        ax.text(0.02, 0.98, f'n = {len(vals)}', transform=ax.transAxes, ha='left', va='top',
                fontsize=7, color='#555')
    for ax in axes.ravel()[len(keys):]:
        ax.axis('off')
    fig.suptitle(f'{title} - green band = indicative normal range (healthy adults)', fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return fig


def _clock_to_hours(text):
    """'23:15:00' -> 23.25; after-midnight times shifted by +24 so the night reads left to right."""
    try:
        hh, mm, ss = (int(x) for x in str(text).split(':'))
    except Exception:
        return np.nan
    h = hh + mm / 60 + ss / 3600
    return h + 24 if h < 12 else h


def plot_lights_overview(wide):
    """Per participant: the in-bed span (lights-off -> lights-on, coloured by the lights source)
    and the sleep period (sleep onset -> last sleep epoch, dark). Clock time on the x axis."""
    if wide is None or not len(wide) or 'lights_off_clock' not in wide.columns:
        return None
    df = wide.copy()
    df['off_h'] = df['lights_off_clock'].map(_clock_to_hours)
    df['on_h'] = df['lights_on_clock'].map(_clock_to_hours)
    df = df[df['off_h'].notna() & df['on_h'].notna()]
    if not len(df):
        return None
    df = df.sort_values('off_h').reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(11, max(2.5, 0.28 * len(df) + 1.2)))
    for i, r in df.iterrows():
        src = r.get('lights_source', 'recording_bounds')
        ax.barh(i, r['on_h'] - r['off_h'], left=r['off_h'], height=0.6,
                color=SOURCE_COLORS.get(src, '#999999'), alpha=0.45, edgecolor='none')
        sol = pd.to_numeric(pd.Series([r.get('sol_min')]), errors='coerce').iloc[0]
        spt = pd.to_numeric(pd.Series([r.get('spt_min')]), errors='coerce').iloc[0]
        if np.isfinite(sol) and np.isfinite(spt):
            ax.barh(i, spt / 60, left=r['off_h'] + sol / 60, height=0.35, color='#37474f')
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['file_id'].astype(str), fontsize=7)
    ax.invert_yaxis()
    lo = np.floor(df['off_h'].min())
    hi = np.ceil(df['on_h'].max())
    ticks = np.arange(lo, hi + 1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{int(t) % 24:02d}:00' for t in ticks], fontsize=8)
    ax.set_xlabel('clock time')
    ax.grid(axis='x', alpha=0.3)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, alpha=0.45) for c in SOURCE_COLORS.values()]
    handles.append(plt.Rectangle((0, 0), 1, 1, color='#37474f'))
    ax.legend(handles, list(SOURCE_COLORS.keys()) + ['sleep period (onset -> last sleep epoch)'],
              fontsize=7, loc='upper left', bbox_to_anchor=(1.01, 1.0), title='lights source',
              title_fontsize=7)
    ax.set_title('Time in bed and sleep period per participant', fontsize=10)
    fig.tight_layout()
    return fig


def plot_comparison(cmp_df):
    """Ours vs reference, one panel per compared metric, identity line, mismatches in red."""
    if cmp_df is None or not len(cmp_df):
        return None
    keys = [k for k in METRIC_KEYS if k in set(cmp_df['metric'])]
    keys = [k for k in keys if cmp_df.loc[cmp_df['metric'] == k, 'status'].isin(['match', 'mismatch']).any()]
    if not keys:
        return None
    n_cols = min(4, len(keys))
    n_rows = int(np.ceil(len(keys) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.4 * n_cols, 3.2 * n_rows), squeeze=False)
    for ax, key in zip(axes.ravel(), keys):
        d = cmp_df[(cmp_df['metric'] == key) & cmp_df['status'].isin(['match', 'mismatch'])]
        col = d['status'].map({'match': '#2e7d32', 'mismatch': '#c62828'})
        ax.scatter(d['reference'], d['ours'], s=16, c=col, alpha=0.8, zorder=3)
        both = np.concatenate([d['reference'].to_numpy(float), d['ours'].to_numpy(float)])
        both = both[np.isfinite(both)]
        if len(both):
            lo, hi = both.min(), both.max()
            pad = 0.05 * (hi - lo if hi > lo else 1)
            ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color='#999', lw=0.8, zorder=1)
        n_mis = int((d['status'] == 'mismatch').sum())
        ax.set_title(f"{METRIC_INFO[key]['label']}\n{len(d) - n_mis} match / {n_mis} mismatch", fontsize=8)
        ax.set_xlabel('reference', fontsize=8)
        ax.set_ylabel('this tool', fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in axes.ravel()[len(keys):]:
        ax.axis('off')
    fig.suptitle('Comparison with the reference table (red = beyond tolerance)', fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def checks_summary_html(checks):
    """Count of participants per (check, status) + the list of participants with a failure."""
    if checks is None or not len(checks):
        return '<p><i>(no checks)</i></p>'
    piv = (checks.groupby(['check', 'status'])['file_id'].nunique().unstack('status')
           .reindex(columns=[s for s in ('ok', 'info', 'warning', 'fail')]).fillna(0).astype(int))
    piv = piv.reset_index()
    html = piv.to_html(index=False, border=0)
    fails = checks[checks['status'] == 'fail']
    if len(fails):
        items = ''.join(f'<li><b>{r.file_id}</b> - {r.check}: {r.detail}</li>' for r in fails.itertuples())
        html += f'<p style="color:#c62828"><b>Failures:</b></p><ul>{items}</ul>'
    warns = checks[checks['status'] == 'warning']
    if len(warns):
        counts = warns.groupby('check')['file_id'].nunique().sort_values(ascending=False)
        html += ('<p><b>Warnings</b> (participants affected): '
                 + ', '.join(f'{c} ({n})' for c, n in counts.items()) + '</p>')
    return html


def participants_html(wide):
    """Compact per-participant table of the headline metrics + provenance."""
    if wide is None or not len(wide):
        return '<p><i>(no participant)</i></p>'
    cols = ['file_id', 'tib_min', 'sol_min', 'tst_min', 'waso_min', 'se_pct', 'lat_rem_min',
            'n1_pct', 'n2_pct', 'n3_pct', 'rem_pct', 'ahi', 'plm_index', 'arousal_index', 'odi',
            't90_pct_tst', 'lights_source', 'n_mt_custom', 'n_unmapped_events', 'n_warnings']
    cols = [c for c in cols if c in wide.columns]
    df = wide[cols].copy()
    for c in cols:
        if c in METRIC_INFO:
            df[c] = pd.to_numeric(df[c], errors='coerce').round(1)
    return df.to_html(index=False, border=0, na_rep='')


def glossary_html():
    rows = ''
    for key, label, unit, group, definition in METRICS:
        rng = REFERENCE_RANGES.get(key)
        rtxt = f'{rng[0]}-{rng[1]} {unit}' + (f' <small>({rng[2]})</small>' if rng[2] else '') if rng else ''
        rows += (f'<tr><td><b>{label}</b><br><small><code>{key}</code></small></td><td>{unit}</td>'
                 f'<td>{definition}</td><td>{rtxt}</td></tr>')
    return ('<p>Indicative ranges are values for healthy adults recorded in a laboratory; '
            'they depend on age, medication, the first-night effect and the scoring rules. Use them '
            'as a sanity check, never as a diagnostic cut-off. Their sources: see <i>References</i>.</p>'
            '<table style="font-size:85%"><tr><th>Metric</th><th>Unit</th><th>Definition</th>'
            f'<th>Indicative range</th></tr>{rows}</table>')


def references_html():
    """Where each indicative value comes from + the reference list (notebook end + database report)."""
    rows = ''.join(f'<tr><td><b>{metric}</b></td><td>{value}</td><td>{source}</td><td>{remark}</td></tr>'
                   for metric, value, source, remark in REFERENCE_SOURCES)
    refs = ''
    for short, citation, doi in REFERENCES:
        link = f' doi:<a href="https://doi.org/{doi}">{doi}</a>' if doi else ''
        refs += f'<li><b>{short}</b> - {citation}{link}</li>'
    return ('<p><b>Read with care.</b> The Ohayon 2004 values are read on the fitted age curves of its '
            'Figure 1 (approximate). They come from healthy, screened subjects recorded in a laboratory, '
            'first (habituation) nights excluded, scored with the Rechtschaffen &amp; Kales rules '
            '(slow-wave sleep = S3 + S4, close to but not identical with AASM N3). Most metrics change with '
            'age: the green band is the young-to-middle-aged adult range, the remark gives the values '
            'at older ages. <i>Convention</i> = commonly used cut-off without a normative study.</p>'
            '<table style="font-size:85%"><tr><th>Metric</th><th>Indicative value</th><th>Source</th>'
            f'<th>Remark</th></tr>{rows}</table>'
            f'<p><b>References</b></p><ol style="font-size:85%">{refs}</ol>')


def guidance_html(cfg, data_dir, reports_dir, has_subject_info=False, has_comparison=False):
    """'Where to look': which file answers which question."""
    join_note = (' The participant-info columns are joined on the right of every row, so the '
                 'table can be regrouped by any population variable (e.g. <code>groupby("group")'
                 '["se_pct"].describe()</code>).' if has_subject_info else
                 ' Select a participant table in Section 1 to have its columns joined on the right '
                 'of every row (grouping by population variables).')
    cmp_note = ('<li><b><code>comparison_with_reference.tsv</code></b> - one row per (participant, '
                'compared metric): <code>ours</code>, <code>reference</code>, <code>diff</code>, the '
                'tolerance and a <code>status</code> (<code>match</code> / <code>mismatch</code> / '
                '<code>missing_ours</code> / <code>missing_reference</code> / '
                '<code>not_in_reference_table</code>). Filter on <code>mismatch</code> first: a '
                'systematic offset on one metric usually means a definition difference (lights, '
                'sleep-onset rule, MT handling), a single outlier a scoring or file problem.</li>'
                if has_comparison else '')
    return f'''
<p><b>Start with</b> <code>{reports_dir / 'global_sleep_metrics.tsv'}</code> - one row per participant, the
table for statistics.</p>
<ul>
<li><b><code>global_sleep_metrics.tsv</code></b> - <code>file_id</code>, then one column per metric
(<code>{', '.join(METRIC_KEYS[:8])}, ...</code>; units in the glossary: minutes, %, events/h, counts),
then the provenance columns (<code>lights_source</code>, <code>lights_off_clock</code>,
<code>lights_on_clock</code>, <code>sleep_onset_rule</code>, <code>n_mt_custom</code>,
<code>event_source</code>, <code>n_unmapped_events</code>, <code>odi_rule</code>, <code>n_warnings</code>
...).{join_note} Before trusting a value: check <code>lights_source</code> (a
<code>recording_bounds</code> row has TIB = whole recording, so SOL and SE are biased), <code>n_warnings</code>
and <code>n_unmapped_events</code> (events absent from <code>event_remap.json</code> are NOT counted in
the AHI / PLM / arousal indices - map them with tool 4 and re-run with Skip unticked).</li>
<li><b><code>global_sleep_metrics_checks.tsv</code></b> - every consistency check of every participant
(<code>check</code>, <code>status</code> ok / info / warning / fail, <code>detail</code>). A
<code>fail</code> means a metric could not be trusted or computed (no sleep, hypnogram/EDF mismatch);
a <code>warning</code> flags a fallback or an unusual situation (lights fallback, MT epochs, unmapped
events, SpO2 artefact ...).</li>
<li><b><code>sleep_metrics_database.xlsx</code></b> - the same tables as sheets
(<code>sleep_metrics</code>, <code>checks</code>, <code>glossary</code>, <code>parameters</code>,
<code>failed</code>); the TSV files stay the source of truth.</li>
<li><b><code>sleep_metrics_failed.tsv</code></b> - participants that produced no metrics and why
(missing hypnogram, length mismatch, unreadable EDF header).</li>
{cmp_note}
<li><b><code>{data_dir}/&lt;subtree&gt;/{{file_id}}_sleep_metrics.tsv</code></b> - the per-participant
long table (<code>metric</code>, <code>value</code>, <code>unit</code>, <code>definition_short</code>)
the global table is rebuilt from; <code>{{file_id}}_sleep_metrics_checks.tsv</code> next to it. Delete
a participant's <code>_sleep_metrics.tsv</code> to have it reprocessed on the next run.</li>
</ul>
<p><b>Parameters of this run</b>: lights source <code>{cfg.get('lights_source')}</code>, sleep-onset rule
<code>{cfg.get('onset_rule')}</code>, epoch length {cfg.get('epoch_s')} s, MT / custom policy
<code>{OTHER_STAGE_POLICY}</code>, event window <code>{EVENT_COUNT_WINDOW}</code>, ODI
{'on (threshold ' + str(cfg.get('odi_thr')) + ' %)' if cfg.get('odi') else 'off'}, T90
{'on (< ' + str(cfg.get('spo2_thr')) + ' %, floor ' + str(cfg.get('spo2_floor')) + ' %)' if cfg.get('t90') else 'off'}.</p>
'''


def _clock(start_dt, epoch_idx, epoch_s):
    return (start_dt + datetime.timedelta(seconds=int(epoch_idx) * epoch_s)).strftime('%H:%M:%S')


def compute_recording(p, cfg, custom_stages, event_remap, subject, hypno_suffix, evt_txt_suffix,
                      evt_csv_suffix, phase=None, stages=None):
    """The computation of tool 5's process_participant for one recording, moved unchanged: EDF header,
    hypnogram, lights-off / lights-on, macrostructure, scored events (+ ODI), SpO2 (T90 / T88), consistency
    checks and provenance. It writes nothing (tool 5 writes the tables, tool 0 shows them).
    `p`: the scanned recording (file_id, edf, header, hypno, summary_txt, evt_xml, light_ch, spo2_ch);
    `cfg`: tool 5's parameters (get_params); `subject`: its participant-table row or None;
    `phase(name, add)`: progress callback (tool 5's progress bar), None = silent; `stages`: a hypnogram
    to use instead of the file (tool 0: the hypnogram in use, rescorings included). A fatal step raises.
    Returns a dict: metrics, checks, provenance, lights, extras, event_info, cls, stages, header."""
    if phase is None:
        phase = lambda name, add: None
    fid = p['file_id']
    checks = []
    epoch_s = cfg['epoch_s']
    # --- header + hypnogram -------------------------------------------------
    phase('reading the EDF header', 0)
    header = p.get('header') or read_edf_header_info(p['edf'])
    phase('loading the hypnogram', COST_HEADER)
    if stages is None:
        if p.get('hypno') is None:
            raise FileNotFoundError(f'no hypnogram {fid}{hypno_suffix} next to the EDF')
        stages = load_hypnogram(p['hypno'])
    stages, chk = reconcile_hypno_length(stages, header['duration_s'], epoch_s)
    checks.append(chk)
    if chk['status'] == 'fail':
        raise ValueError(chk['detail'])
    custom = custom_stages
    cls = classify_stages(stages, custom)
    n_hyp = len(stages)

    # --- lights ---------------------------------------------------------------
    phase('resolving lights-off / lights-on', COST_HYPNO)
    light_ch = cfg['light_ch'] if cfg['light_ch'] and cfg['light_ch'] in header['ch_names'] else p.get('light_ch')
    lights = resolve_lights(cfg, {**p, 'light_ch': light_ch}, header, n_hyp, cls['sleep_mask'],
                            subject)
    checks += lights['checks']
    off_idx, on_idx = lights['off_idx'], lights['on_idx']

    # --- macrostructure -------------------------------------------------------
    metrics, mchecks, extras = compute_macrostructure(stages, cls, off_idx, on_idx, epoch_s, cfg['onset_rule'])
    checks += mchecks
    tst = metrics.get('tst_min', np.nan)

    # --- events ---------------------------------------------------------------
    phase('counting the scored events', COST_LIGHTS)
    event_info = {'enabled': cfg['events'] or cfg['odi'], 'source': None, 'n_total': 0,
                  'n_outside': 0, 'n_bad_onset': 0, 'n_unmapped': 0}
    odi_rule = ''
    if cfg['events'] or cfg['odi']:
        remap = event_remap
        events_df, source = load_events(p['edf'], evt_txt_suffix, evt_csv_suffix)
        event_info['source'] = source
        if events_df is not None:
            event_info['n_total'] = int(len(events_df))
            win_df, n_out, n_bad = select_events_in_window(events_df, cls, off_idx, on_idx, epoch_s,
                                                           EVENT_COUNT_WINDOW)
            event_info.update({'n_outside': n_out, 'n_bad_onset': n_bad})
            if remap is None:
                checks.append({'check': 'event_remap', 'status': 'warning',
                               'detail': 'no event_remap.json - the scored events could not be '
                                         'classified (run tool 4); event indices are NaN'})
            else:
                if cfg['events']:
                    em, ec, n_unm = compute_event_indices(win_df, remap, tst)
                    metrics.update(em)
                    checks += ec
                    event_info['n_unmapped'] = n_unm
                if cfg['odi']:
                    om, oc, odi_rule = compute_odi(win_df, remap, tst, cfg['odi_thr'], p.get('evt_xml'),
                                                   cls, off_idx, on_idx, epoch_s, EVENT_COUNT_WINDOW)
                    metrics.update(om)
                    checks += oc

    # --- SpO2 -----------------------------------------------------------------
    spo2_ch = ''
    if cfg['t90']:
        phase('reading the SpO2 signal', COST_EVENTS if (cfg['events'] or cfg['odi']) else 0)
        spo2_ch = cfg['spo2_ch'] if cfg['spo2_ch'] and cfg['spo2_ch'] in header['ch_names'] else p.get('spo2_ch')
        if not spo2_ch:
            checks.append({'check': 'spo2_channel', 'status': 'warning',
                           'detail': 'no SpO2 channel found - T90 not computed'})
            spo2_ch = ''
        else:
            try:
                x, sf, rescaled = load_spo2_signal(p['edf'], spo2_ch)
                if rescaled:
                    checks.append({'check': 'spo2_units', 'status': 'warning',
                                   'detail': 'SpO2 stored as a fraction (max <= 1) - rescaled to %'})
                tm, tc = compute_t90(x, sf, off_idx, on_idx, epoch_s, tst, cfg['spo2_thr'], cfg['spo2_floor'],
                                     cls['sleep_mask'])
                metrics.update(tm)
                checks += tc
                del x
            except Exception as exc:
                checks.append({'check': 'spo2_channel', 'status': 'warning',
                               'detail': f'SpO2 signal unreadable ({exc}) - T90 not computed'})

    # --- checks, provenance, write --------------------------------------------
    phase('checks and writing', (COST_SPO2 if cfg['t90'] else 0)
          + (0 if cfg['t90'] else (COST_EVENTS if (cfg['events'] or cfg['odi']) else 0)))
    checks += run_checks(metrics, extras, lights, n_hyp, cls, event_info)
    n_warn = sum(1 for c in checks if c['status'] in ('warning', 'fail'))
    provenance = {
        'lights_source': lights['source'],
        'lights_off_clock': _clock(header['start_dt'], off_idx, epoch_s),
        'lights_on_clock': _clock(header['start_dt'], on_idx, epoch_s),
        'lights_off_epoch': off_idx, 'lights_on_epoch': on_idx,
        'sleep_onset_rule': cfg['onset_rule'],
        'sleep_onset_epoch': extras['onset_idx_abs'] if extras['onset_idx_abs'] is not None else '',
        'other_stage_policy': OTHER_STAGE_POLICY, 'event_count_window': EVENT_COUNT_WINDOW,
        'epoch_length_s': epoch_s, 'n_epochs_hypno': n_hyp, 'n_epochs_in_bed': extras['n_epochs_in_bed'],
        'n_mt_custom': int(cls['other_mask'].sum()),
        'event_source': event_info['source'] or '', 'n_events_total': event_info['n_total'],
        'n_unmapped_events': event_info['n_unmapped'], 'n_events_outside_window': event_info['n_outside'],
        'has_spo2': bool(spo2_ch), 'odi_threshold_pct': cfg['odi_thr'] if cfg['odi'] else '',
        'odi_rule': odi_rule, 'spo2_channel': spo2_ch, 'light_channel': lights.get('light_channel', ''),
        'recording_start': header['start_dt'].strftime('%Y-%m-%d %H:%M:%S'),
        'edf_duration_min': round(header['duration_s'] / 60, 2), 'n_warnings': n_warn}
    return {'metrics': metrics, 'checks': checks, 'provenance': provenance, 'lights': lights,
            'extras': extras, 'event_info': event_info, 'cls': cls, 'stages': stages, 'header': header}
