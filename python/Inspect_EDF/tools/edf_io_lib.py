"""
Shared EDF reading helpers: header, channel types, configuration files, hypnogram, scored events,
context channels. EDF-specific (the Curry twins have tools_curry/curry_header.py + curry_io.py).

Part of the tool-0 redesign (tools/plan_tool0_redesign.md): tool 0 walks ONE participant through the whole
pipeline and must read a recording exactly like the batch tools do, so the reading code lives here once.
Every function below is MOVED VERBATIM from the tool named in its section header; when the code was
written inline in a notebook loop, it is lifted unchanged into a small function.

Consumed by:
  - 0_live_explore_1file_voila.ipynb
  - (phase 8 of the plan) the batch tools named below switch to these functions; until then their own copy
    is the reference, and tests/test_tool0.py compares tool 0's results with theirs.
"""
import datetime
import json
import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import chardet
import numpy as np
import pandas as pd
import mne


# ---------------------------------------------------------------------------------------------------------
# EDF header (from 1_inspect_edf_voila; same parser in 1bis / 2 / 0)
# ---------------------------------------------------------------------------------------------------------
def detect_encoding(byte_string, min_confidence=0.6, fallback=None):
    """Encoding of the 256-byte EDF header. Tool 1 raises when chardet is not confident; tool 0 passes
    fallback='latin-1' instead (latin-1 maps every byte, so an odd header can still be displayed)."""
    result = chardet.detect(byte_string)
    encoding = result['encoding']
    confidence = result['confidence']
    if encoding is None or confidence < min_confidence:
        if fallback is not None:
            return fallback
        raise UnicodeDecodeError("chardet", byte_string, 0, len(byte_string),
                                 f"\tUnable to reliably detect encoding. Detected: {encoding} with confidence {confidence}")
    return encoding


def read_edf_header_custom(file_path, fallback_encoding=None):
    """Read the EDF header octet by octet (never the signal). Returns a dict: global fields + one list per
    channel field. `sampling_frequency` = samples_per_record / duration_data_record (the 8-byte field is
    the number of samples PER DATA RECORD, not the rate; the two coincide only for 1 s records)."""
    with open(file_path, 'rb') as f:
        header = {}
        raw_header = f.read(256)
        encoding = detect_encoding(raw_header, fallback=fallback_encoding)
        f.seek(0)
        header['version'] = f.read(8).decode(encoding).strip()
        header['patient_id'] = f.read(80).decode(encoding).strip()
        header['recording_id'] = f.read(80).decode(encoding).strip()
        header['start_date'] = f.read(8).decode(encoding).strip()
        header['start_time'] = f.read(8).decode(encoding).strip()
        header['header_bytes'] = int(f.read(8).decode(encoding).strip())
        header['reserved'] = f.read(44).decode(encoding).strip()
        header['n_data_records'] = int(f.read(8).decode(encoding).strip())
        header['duration_data_record'] = float(f.read(8).decode(encoding).strip())
        header['n_channels'] = int(f.read(4).decode(encoding).strip())

        n = header['n_channels']
        channel_fields = {
            'channel': [],
            'transducer_type': [],
            'dimension': [],
            'physical_min': [],
            'physical_max': [],
            'digital_min': [],
            'digital_max': [],
            'prefiltering': [],
            'samples_per_record': [],
            'reserved': [],
        }
        for key in channel_fields:
            length = {
                'channel': 16,
                'transducer_type': 80,
                'dimension': 8,
                'physical_min': 8,
                'physical_max': 8,
                'digital_min': 8,
                'digital_max': 8,
                'prefiltering': 80,
                'samples_per_record': 8,
                'reserved': 32,
            }[key]
            channel_fields[key] = [f.read(length).decode(encoding).strip() for _ in range(n)]
        header.update(channel_fields)

        dur = header['duration_data_record']
        sf_list = []
        for s in header['samples_per_record']:
            try:
                sf = float(s) / dur
                sf_list.append(str(int(sf)) if sf == int(sf) else str(sf))
            except (ValueError, ZeroDivisionError):
                sf_list.append(s)  # fall back to the raw field if non-numeric or zero duration
        header['sampling_frequency'] = sf_list
    return header


def extract_filter_value(s, tag):
    """Value of a hardware filter written in the EDF 'prefiltering' field, e.g. 'HP:0.3Hz LP:35Hz'."""
    if pd.isna(s):
        return None
    match = re.search(rf'{tag}[:\s]*([\d\.]+)\s*', s, re.IGNORECASE)
    return float(match.group(1)) if match else None


def check_anonymization(patient_id, file_stem):
    """
    Inspect the EDF+ 'Local Patient ID' field for a real patient name.
    EDF+ format: 'code sex birthdate name' (space-separated).
    Compumedics encodes the name as 'LASTNAME_FIRSTNAME'; anonymized = 'X_X'.
    The name subfield (4th token onward) is isolated, placeholder chars
    (X/x/_ and whitespace) are stripped; anything remaining is a real name.
    Filename cross-check: each name token with >=3 chars is searched in the
    file stem (case-insensitive) to detect PII leaking into the file name.
    """
    parts = str(patient_id).split()
    if len(parts) >= 4:
        name_subfield = ' '.join(parts[3:])
        id_format = 'edf+'
    else:
        name_subfield = str(patient_id).strip()
        id_format = 'non-standard'
    cleaned = re.sub(r'[Xx_,;\s]', '', name_subfield)
    header_anonymized = (cleaned == '')
    name_in_filename = False
    if not header_anonymized:
        stem_low = file_stem.lower()
        for token in re.split(r'[_\s]+', name_subfield):
            if len(token) >= 3 and token.lower() in stem_low:
                name_in_filename = True
                break
    warning = ''
    if not header_anonymized:
        if name_in_filename:
            warning = 'PII in header AND file name'
        else:
            warning = 'header NOT anonymized (file name looks clean)'
    return {
        'patient_id'        : str(patient_id).strip(),
        'name_subfield'     : name_subfield,
        'patient_id_format' : id_format,
        'header_anonymized' : header_anonymized,
        'name_in_filename'  : name_in_filename,
        'anon_warning'      : warning,
    }


# Tool 1 flags an EEG channel whose exported dynamic range (|physical min| + |physical max|, µV) is at or
# below this value: a ±250 µV range clips large slow waves and movement artefacts.
DYNAMIC_RANGE_THRESHOLD_UV = 500.0


def header_channel_table(header):
    """One row per channel of the header, with tool 1's derived columns (lifted from its inspection loop):
    res_theoretical = (|physical_min| + |physical_max|) / 2^16 (EDF stores 16-bit samples), in µV when the
    dimension is mV; highpass / lowpass / notch read from the prefiltering field."""
    df = pd.DataFrame({k: header[k] for k in ('channel', 'transducer_type', 'dimension', 'physical_min',
                                              'physical_max', 'digital_min', 'digital_max', 'prefiltering',
                                              'samples_per_record', 'sampling_frequency')})
    df['res_theoretical'] = (abs(pd.to_numeric(df['physical_min'])) + abs(pd.to_numeric(df['physical_max']))) / pow(2, 16)
    df.loc[df['dimension'].str.contains('mv', case=False, na=False), 'res_theoretical'] *= 1000
    df['lowpass'] = df['prefiltering'].apply(lambda x: extract_filter_value(x, 'LP'))
    df['highpass'] = df['prefiltering'].apply(lambda x: extract_filter_value(x, 'HP'))
    df['notch'] = df['prefiltering'].apply(lambda x: extract_filter_value(x, 'NOTCH'))
    return df


# ---------------------------------------------------------------------------------------------------------
# Channel types (detect_channel_types from 0_live_explore_1file; detect_context_channels from
# 2_select&remap_channels_edf, Section 2bis)
# ---------------------------------------------------------------------------------------------------------
# Anchored full 10-10 system + mastoids + older 10-20 temporal labels (T3/T4/T5/T6) + literal 'EEG'.
KNOWN_EEG_CHANNEL_RE = (
    r'\bFp1\b|\bFpz\b|\bFp2\b|\bAF7\b|\bAF3\b|\bAFz\b|\bAF4\b|\bAF8\b|\bF7\b|\bF5\b|\bF3\b|\bF1\b'
    r'|\bFz\b|\bF2\b|\bF4\b|\bF6\b|\bF8\b|\bFT7\b|\bFC5\b|\bFC3\b|\bFC1\b|\bFCz\b|\bFC2\b|\bFC4\b'
    r'|\bFC6\b|\bFT8\b|\bT7\b|\bC5\b|\bC3\b|\bC1\b|\bCz\b|\bC2\b|\bC4\b|\bC6\b|\bT8\b|\bTP7\b|\bCP5\b'
    r'|\bCP3\b|\bCP1\b|\bCPz\b|\bCP2\b|\bCP4\b|\bCP6\b|\bTP8\b|\bP7\b|\bP5\b|\bP3\b|\bP1\b|\bPz\b'
    r'|\bP2\b|\bP4\b|\bP6\b|\bP8\b|\bPO7\b|\bPO5\b|\bPO3\b|\bPOz\b|\bPO4\b|\bPO6\b|\bPO8\b|\bO1\b'
    r'|\bOz\b|\bO2\b|\bM1\b|\bM2\b|\bT3\b|\bT4\b|\bT5\b|\bT6\b|EEG'
)


def detect_channel_types(channels, transducers):
    """Classify each channel as EEG / EOG / EMG / ECG / other from name + transducer type.
    EOG: transducer~EOG or name~eog. EMG: transducer~EMG or name~emg|chin|menton.
    ECG: transducer~ECG or name~ecg. EEG: known 10-10 name or transducer EEG/AGAGCL,
    excluding anything already matched as eog/ecg/emg. Returns dict name -> type."""
    types = {}
    for ch, td in zip(channels, transducers):
        td = td or ''
        is_eog = bool(re.search(r'EOG', td, re.I)) or bool(re.search(r'eog', ch, re.I))
        is_emg = bool(re.search(r'EMG', td, re.I)) or bool(re.search(r'emg|chin|menton', ch, re.I))
        is_ecg = bool(re.search(r'ECG|EKG', td, re.I)) or bool(re.search(r'ecg|ekg', ch, re.I))
        is_eeg = (bool(re.search(r'\bEEG\b|\bAGAGCL ELECTRODE\b', td, re.I))
                  or bool(re.search(KNOWN_EEG_CHANNEL_RE, ch, re.I)))
        if is_eog:
            types[ch] = 'EOG'
        elif is_emg:
            types[ch] = 'EMG'
        elif is_ecg:
            types[ch] = 'ECG'
        elif is_eeg:
            types[ch] = 'EEG'
        else:
            types[ch] = 'other'
    return types


def detect_context_channels(df_sub):
    """Guess EOG-left/right, EMG and ECG channels from the scanned header dataframe of one
    configuration's participants. Follows the toolkit convention 'transducer type OR channel
    name' but works on the header scan (tool 2 never loads signal). Tolerant of a missing
    'transducer_type' column so the generated Curry twin keeps working (name-only there)."""
    chans = list(dict.fromkeys(df_sub['channel'].tolist()))
    if 'transducer_type' in df_sub.columns:
        trans = dict(zip(df_sub['channel'], df_sub['transducer_type'].astype(str)))
    else:
        trans = {}
    def td(ch):
        return trans.get(ch, '') or ''
    eog = [c for c in chans
           if re.search(r'\bEOG\b', td(c), re.I)
           or re.search(r'LOC|ROC|\bE1\b|\bE2\b|eog', c, re.I)]
    eog_left  = next((c for c in eog if re.search(r'E1|LOC|left|gauche', c, re.I)), None)
    eog_right = next((c for c in eog if re.search(r'E2|ROC|right|droit', c, re.I)), None)
    for c in eog:
        if c in (eog_left, eog_right):
            continue
        if eog_left is None:
            eog_left = c
        elif eog_right is None:
            eog_right = c
    emg = next((c for c in chans
                if re.search(r'\bEMG\b', td(c), re.I)
                or re.search(r'emg|chin|menton', c, re.I)), None)
    ecg = next((c for c in chans
                if re.search(r'\bE[CK]G\b', td(c), re.I)
                or re.search(r'ecg|ekg', c, re.I)), None)
    return {'eog_left': eog_left, 'eog_right': eog_right, 'emg': emg, 'ecg': ecg}


# ---------------------------------------------------------------------------------------------------------
# Configuration files (lenient JSON loader shared by tools 0, 2, 4)
# ---------------------------------------------------------------------------------------------------------
def load_json_lenient(path):
    """Parse JSON strictly; on failure, repair a single trailing comma before }/] and retry."""
    text = Path(path).read_text(encoding='utf-8')
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        repaired = re.sub(r',(\s*[}\]])', r'\1', text)
        return json.loads(repaired)


def participant_config(config, file_id):
    """The entry of remap_reref_persubject.json (tool 2) for one recording, matched on the file name stem
    with os.path.normcase on both sides (tool 0 rule). None when the participant is not in the file."""
    key = next((k for k in config if os.path.normcase(k) == os.path.normcase(file_id)), None)
    return None if key is None else config[key]


# ---------------------------------------------------------------------------------------------------------
# MNE reading helpers (from 7_preprocessing_voila / 6_quality_overview_voila)
# ---------------------------------------------------------------------------------------------------------
def drop_suffix_duplicates(raw):
    """Keep only the -0 variant when MNE creates -0/-1 duplicates for repeated channel names."""
    groups = {}
    for ch in raw.ch_names:
        if ch.endswith('-0') or ch.endswith('-1'):
            base = ch.rsplit('-', 1)[0]
            groups.setdefault(base, []).append(ch)
    to_drop = []
    for base, ch_list in groups.items():
        if any(c.endswith('-0') for c in ch_list):
            to_drop.extend(c for c in ch_list if not c.endswith('-0'))
    if to_drop:
        raw.drop_channels(to_drop)
    return raw, to_drop


def adapt_remap_dict_to_suffixes(raw, remap_dict):
    """Handle MNE's -0 suffix when the remap config uses the base channel name."""
    ch_set = set(raw.ch_names)
    new_remap = {}
    for base_label, target in remap_dict.items():
        if base_label in ch_set:
            new_remap[base_label] = target
        elif f'{base_label}-0' in ch_set:
            new_remap[f'{base_label}-0'] = target
    return new_remap


def get_phys_bounds_uV(extras, ch_idx):
    """Reconstruct physical_min/max (µV) from MNE _raw_extras (MNE 1.9 stores physical_max + offset).
    physical_max/offset are in the channel's native EDF unit (uV, mV, V...); extras['units'] is the
    multiplier to volts (1e-6 uV, 1e-3 mV, 1.0 V), so * units * 1e6 converts to uV to match get_data()*1e6.
    Without it an mV channel (EOG/EMG/ECG, physical_max=1.0 mV) gets a 1 uV bound and flags ~100% of samples."""
    to_uV = float(extras['units'][ch_idx]) * 1e6
    phys_max = float(extras['physical_max'][ch_idx]) * to_uV
    offset = float(extras['offsets'][ch_idx]) * to_uV
    return -phys_max + 2 * offset, phys_max


def read_eeg_raw(edf_path, remap):
    """Tool 7's EEG read ([load] + [B]): `include=` evaluated AT READ TIME with the original names of the
    tool-2 remap (excludes faster off-montage channels, which avoids MNE's partial-read AssertionError and
    keeps the native EEG rate), suffix duplicates dropped, channels renamed to their harmonised names, then
    the data loaded. Returns the preloaded Raw (no re-reference, no filter)."""
    raw = mne.io.read_raw_edf(str(edf_path), preload=False, encoding='latin-1',
                              include=list(remap.keys()), verbose=False)
    raw, _ = drop_suffix_duplicates(raw)
    raw.rename_channels(adapt_remap_dict_to_suffixes(raw, remap))
    raw.load_data()
    return raw


# ---------------------------------------------------------------------------------------------------------
# Recordings and companions (detection rules of 7_preprocessing_voila, lifted from _detect_hypno_suffixes)
# ---------------------------------------------------------------------------------------------------------
def list_edf_files(folder):
    """Every .edf under `folder` (recursive), macOS '._' shadow files excluded, sorted."""
    return [f for f in sorted(Path(folder).rglob('*')) if f.suffix.lower() == '.edf' and not f.name.startswith('._')]


def companion_suffix_counts(edf_files, files):
    """{suffix: number of recordings having a file named <recording stem><suffix>} (normcase prefix)."""
    counts = {}
    for edf in edf_files:
        for f in files:
            if os.path.normcase(f.name).startswith(os.path.normcase(edf.stem)):
                suf = f.name[len(edf.stem):]
                counts[suf] = counts.get(suf, 0) + 1
    return counts


def select_hypno_suffix(suffix_counts):
    """Tool 7's hypnogram rule: among the .txt suffixes found for >= 50 % of the maximum count, prefer the
    LONGEST (the more specific remapped/processed version), ties broken by count. Suffixes containing
    'event' are skipped (Curry/Profusion write *_ScoredEvents_Export.txt next to the hypnograms; its
    longer suffix would otherwise win) unless nothing else exists. Returns (suffix, count) or (None, 0)."""
    if not suffix_counts:
        return None, 0
    sel_counts = {s: c for s, c in suffix_counts.items() if 'event' not in s.lower()}
    if not sel_counts:
        sel_counts = suffix_counts  # only event exports found: fall back rather than fail
    max_count = max(sel_counts.values())
    candidates = {s: c for s, c in sel_counts.items() if c >= max_count * 0.5}
    return max(candidates.items(), key=lambda x: (len(x[0]), x[1]))


def select_event_suffix(suffix_counts):
    """Tool 7's event TXT / CSV rule: the most frequent suffix, the shortest on ties. (None, 0) if none."""
    if not suffix_counts:
        return None, 0
    return min(suffix_counts.items(), key=lambda x: (-x[1], len(x[0])))


# ---------------------------------------------------------------------------------------------------------
# Hypnogram (7_preprocessing_voila)
# ---------------------------------------------------------------------------------------------------------
def load_hypno_labels(path):
    """One stage label per 30 s epoch, as tool 7 reads it."""
    return np.loadtxt(str(path), dtype=str).astype('<U10')


def hypno_length_status(n_hyp, n_epochs):
    """Tool 7's length rule between a hypnogram and the complete 30 s epochs of the signal:
    'ok' (equal), 'trim' (one extra score: the recording ends mid-epoch, the last score is dropped),
    'mismatch' (anything else: tool 7 skips the participant)."""
    n_over = n_hyp - n_epochs
    if n_over == 0:
        return 'ok'
    if 0 < n_over <= 1:
        return 'trim'
    return 'mismatch'


def mid_recording_unscored(labels):
    """Tool 3's Scan rule (3_remap_hypno_voila): the unscored '?' epochs more than 10 epochs away from both
    ends of the night. A '?' at the very start or end is lights-off / lights-on; one in the middle of the
    night is an epoch the scorer skipped, to look at. Returns their 0-based indices."""
    hypno = np.asarray(labels)
    q_idx = np.where(hypno == '?')[0]
    return [int(i) for i in q_idx if i > 10 and i < len(hypno) - 10]


# ---------------------------------------------------------------------------------------------------------
# Scored events (from 5_compute_sleep_macrostructure_voila: tool 7's TXT -> CSV -> XML chain, plus the
# additive Desaturation column of the XML parser)
# ---------------------------------------------------------------------------------------------------------
_TXT_DUR_RE = re.compile(r'^(\d+):(\d+(?:\.\d+)?)$')   # duration column "M:SS" or "M:SS.s"


def read_edf_start_datetime(edf_path):
    """Read the EDF recording-start datetime from the fixed header (offset 168 = date
    'dd.mm.yy', 176 = time 'hh.mm.ss'), applying the EDF 2-digit-year clipping
    (00-84 -> 20xx, 85-99 -> 19xx). Header-only read; returns a datetime or None on failure."""
    try:
        with open(edf_path, 'rb') as f:
            f.seek(168)
            date_str = f.read(8).decode('ascii', 'replace').strip()   # dd.mm.yy
            time_str = f.read(8).decode('ascii', 'replace').strip()   # hh.mm.ss
        dd, mm, yy = (int(x) for x in date_str.split('.'))
        hh, mi, ss = (int(x) for x in time_str.split('.'))
        year = 2000 + yy if yy <= 84 else 1900 + yy
        return datetime.datetime(year, mm, dd, hh, mi, ss)
    except Exception:
        return None


def event_companion_paths(edf_path, txt_suffix='_ScoredEvents_Export.txt', csv_suffix='_event_xml.csv'):
    """Return (txt_path_or_None, csv_path_or_None, xml_path_or_None) for an EDF stem."""
    edf_path = Path(edf_path)
    txt = edf_path.with_name(f'{edf_path.stem}{txt_suffix}')
    txt = txt if txt.exists() else None
    csv = edf_path.with_name(f'{edf_path.stem}{csv_suffix}')
    csv = csv if csv.exists() else None
    xml = None
    for cand in (f'{edf_path.name}.XML', f'{edf_path.name}.xml'):
        p = edf_path.with_name(cand)
        if p.exists():
            xml = p
            break
    return txt, csv, xml


def _events_df_from_txt(txt_path, rec_start):
    """Parse a Compumedics/Curry French text export (*_ScoredEvents_Export.txt)
    -> DataFrame Name/Start/Duration (seconds), or None. Comma-separated, no header, columns:
        HH:MM:SS , epoch# , stage_FR , event_label_FR , M:SS[.s] , - , - , position
    Encoding varies (UTF-16 with BOM on some exports, UTF-8/ANSI on others) -> the BOM is sniffed.
    Clock times are converted to seconds-from-recording-start (with midnight rollover) using rec_start."""
    raw = open(txt_path, 'rb').read()
    if raw[:2] in (b'\xff\xfe', b'\xfe\xff'):
        text = raw.decode('utf-16')
    else:
        text = raw.decode('utf-8', errors='replace')
    rec_date = rec_start.date()
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = [p.strip() for p in line.split(',')]
        if len(parts) < 5:
            continue
        name = parts[3]
        try:
            hh, mm, ss = parts[0].split(':')   # hours may be single-digit
            clock_t = datetime.time(int(hh), int(mm), int(ss))
        except (ValueError, TypeError):
            continue
        event_dt = datetime.datetime.combine(rec_date, clock_t)
        # events recorded after midnight fall on the next calendar day
        if (rec_start - event_dt).total_seconds() > 3600:
            event_dt += datetime.timedelta(days=1)
        start_sec = (event_dt - rec_start).total_seconds()
        m = _TXT_DUR_RE.match(parts[4])
        dur_sec = int(m.group(1)) * 60 + float(m.group(2)) if m else 0.0
        rows.append((name, start_sec, dur_sec))
    if not rows:
        return None
    return pd.DataFrame(rows, columns=['Name', 'Start', 'Duration'])


def _events_df_from_csv(csv_path):
    """Parse a Compumedics event CSV -> DataFrame Name/Start/Duration, or None."""
    df = pd.read_csv(str(csv_path))
    cols = {c.lower(): c for c in df.columns}
    if 'name' in cols and 'start' in cols and 'duration' in cols:
        return df.rename(columns={cols['name']: 'Name', cols['start']: 'Start',
                                  cols['duration']: 'Duration'})[['Name', 'Start', 'Duration']]
    return None


def _events_df_from_xml(xml_path):
    """Parse the <ScoredEvents> of a Profusion CMPStudyConfig .edf.XML
    -> DataFrame Name/Start/Duration (seconds) + Desaturation (depth in %, NaN when the
    event has no <Desaturation> child - only SpO2 desaturations carry one), or None.
    Tool 7's parser plus the additive Desaturation column (needed by the ODI threshold)."""
    root = ET.parse(str(xml_path)).getroot()
    rows = []
    for se in root.iter('ScoredEvent'):
        name_el = se.find('Name')
        if name_el is None or name_el.text is None:
            continue
        start_el = se.find('Start')
        dur_el = se.find('Duration')
        desat_el = se.find('Desaturation')
        start = float(start_el.text) if (start_el is not None and start_el.text) else np.nan
        dur = float(dur_el.text) if (dur_el is not None and dur_el.text) else np.nan
        try:
            desat = float(desat_el.text) if (desat_el is not None and desat_el.text) else np.nan
        except ValueError:
            desat = np.nan
        rows.append((name_el.text.strip(), start, dur, desat))
    if not rows:
        return None
    return pd.DataFrame(rows, columns=['Name', 'Start', 'Duration', 'Desaturation'])


def _events_df_from_manual(manual_path):
    """Read a tool-8 manual-annotation table ({stem}_manual_events.tsv: type / onset_s / duration_s,
    written by the reviewer in tool 8's navigator) -> the same Name/Start/Duration DataFrame as the
    scored companions, or None. These labels are ALREADY canonical (see augment_remap_for_manual)."""
    df = pd.read_csv(manual_path, sep='\t')
    if not {'type', 'onset_s'}.issubset(df.columns):
        return None
    out = pd.DataFrame({'Name': df['type'].astype(str),
                        'Start': pd.to_numeric(df['onset_s'], errors='coerce'),
                        'Duration': (pd.to_numeric(df['duration_s'], errors='coerce')
                                     if 'duration_s' in df.columns else np.nan),
                        'Source': 'manual'})   # tags the rows augment_remap_for_manual may map by identity
    out = out[out['Start'].notna()]
    return out if len(out) else None


def augment_remap_for_manual(event_remap, events_df):
    """Manual annotations already carry a CANONICAL label (tool 8 offers the canonical vocabulary), but
    compute_event_epoch_mask maps every name through event_remap.get(name) and would drop them as
    unknown. Return a copy of the remap with an identity entry per manual label. ONLY the rows tagged
    Source == 'manual' are treated this way: a SCORED label missing from the remap stays unmapped on
    purpose (the user chose not to map it in tool 4). The event_remap.json file itself is never
    rewritten — this lives in memory, for this run only."""
    if events_df is None or not len(events_df) or 'Source' not in events_df.columns:
        return event_remap
    out = dict(event_remap)
    manual = events_df[events_df['Source'].astype(str) == 'manual']
    for name in manual['Name'].astype(str).unique():
        out.setdefault(name, name)
    return out


def load_events(edf_path, txt_suffix='_ScoredEvents_Export.txt', csv_suffix='_event_xml.csv',
                manual_suffix='_manual_events.tsv', include_manual=False):
    """Load scored events next to the EDF, TXT-first then CSV then *.edf.XML fallback.
    Returns (DataFrame Name/Start/Duration in seconds, source) with source in
    {'txt', 'csv', 'xml'}, or (None, None) when no usable event companion is found.
    The .txt onsets need the recording-start datetime (header offsets 168/176), so a .txt present
    without a readable start datetime is skipped in favour of the CSV/XML companions.

    With include_manual=True the tool-8 manual annotations ({stem}_manual_events.tsv) are APPENDED to
    whichever scored source won the priority chain (source tag e.g. 'csv+manual'), or used alone when
    there is no scored companion at all ('manual'). They are an ADDITION to the scored corpus, never an
    alternative source, so the priority chain above is untouched. Default False -> byte-identical to the
    pre-annotation tool."""
    txt, csv, xml = event_companion_paths(edf_path, txt_suffix, csv_suffix)
    df, source = None, None
    if txt is not None:
        rec_start = read_edf_start_datetime(edf_path)
        if rec_start is not None:
            try:
                d = _events_df_from_txt(txt, rec_start)
                if d is not None:
                    df, source = d, 'txt'
            except Exception:
                pass  # fall through to the CSV
    if df is None and csv is not None:
        try:
            d = _events_df_from_csv(csv)
            if d is not None:
                df, source = d, 'csv'
        except Exception:
            pass  # fall through to the XML
    if df is None and xml is not None:
        try:
            d = _events_df_from_xml(xml)
            if d is not None:
                df, source = d, 'xml'
        except Exception:
            pass
    if include_manual:
        # Tool-8 annotations are ADDED to whatever the chain above found (never a replacement).
        manual_path = Path(edf_path).with_name(f'{Path(edf_path).stem}{manual_suffix}')
        if manual_path.exists():
            try:
                m = _events_df_from_manual(manual_path)
                if m is not None:
                    df = m if df is None else pd.concat([df, m], ignore_index=True)
                    source = 'manual' if source is None else f'{source}+manual'
            except Exception:
                pass  # a broken annotation file must never cost us the scored events
    if df is None:
        return None, None
    return df, source


def load_event_remap(path):
    """Lenient loader for event_remap.json (strict parse, then repair one trailing comma
    before a closing } or ]). Returns {raw_label: canonical_label_or_None}."""
    with open(path, encoding='utf-8') as f:
        text = f.read()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return json.loads(re.sub(r',(\s*[}\]])', r'\1', text))


# ---------------------------------------------------------------------------------------------------------
# Context channels (7_preprocessing_voila, block [H]: the EOG / EMG / ECG companion of the EEG epochs)
# ---------------------------------------------------------------------------------------------------------
CONTEXT_ROLES = ('eog_left', 'eog_right', 'emg', 'ecg')
CONTEXT_ROLE_TYPE = {'eog_left': 'eog', 'eog_right': 'eog', 'emg': 'emg', 'ecg': 'ecg'}
CONTEXT_ROLE_LABEL = {'eog_left': 'EOG-L', 'eog_right': 'EOG-R', 'emg': 'EMG', 'ecg': 'ECG'}


def read_context_raw(edf_path, ctx_map, eeg_sfreq):
    """The declared context channels ({role: original EDF name}, tool 2's `context_channels` block), read
    like tool 7 does: include= at read time, renamed to their role labels (EOG-L / EOG-R / EMG / ECG),
    resampled to the EEG rate, then filtered for DISPLAY per role (AASM-like: EOG band-pass 0.3-35 Hz,
    EMG high-pass 10 Hz, ECG band-pass 0.5-40 Hz). Returns the Raw, or None when no declared channel is
    in the file."""
    orig_to_role = {}
    for _role in CONTEXT_ROLES:
        _orig = (ctx_map or {}).get(_role)
        if _orig:
            orig_to_role[_orig] = _role
    if not orig_to_role:
        return None
    ctx_probe = mne.io.read_raw_edf(str(edf_path), preload=False, encoding='latin-1',
                                    include=list(orig_to_role.keys()), verbose=False)
    ctx_probe, _ = drop_suffix_duplicates(ctx_probe)
    present = [o for o in orig_to_role if o in ctx_probe.ch_names]
    if not present:
        return None
    ctx_raw = ctx_probe.pick(present)
    ctx_raw.load_data()
    ctx_raw.rename_channels({o: CONTEXT_ROLE_LABEL[orig_to_role[o]] for o in present})
    ctx_raw.set_channel_types(
        {CONTEXT_ROLE_LABEL[orig_to_role[o]]: CONTEXT_ROLE_TYPE[orig_to_role[o]] for o in present},
        verbose=False)
    if not np.isclose(ctx_raw.info['sfreq'], eeg_sfreq):
        ctx_raw.resample(eeg_sfreq, npad='auto', verbose=False)
    _nyq = ctx_raw.info['sfreq'] / 2.0
    _eog = [CONTEXT_ROLE_LABEL[orig_to_role[o]] for o in present if orig_to_role[o] in ('eog_left', 'eog_right')]
    _emg = [CONTEXT_ROLE_LABEL[orig_to_role[o]] for o in present if orig_to_role[o] == 'emg']
    _ecg = [CONTEXT_ROLE_LABEL[orig_to_role[o]] for o in present if orig_to_role[o] == 'ecg']
    if _eog:
        ctx_raw.filter(0.3, min(35.0, _nyq - 0.5), picks=_eog, verbose=False)
    if _emg:
        ctx_raw.filter(10.0, None, picks=_emg, verbose=False)
    if _ecg:
        ctx_raw.filter(0.5, min(40.0, _nyq - 0.5), picks=_ecg, verbose=False)
    return ctx_raw
