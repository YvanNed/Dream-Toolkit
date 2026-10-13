"""
Automatic sleep scoring with GSSC: the computation and the figures of tool 3bis
(3bis_autoscoring_gssc_voila.ipynb), shared with tool 0 (tool-0 redesign, tools/plan_tool0_redesign.md).

The functions below are MOVED VERBATIM from tool 3bis (its cell 2), plus its guarded GSSC import (cell 1) and
`score_night`, the scoring part of its `score_one` lifted unchanged so both tools chain the same steps.
Tool 3bis keeps its user interface and everything it writes; tool 0 calls the same functions on one night.

GSSC is pinned to 0.0.9: `gssc_combination_logits` reproduces gssc's own scoring loop (see SPEC §3bis). Without
torch / gssc the module still imports (GSSC_OK = False, GSSC_MSG says how to install them), so the plots and
the comparison remain usable on a scoring read from disk.
"""
import os
import re
import json
import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mne
from sklearn.metrics import accuracy_score, cohen_kappa_score, f1_score, confusion_matrix

# GSSC and torch are imported apart: without them nothing can be scored, but the rest still works.
# torch >= 2.6 loads model files with weights_only=True by default, which refuses the weights shipped
# inside the gssc package. They come with the installed package (not from an untrusted download), so
# the former behaviour is restored for this process before GSSC loads them.
os.environ['TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD'] = '1'
GSSC_OK = False
GSSC_MSG = ''
try:
    import torch
    from gssc.infer import EEGInfer
    from gssc.utils import (prepare_inst as gssc_prepare_inst, permute_sigs as gssc_permute_sigs,
                            epo_arr_zscore as gssc_zscore, check_flip_chan as gssc_check_flip_chan,
                            loudest_vote as gssc_loudest_vote)
    from importlib.metadata import version as _pkg_version
    GSSC_VERSION = _pkg_version('gssc')
    TORCH_VERSION = torch.__version__
    GSSC_OK = True
except Exception as e:
    GSSC_VERSION, TORCH_VERSION = None, None
    GSSC_MSG = (
        '<div style="border:2px solid #c0392b;background:#fdecea;padding:8px 12px;border-radius:4px;">'
        f'<b>GSSC is not available</b> ({e}). Install it in the <code>inspect_edf</code> environment, '
        'then restart this page:<br><code>conda install -n inspect_edf -c conda-forge pytorch-cpu</code><br>'
        '<code>python -m pip install gssc==0.0.9</code></div>')


EPOCH_SEC = 30
STAGES = ['W', 'N1', 'N2', 'N3', 'R']
GSSC_STAGE = {0: 'W', 1: 'N1', 2: 'N2', 3: 'N3', 4: 'R'}       # GSSC's integer coding
# Clinical hypnogram layout: W on top, then R, then NREM getting deeper downwards.
HYPNO_Y = {'W': 0, 'R': 1, 'N1': 2, 'N2': 3, 'N3': 4}
# Stage colours of the rest of the toolkit (BASE_STAGE_COLORS of tools 0/6/7/9 and qc_rejected_epochs_lib).
STAGE_COLORS = {'W': '#969696', 'N1': '#9e9ac8', 'N2': '#807dba', 'N3': '#6a51a3', 'R': '#c994c7'}
# Hypnogram drawing of tools 6/7: grey step line, REM highlighted in red.
HYPNO_LINE_COLOR = '#555555'
HYPNO_REM_COLOR = '#c0392b'
# Custom (non-AASM) stages declared in config_param/custom_stages.json (written by tool 3): same helper and
# palette as tools 6/7/9 (red is reserved for REM on the hypnogram).
CUSTOM_STAGE_PALETTE = ['#8dd3c7', '#ffffb3', '#bebada', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9']


def load_custom_stages(folder):
    """Read config_param/custom_stages.json; return list of kept non-AASM labels ([] if absent/unreadable)."""
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


def non_aasm_group(label, custom_stages):
    """How an expert label outside W/N1/N2/N3/R is reported: a declared custom stage or MT under its own
    name, anything else (an unexpected label tool 3 did not correct) grouped as 'other'."""
    return label if (label in custom_stages or label == 'MT') else 'other'


def non_aasm_order(groups, custom_stages):
    """Display order of the non-AASM groups: the declared custom stages (in custom_stages.json order),
    then MT, then other."""
    def rank(g):
        if g in custom_stages:
            return (0, custom_stages.index(g))
        return (1, 0) if g == 'MT' else (2, 0)
    return sorted(dict.fromkeys(groups), key=rank)


def non_aasm_color(group, custom_stages):
    """A custom stage keeps the colour tools 6/7 give it (its rank in custom_stages.json); MT and other grey."""
    if group in custom_stages:
        return CUSTOM_STAGE_PALETTE[custom_stages.index(group) % len(CUSTOM_STAGE_PALETTE)]
    return '#bdbdbd'


# Expert hypnogram drawn over a hypnodensity (probability axis 0-1): W on top, N3 at the bottom.
OVERLAY_Y = {'W': 0.9, 'R': 0.7, 'N1': 0.5, 'N2': 0.3, 'N3': 0.1}
CONSENSUS_RULE = 'loudest_vote'   # GSSC's native rule: the most confident combination wins, per epoch
# EOG names that tool 2 may have kept in the remap: never offered as an EEG electrode.
EOG_NAME_RE = re.compile(r'(EOG|^LOC$|^ROC$|^E1$|^E2$)', re.IGNORECASE)


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


def load_json_lenient(path):
    """Strict parse, then repair one trailing comma before a closing } or ] (a common hand edit)."""
    with open(path, encoding='utf-8') as f:
        text = f.read()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return json.loads(re.sub(r',(\s*[}\]])', r'\1', text))


def list_edf(folder):
    return [f for f in sorted(Path(folder).rglob('*'))
            if f.suffix.lower() == '.edf' and not f.name.startswith('._')]


def eeg_options(cfg_entry):
    """Harmonized EEG names usable for scoring: the remap targets minus the reference and the EOG."""
    ref = cfg_entry.get('ref_channels', [])
    ref = ref if isinstance(ref, list) else []
    return [ch for ch in cfg_entry.get('remap', {}).values()
            if ch not in ref and not EOG_NAME_RE.search(ch)]


def load_eeg(edf_path, cfg_entry, wanted):
    """EEG channels `wanted` (harmonized names), re-referenced as declared in tool 2.
    Returns (raw, present) where `present` lists the wanted channels found."""
    raw = mne.io.read_raw_edf(str(edf_path), preload=False, encoding='latin-1',
                              include=list(cfg_entry.get('remap', {}).keys()), verbose=False)
    raw, _ = drop_suffix_duplicates(raw)
    raw.rename_channels(adapt_remap_dict_to_suffixes(raw, cfg_entry.get('remap', {})))
    ref = cfg_entry.get('ref_channels', [])
    present = [ch for ch in wanted if ch in raw.ch_names]
    if not present:
        return None, []
    ref_present = [ch for ch in ref if ch in raw.ch_names] if isinstance(ref, list) else []
    if ref == 'average':
        keep = [ch for ch in raw.ch_names if not EOG_NAME_RE.search(ch)]   # average of every EEG
    else:
        keep = present + [ch for ch in ref_present if ch not in present]
    raw.pick(keep)
    raw.load_data()
    raw.set_channel_types({ch: 'eeg' for ch in raw.ch_names}, verbose=False)
    # GSSC was trained on referenced derivations (e.g. C3-M2): apply tool 2's reference.
    if ref == 'average':
        raw.set_eeg_reference(ref_channels='average', verbose=False)
    elif ref_present:
        raw.set_eeg_reference(ref_channels=ref_present, verbose=False)
    raw.pick(present)
    return raw, present


def load_eog(edf_path, cfg_entry, sfreq, n_times, mode):
    """EOG from tool 2's context channels, at the EEG rate and length.
    Returns (data (n_ch, n_times) in V, names, note)."""
    ctx = cfg_entry.get('context_channels') or {}
    sides = {'EOG-L': ctx.get('eog_left'), 'EOG-R': ctx.get('eog_right')}
    sides = {k: v for k, v in sides.items() if v}
    if mode == 'None':
        return None, [], 'EOG not used (EOG mode = None)'
    if not sides:
        return None, [], 'no EOG declared in tool 2 (context channels): scored on EEG only'
    # Read the EOG apart from the EEG: one read with both would upsample everything to the fastest rate.
    eog = mne.io.read_raw_edf(str(edf_path), preload=False, encoding='latin-1',
                              include=list(sides.values()), verbose=False)
    eog, _ = drop_suffix_duplicates(eog)
    found = {k: v if v in eog.ch_names else f'{v}-0' for k, v in sides.items()}
    found = {k: v for k, v in found.items() if v in eog.ch_names}
    if not found:
        return None, [], f'declared EOG {list(sides.values())} not found in the EDF: scored on EEG only'
    eog.pick(list(found.values()))
    eog.load_data()
    if not np.isclose(eog.info['sfreq'], sfreq):
        eog.resample(sfreq, npad='auto', verbose=False)
    data = eog.get_data()
    # Same length as the EEG (a resample can be off by a sample or two): pad the end with zeros.
    if data.shape[1] >= n_times:
        data = data[:, :n_times]
    else:
        data = np.pad(data, ((0, 0), (0, n_times - data.shape[1])))
    rows = {k: data[i] for i, k in enumerate(found)}
    if mode == 'L − R bipolar' and len(rows) == 2:
        # The L−R derivation doubles the conjugate eye movements and cancels the EEG common to both.
        return (rows['EOG-L'] - rows['EOG-R'])[None, :], ['EOG L-R'], ''
    note = '' if len(rows) == 2 else f'only {list(rows)[0]} found: used alone'
    return np.vstack(list(rows.values())), list(rows), note


def gssc_combination_logits(infer, inst, eeg_names, eog_names):
    """The scoring loop of gssc 0.0.9's EEGInfer.mne_infer, step for step, keeping what mne_infer
    throws away: the network output of EVERY signal combination (mne_infer only returns the final
    stages, not the probabilities). Pinned to gssc 0.0.9: re-check this loop before upgrading gssc.
    Returns (logits (n_combinations, n_epochs, 5), combination names, epoch start times in s)."""
    # GSSC was trained on 0.3-30 Hz data: filter where the data are wider (as mne_infer does)
    band = [0.3 if round(inst.info['highpass'], 2) < 0.3 else None,
            30. if round(inst.info['lowpass'], 2) > 30. else None]
    inst.filter(*band, verbose=False)
    signals = {'eeg': {'chans': eeg_names, 'drop': True, 'flip': False},
               'eog': {'chans': eog_names, 'drop': True, 'flip': False}}
    epo, start_time = gssc_prepare_inst(inst, infer.sig_len, infer.cut)   # 2560-sample (30 s) epochs
    sig_combs, perm_matrix, all_chans, _ = gssc_permute_sigs(epo, signals)
    data = gssc_zscore(epo.get_data(picks=all_chans) * 1e6)
    device = 'cuda' if infer.use_cuda else 'cpu'
    logits, combos = [], []
    for perm in perm_matrix:
        sigs, used = {}, {'eeg': None, 'eog': None}
        for sig_idx, (sig_name, chans) in enumerate(sig_combs.items()):
            if not len(chans[perm[sig_idx]]):
                continue                                  # this signal left out of the combination
            ch, coef = gssc_check_flip_chan(chans[perm[sig_idx]][0])
            sig = coef * data[:, all_chans.index(ch), :infer.sig_len]
            sigs[sig_name] = torch.tensor(sig, dtype=torch.float32).reshape(-1, 1, sig.shape[-1]).to(device)
            used[sig_name] = ch
        with torch.no_grad():
            reps = infer.net(sigs, rep_output='rep_only').swapaxes(-1, 1)
            y, _ = infer.con_net(reps, torch.zeros(10, 1, 256).to(device))   # context network over the night
        logits.append(y[:, 0, ].float().cpu().numpy())
        combos.append(' + '.join(v for v in used.values() if v))
    n_ep = logits[0].shape[0]
    return np.array(logits), combos, np.arange(start_time, n_ep * 30, 30)


def run_gssc(infer, raw_eeg, eog_data, eog_names):
    """One GSSC scoring. Returns (stages (str array), onsets_s, probs (n_epochs, 5),
    winning combination per epoch, number of combinations)."""
    eeg_names = list(raw_eeg.ch_names)
    data = raw_eeg.get_data()
    names, types = list(eeg_names), ['eeg'] * len(eeg_names)
    if eog_data is not None:
        data = np.vstack([data, eog_data])
        names += eog_names
        types += ['eog'] * len(eog_names)
    raw_in = mne.io.RawArray(data, mne.create_info(names, raw_eeg.info['sfreq'], types), verbose=False)
    logits, combos, onsets = gssc_combination_logits(infer, raw_in, eeg_names, eog_names)
    # Consensus = gssc's own loudest_vote (the same call as in mne_infer): per epoch, the combination
    # whose output is the most confident (lowest entropy) wins.
    stage_int = np.asarray(gssc_loudest_vote(logits))
    # Which combination won: the network outputs log-probabilities, so the lowest entropy of
    # loudest_vote is the highest log-probability of the predicted stage.
    win = np.argmax(logits.max(axis=2), axis=0)
    win_logits = logits[win, np.arange(logits.shape[1])]
    e = np.exp(win_logits - win_logits.max(axis=1, keepdims=True))     # softmax -> probabilities
    probs = e / e.sum(axis=1, keepdims=True)
    if not np.array_equal(probs.argmax(axis=1), stage_int):
        raise RuntimeError('hypnodensity inconsistent with the GSSC stages (gssc version changed?)')
    stages = np.array([GSSC_STAGE[int(s)] for s in stage_int])
    return stages, onsets.astype(float), probs, np.array(combos)[win], len(combos)


def load_expert(path):
    """Expert hypnogram, one label per line (tool 3 format)."""
    with open(path, encoding='utf-8', errors='replace') as f:
        return np.array([ln.strip() for ln in f if ln.strip()])


def compare_scoring(expert, gssc):
    """Agreement on the epochs both scored, restricted to expert labels W/N1/N2/N3/R
    (MT, custom stages and '?' are left out and counted)."""
    n = min(len(expert), len(gssc))
    e, g = expert[:n], gssc[:n]
    valid = np.isin(e, STAGES)
    e, g = e[valid], g[valid]
    res = {'n_compared': int(valid.sum()), 'n_excluded_expert': int((~valid).sum())}
    if len(e) == 0:
        res.update({'accuracy': np.nan, 'kappa': np.nan, **{f'f1_{s}': np.nan for s in STAGES}})
        cm = np.zeros((5, 5), dtype=int)
    else:
        res['accuracy'] = accuracy_score(e, g)
        res['kappa'] = cohen_kappa_score(e, g, labels=STAGES)
        f1 = f1_score(e, g, labels=STAGES, average=None, zero_division=np.nan)
        res.update({f'f1_{s}': (f1[i] if s in e or s in g else np.nan) for i, s in enumerate(STAGES)})
        cm = confusion_matrix(e, g, labels=STAGES)
    cm_long = pd.DataFrame([{'expert': STAGES[i], 'gssc': STAGES[j], 'n': int(cm[i, j])}
                            for i in range(5) for j in range(5)])
    return res, cm, cm_long


def _hours(n):
    return np.arange(n + 1) * EPOCH_SEC / 3600


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


def recording_start(edf_path):
    """Recording start usable as a clock time, or None. 01.01.85 00.00.00 is the EDF+ placeholder of an
    unknown (anonymised) start date and time, not a night starting at midnight: treated as unknown."""
    dt = read_edf_start_datetime(edf_path)
    if dt is None or dt == datetime.datetime(1985, 1, 1, 0, 0, 0):
        return None
    return dt


def clock_hours(dt):
    """Clock time as hours on a noon-to-noon scale (12 -> 36): a night that crosses midnight stays
    continuous (23:30 -> 23.5, 01:00 -> 25.0)."""
    h = dt.hour + dt.minute / 60 + dt.second / 3600
    return h + 24 if h < 12 else h


def clock_label(h):
    """Hours on the noon-to-noon scale -> 'HH:MM'."""
    total_min = int(round(h * 60))
    return f'{(total_min // 60) % 24:02d}:{total_min % 60:02d}'


def clock_ticks(x_min, x_max):
    """Full-hour ticks between x_min and x_max (every hour up to 12 h shown, else every 2 h)."""
    step = 1 if x_max - x_min <= 12 else 2
    first = int(np.ceil(x_min / step) * step)
    return list(range(first, int(np.floor(x_max)) + 1, step))


def add_clock_axis(ax, start_dt, x_max_h):
    """Second x axis under a relative-time axis (hours from recording start) giving the clock time."""
    if start_dt is None:
        return
    s = clock_hours(start_dt)
    ticks = clock_ticks(s, s + x_max_h)
    # placed ~0.65 inch under the axis whatever its height (the location is a fraction of the axis height)
    ax_height_in = ax.get_position().height * ax.figure.get_figheight()
    sec = ax.secondary_xaxis(-0.65 / ax_height_in)
    sec.set_xticks([t - s for t in ticks])
    sec.set_xticklabels([clock_label(t) for t in ticks], fontsize=8)
    sec.set_xlabel(f'Clock time (recording start {start_dt:%H:%M}, EDF header)', fontsize=8)


def _hypno_rows(extra):
    """Stage -> row (0 = W at the top ... 4 = N3), then the extra non-AASM labels below N3."""
    rows = dict(HYPNO_Y)
    for i, s in enumerate(extra):
        rows[s] = 5 + i
    return rows


def _step_hypno(ax, labels, title, extra=(), custom_stages=()):
    """Hypnogram as in tools 6/7 (grey line, REM red); `extra` = non-AASM labels drawn as their own rows
    below N3, each highlighted in its custom-stage colour. REM and the extra rows at the line's own width
    (the toolkit's discreet convention)."""
    rows = _hypno_rows(extra)
    y = np.array([rows.get(s, np.nan) for s in labels], dtype=float)
    ax.step(_hours(len(y)), np.append(y, y[-1] if len(y) else np.nan), where='post',
            color=HYPNO_LINE_COLOR, lw=1.2)
    rem = y == HYPNO_Y['R']
    for i in np.where(rem)[0]:
        ax.plot(_hours(len(y))[i:i + 2], [y[i], y[i]], color=HYPNO_REM_COLOR, lw=1.2,
                solid_capstyle='butt')
    for s in extra:
        for i in np.where(y == rows[s])[0]:
            ax.plot(_hours(len(y))[i:i + 2], [y[i], y[i]], lw=1.2, solid_capstyle='butt',
                    color=non_aasm_color(s, list(custom_stages)))
    ax.set_yticks(list(rows.values()))
    ax.set_yticklabels(list(rows.keys()))
    ax.set_ylim(max(rows.values()) + 0.5, -0.5)
    ax.set_title(title, loc='left', fontsize=10)
    ax.grid(axis='y', alpha=0.3)


def plot_hypnograms(gssc, expert, file_id, custom_stages=(), start_dt=None):
    """Expert (if any) above GSSC. The expert's non-AASM labels (custom stages, MT, other) get their own rows
    below N3 on BOTH panels, so the two hypnograms keep the same vertical layout."""
    custom_stages = list(custom_stages)
    extra = []
    if expert is not None:
        expert = np.array([s if s in STAGES else non_aasm_group(s, custom_stages) for s in expert])
        extra = non_aasm_order([s for s in expert if s not in STAGES], custom_stages)
    n_rows = 2 if expert is not None else 1
    fig, axes = plt.subplots(n_rows, 1, figsize=(12, (2.0 + 0.3 * len(extra)) * n_rows + 0.4
                                                 + (0.6 if start_dt is not None else 0)),
                             sharex=True, squeeze=False)
    axes = axes[:, 0]
    if expert is not None:
        _step_hypno(axes[0], expert, 'Expert hypnogram', extra, custom_stages)
    _step_hypno(axes[-1], gssc, 'GSSC hypnogram (automatic)', extra, custom_stages)
    axes[-1].set_xlabel('Time from recording start (h)')
    add_clock_axis(axes[-1], start_dt, len(gssc) * EPOCH_SEC / 3600)
    fig.suptitle(file_id, fontsize=11)
    fig.tight_layout()
    return fig


def _overlay_expert(ax, expert, n_epochs, lw=0.9, x0=0.0):
    """Expert hypnogram as a black step line over a hypnodensity (W on top, N3 at the bottom), REM in red at
    the line's width. Epochs whose expert label is not W/N1/N2/N3/R (MT, custom, '?') leave a gap. Returns
    True if drawn."""
    if expert is None:
        return False
    labels = np.asarray(expert)[:n_epochs]
    if not np.isin(labels, STAGES).any():
        return False
    y = np.array([OVERLAY_Y.get(s, np.nan) for s in labels], dtype=float)
    x = x0 + _hours(len(y))
    ax.step(x, np.append(y, y[-1]), where='post', color='black', lw=lw)
    for i in np.where(labels == 'R')[0]:
        ax.plot(x[i:i + 2], [y[i], y[i]], color=HYPNO_REM_COLOR, lw=lw, solid_capstyle='butt')
    return True


def plot_hypnodensity(probs, file_id, expert=None, start_dt=None):
    fig, ax = plt.subplots(figsize=(12, 2.8 + (0.6 if start_dt is not None else 0)))
    x = _hours(len(probs))[:-1]
    ax.stackplot(x, probs.T, labels=STAGES, colors=[STAGE_COLORS[s] for s in STAGES],
                 step='post', linewidth=0)
    if _overlay_expert(ax, expert, len(probs)):
        ax_st = ax.twinx()                       # stage names of the expert line, on the right
        ax_st.set_ylim(0, 1)
        ax_st.set_yticks(list(OVERLAY_Y.values()))
        ax_st.set_yticklabels(list(OVERLAY_Y.keys()), fontsize=8)
        ax_st.set_ylabel('Expert (black line)', fontsize=8)
    ax.set_xlim(0, _hours(len(probs))[-1])
    ax.set_ylim(0, 1)
    ax.set_ylabel('Probability')
    ax.set_xlabel('Time from recording start (h)')
    ax.legend(loc='upper left', bbox_to_anchor=(1.06, 1.0), fontsize=9, frameon=False)
    ax.set_title(f'{file_id}: hypnodensity (probability of each stage per epoch)', loc='left', fontsize=10)
    add_clock_axis(ax, start_dt, len(probs) * EPOCH_SEC / 3600)
    fig.tight_layout()
    return fig


HYPNODENSITY_ROWS_PER_FIG = 20   # participants per stacked figure in the database report


def plot_hypnodensity_stack(nights, x_range, title, clock=False):
    """One hypnodensity per row, all on the same time axis, for a quick visual comparison across
    participants. nights = list of (row label, probs (n_epochs, 5), expert labels or None, start in
    noon-to-noon clock hours or None). clock=True: rows placed at their real clock time (x_range in
    clock hours, midnight marked), so early and late sleepers stand out; clock=False: hours from each
    recording start. The expert hypnogram, when available, is drawn over the row as a black step line
    (W on top, N3 at the bottom)."""
    n = len(nights)
    height = 0.6 * n + 1.0                       # inches: fixed margins below, whatever n
    fig, axes = plt.subplots(n, 1, figsize=(12, height), sharex=True, squeeze=False)
    axes = axes[:, 0]
    has_expert = False
    for ax, (label, probs, expert, start_h) in zip(axes, nights):
        x0 = start_h if clock else 0.0
        x = x0 + _hours(len(probs))[:-1]
        ax.stackplot(x, probs.T, colors=[STAGE_COLORS[s] for s in STAGES], step='post', linewidth=0)
        has_expert = _overlay_expert(ax, expert, len(probs), lw=0.7, x0=x0) or has_expert
        if clock and x_range[0] <= 24 <= x_range[1]:
            ax.axvline(24, color='#444444', lw=0.6, ls='--')          # midnight
        ax.set_xlim(*x_range)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_ylabel(label, rotation=0, ha='right', va='center', fontsize=8)
        for side in ('top', 'right', 'left'):
            ax.spines[side].set_visible(False)
    if clock:
        ticks = clock_ticks(*x_range)
        axes[-1].set_xticks(ticks)
        axes[-1].set_xticklabels([clock_label(t) for t in ticks])
        axes[-1].set_xlabel('Clock time (recording start from the EDF header; dashed line = midnight)')
    else:
        axes[-1].set_xlabel('Time from recording start (h)')
    handles = [plt.Rectangle((0, 0), 1, 1, color=STAGE_COLORS[s]) for s in STAGES]
    fig.legend(handles, STAGES, loc='upper right', bbox_to_anchor=(0.99, 1.0), ncol=5, fontsize=8,
               frameon=False)
    if has_expert:
        title += ' (black line = expert hypnogram, W top to N3 bottom)'
    fig.suptitle(title, x=0.01, ha='left', fontsize=10)
    fig.subplots_adjust(left=0.1, right=0.99, top=1 - 0.45 / height, bottom=0.5 / height, hspace=0.25)
    return fig


def plot_confusion(cm, title):
    """Rows = expert, columns = GSSC. Colour and bold text = % of the expert row, (n) = epoch count."""
    row_sum = cm.sum(axis=1, keepdims=True)
    pct = np.where(row_sum > 0, cm / np.where(row_sum > 0, row_sum, 1) * 100, np.nan)
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    # Full ColorBrewer Greens (white to #00441b); the text turns white on the dark cells
    greens = plt.cm.Greens
    im = ax.imshow(pct, cmap=greens, vmin=0, vmax=100)
    for i in range(5):
        for j in range(5):
            if row_sum[i, 0] > 0:
                # the % of the expert row is the main reading (comparable across stages and nights);
                # the epoch count below says how many epochs that % rests on.
                # White text where the cell colour is dark (luma of its RGB below 0.5).
                r, g, b, _ = greens(pct[i, j] / 100)
                color = 'white' if 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.5 else '#222222'
                ax.text(j, i - 0.1, f'{pct[i, j]:.0f} %', ha='center', va='center', fontsize=10,
                        fontweight='bold', color=color)
                ax.text(j, i + 0.25, f'({cm[i, j]})', ha='center', va='center', fontsize=7, color=color)
    ax.set_xticks(range(5))
    ax.set_xticklabels(STAGES)
    ax.set_yticks(range(5))
    ax.set_yticklabels(STAGES)
    ax.set_xlabel('GSSC')
    ax.set_ylabel('Expert')
    ax.set_title(title, fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.046, label='% of the expert stage (row)')
    fig.tight_layout()
    return fig


PROFILE_COLS = [f'p_{s}' for s in STAGES] + [f'pct_gssc_{s}' for s in STAGES]


def non_aasm_profile(expert, gssc, probs, custom_stages):
    """What GSSC makes of the expert epochs outside W/N1/N2/N3/R (custom stages, MT, other): per group, the
    number of epochs, the mean hypnodensity (p_W..p_R) and the % of epochs per GSSC stage. The 5 AASM expert
    stages are added as reference rows (kind 'aasm_reference'), to see which stage a custom one resembles.
    Returns an empty DataFrame when the night has no such epoch."""
    n = min(len(expert), len(gssc))
    e, g, pr = np.asarray(expert[:n]), np.asarray(gssc[:n]), probs[:n]
    if np.isin(e, STAGES).all():
        return pd.DataFrame()
    groups = np.array([s if s in STAGES else non_aasm_group(s, custom_stages) for s in e])
    order = list(STAGES) + non_aasm_order([x for x in groups if x not in STAGES], list(custom_stages))
    rows = []
    for st in order:
        m = groups == st
        if not m.any():
            continue
        kind = ('aasm_reference' if st in STAGES else 'custom' if st in custom_stages
                else 'MT' if st == 'MT' else 'other')
        row = {'stage_expert': st, 'kind': kind, 'n_epochs': int(m.sum())}
        for i, s in enumerate(STAGES):
            row[f'p_{s}'] = float(pr[m, i].mean())
        for s in STAGES:
            row[f'pct_gssc_{s}'] = float((g[m] == s).mean() * 100)
        if kind == 'other':
            row['labels'] = ', '.join(sorted(set(e[m])))
        rows.append(row)
    return pd.DataFrame(rows)


def pool_profiles(df, custom_stages=()):
    """Database profile: per expert stage, the per-night means weighted by the number of epochs."""
    rows = []
    for st, d in df.groupby('stage_expert', sort=False):
        w = d['n_epochs'].astype(float)
        row = {'stage_expert': st, 'kind': d['kind'].iloc[0], 'n_epochs': int(w.sum()),
               'n_nights': int(d['file_id'].nunique())}
        for c in PROFILE_COLS:
            row[c] = float((d[c].astype(float) * w).sum() / w.sum())
        rows.append(row)
    out = pd.DataFrame(rows)
    order = list(STAGES) + non_aasm_order([s for s in out['stage_expert'] if s not in STAGES],
                                          list(custom_stages))
    out['_o'] = [order.index(s) for s in out['stage_expert']]
    return out.sort_values('_o').drop(columns='_o').reset_index(drop=True)


def plot_non_aasm_profile(df, title):
    """Two panels of horizontal stacked bars, one bar per expert stage: the mean hypnodensity (left) and
    the % of epochs per GSSC stage (right). AASM reference bars on top, the non-AASM stages below a line."""
    labels = [f'{r.stage_expert}  (n={r.n_epochs})' for r in df.itertuples()]
    y = np.arange(len(df))
    fig, axes = plt.subplots(1, 2, figsize=(12, 0.45 * len(df) + 1.6), sharey=True)
    for ax, prefix, scale, sub in [(axes[0], 'p_', 1.0, 'Mean hypnodensity (probability)'),
                                   (axes[1], 'pct_gssc_', 100.0, 'Stage given by GSSC (% of epochs)')]:
        left = np.zeros(len(df))
        for s in STAGES:
            v = df[f'{prefix}{s}'].astype(float).values / scale
            ax.barh(y, v, left=left, color=STAGE_COLORS[s], label=s, height=0.7)
            for yi, (l0, vi) in enumerate(zip(left, v)):
                if vi >= 0.12:
                    ax.text(l0 + vi / 2, yi, f'{vi * 100:.0f}%', ha='center', va='center', fontsize=7,
                            color='#222222')
            left += v
        n_ref = int((df['kind'] == 'aasm_reference').sum())
        if 0 < n_ref < len(df):
            ax.axhline(n_ref - 0.5, color='#222222', lw=1)
        ax.set_xlim(0, 1)
        ax.set_title(sub, fontsize=9)
        ax.set_xticks([0, 0.5, 1])
        ax.set_xticklabels(['0', '50 %', '100 %'])
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(labels, fontsize=8)
    axes[0].invert_yaxis()
    axes[0].set_ylabel('Expert stage')
    axes[1].legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=8, frameon=False)
    fig.suptitle(title, x=0.01, ha='left', fontsize=10)
    fig.tight_layout()
    return fig


NON_AASM_HELP = ('<p style="font-size:85%;color:#555">GSSC only produces W, N1, N2, N3 and R, so the expert '
                 'epochs scored with another label (custom stages declared in tool 3, MT, other) cannot enter '
                 'the kappa, accuracy, F1 or confusion matrix. This section shows what GSSC makes of them: '
                 'their mean hypnodensity and the stage GSSC gives them, next to the same bars for the five '
                 'AASM expert stages (above the line) as a reference. Example: a custom N4 whose bar looks '
                 'like the N3 bar is read as deep N3 sleep by GSSC.</p>')


def metrics_html(res):
    def f(v, pct=False):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return '–'
        return f'{v * 100:.1f} %' if pct else f'{v:.3f}'
    rows = [('Epochs compared', res['n_compared']),
            ('Epochs left out (expert label not W/N1/N2/N3/R)', res['n_excluded_expert']),
            ('Accuracy', f(res['accuracy'], pct=True)), ("Cohen's kappa", f(res['kappa']))]
    rows += [(f'F1 {s}', f(res[f'f1_{s}'])) for s in STAGES]
    body = ''.join(f'<tr><td style="padding:2px 14px 2px 0">{k}</td><td><b>{v}</b></td></tr>'
                   for k, v in rows)
    return (f'<table style="font-size:90%">{body}</table>'
            '<p style="font-size:85%;color:#555">Kappa: agreement corrected for chance (above 0.6 = '
            'substantial, above 0.8 = almost perfect; two human scorers typically reach 0.7-0.8). '
            'F1 per stage: balance of the epochs of that stage found and correctly labelled.</p>')


def params_html(params):
    rows = ''.join(f'<tr><td style="padding:2px 12px 2px 0"><b>{k}</b></td><td>{v}</td></tr>'
                   for k, v in params.items())
    return f'<table style="font-size:90%">{rows}</table>'


def save_report_html(path, title, items):
    """mne.Report from an ORDERED list of ('fig'|'html', title, payload); None payloads skipped."""
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
    for kind, _, payload in items:
        if kind == 'fig' and payload is not None:
            plt.close(payload)


def score_night(infer, edf_path, cfg_entry, wanted, one_mode, eog_mode):
    """The scoring part of tool 3bis's score_one, unchanged: EEG `wanted` re-referenced as declared in tool 2,
    EOG as declared (EOG mode), then GSSC. One-electrode mode never substitutes a missing electrode.
    Raises on a fatal problem. Returns (stages, onsets, probs, winner, n_comb, present, eog_names, notes)."""
    raw, present = load_eeg(edf_path, cfg_entry, wanted)
    missing = [ch for ch in wanted if ch not in present]
    if one_mode and missing:
        raise RuntimeError(f'electrode {wanted[0]} absent (not substituted)')
    if not present:
        raise RuntimeError(f'none of the chosen electrodes present ({", ".join(wanted)})')
    notes = []
    if missing:
        notes.append(f'electrode(s) absent, scored without them: {", ".join(missing)}')
    eog_data, eog_names, eog_note = load_eog(edf_path, cfg_entry, raw.info['sfreq'], raw.n_times, eog_mode)
    if eog_note:
        notes.append(eog_note)
    stages, onsets, probs, winner, n_comb = run_gssc(infer, raw, eog_data, eog_names)
    return stages, onsets, probs, winner, n_comb, present, eog_names, notes
