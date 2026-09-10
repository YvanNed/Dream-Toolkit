"""
Shared analysis + report helpers for Tool 7 — Interactive QC of rejected epochs (Phase 2b).

Consumed by:
  - 7_reject_manually_voila.ipynb  (interactive, per-participant QC + manual override)
  - 7_reject_manually_batch.py      (database-level Section-2 report, no interaction)

Reads the outputs written by 6_preprocessing_voila for one participant:
  - {file_id}_all-epo.fif                    all epochs + per-epoch rejection metadata (MNE)
  - {file_id}_preprocessing_params.json      per-stage thresholds + preprocessing params (optional)

Design notes
------------
* The per-EPOCH rejection decision is authoritative from epochs.metadata (reject_flag,
  reject_method, flag_<method>). Tool 6 does NOT persist a per-(epoch, channel) mask, so the
  per-CHANNEL attribution shown here (which channel drove the flag, margins) is RECOMPUTED from
  the signal with the same formulas + the persisted thresholds. This stays faithful because the
  formulas and thresholds are shared; only the display attribution is reconstructed.
* EEG only: whatever channels the .fif contains are used as-is (that is exactly the channel set
  tool 6 flagged over). The raw EDF is never reloaded. The optional EOG/EMG/ECG "context" traces
  shown under the per-epoch montage are NOT read from the raw EDF either: tool 6 persists them as a
  {file_id}_context-epo.fif companion (epoched identically), which this module loads on demand.
* Kept in sync with 6_preprocessing_voila: METHOD_ORDER / colours / custom-stage helpers /
  Welch-PSD config / specparam 1/f fit.
"""
import os
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

import mne
from scipy.signal import spectrogram as sp_spectrogram

try:
    from specparam import SpectralModel
    HAS_SPECPARAM = True
except Exception:
    HAS_SPECPARAM = False

# Optional PSD smoothing before the 1/f fit (honours tool 6's psd_smoothing sidecar block).
try:
    from statsmodels.nonparametric.smoothers_lowess import lowess
    HAS_LOWESS = True
except Exception:
    HAS_LOWESS = False
    lowess = None

mne.set_log_level('ERROR')

# ---- Rejection-method registry (single source of truth, kept in sync with tool 6) ----
METHOD_ORDER  = ['amplitude', 'gradient', 'flat', '1f_r2', '1f_error', 'event']
METHOD_CODE   = {m: i + 1 for i, m in enumerate(METHOD_ORDER)}   # 1..6  (0 = none)
MULTIPLE_CODE = len(METHOD_ORDER) + 1                            # 7 = multiple
METHOD_LABEL  = {'amplitude': 'Amplitude', 'gradient': 'Gradient', 'flat': 'Flat',
                 '1f_r2': '1/f R²', '1f_error': '1/f error', 'event': 'Event'}
# Heatmap / per-method colours (index = method code; SINGLE SOURCE, duplicated verbatim in tool 6's
# plot_rejection_heatmap and tools 8/8-voila — keep in sync). CVD-validated qualitative palette
# (dataviz skill: 6 method hues chosen to maximise colour-blind separation, worst adjacent protan/deutan
# ΔE ~7; index 0 'none' stays dark so heatmap marks pop; 'multiple' = cyan, distinct from all 6 methods
# and from black). Identity is never colour-alone: every plot carries a legend / colour-coded labels.
HEATMAP_COLORS = ['#1c0a3b', '#e34948', '#4a3aa7', '#2a78d6', '#008300', '#eda100', '#e87ba4', '#22d3d3']
HEATMAP_LABELS = ['none', 'amplitude', 'gradient', 'flat', '1/f R2', '1/f error', 'event', 'multiple']
# Per-method overlay colour keyed by METHOD_ORDER label.
METHOD_COLOR   = {m: HEATMAP_COLORS[METHOD_CODE[m]] for m in METHOD_ORDER}

# Per-epoch montage display scales: µV peak-to-peak allocated to ONE channel row, per channel type.
# FIXED (never auto-scaled to the displayed window) so amplitudes stay comparable from epoch to epoch —
# a 75 µV slow wave always draws the same height, which is what makes visual scoring criteria usable.
# Values follow clinical PSG display conventions (AASM ~7 µV/mm on a ~2 cm row for EEG/EOG, high-gain
# chin EMG, 1 mV for ECG); traces larger than their row overflow into the neighbouring one, exactly as in
# a clinical viewer — the amplitude is never clipped or hidden. Editable per run in the tool-7 navigator.
DISPLAY_SCALE_UV = {'EEG': 150.0, 'EOG': 300.0, 'EMG': 100.0, 'ECG': 1000.0}

# Context-channel (EOG/EMG/ECG) trace colours for the per-epoch montage (role labels written by tool 6).
# Kept OUTSIDE the six method hues (teal / olive / sienna) so a context trace never impersonates a
# rejection-method colour in the montage.
CTX_COLOR = {'EOG-L': '#0a8f8f', 'EOG-R': '#0a8f8f', 'EMG': '#8a6d1f', 'ECG': '#a0522d'}

# ---- High-density montage support (32-64 channel Curry / HD-EEG montages) ----
# Everything below adapts to the CHANNEL COUNT, never to the file format: a dense EDF montage gets the
# same treatment and a sparse Curry one keeps the classic layout. At or below HD_CHANNEL_THRESHOLD the
# rendering is IDENTICAL to the pre-high-density version, so the 3-6 channel PSG tools are unchanged.
#
# Why this is needed at all: with 3 channels "reject the epoch if ANY channel is flagged" is sound; with
# 32 it is not. Measured on a real 32-channel Curry night (o_S007, 1314 epochs): 72.6% of the epochs are
# rejected, and 451 of them (34% of the file) are flagged on a SINGLE electrode out of 32 - where the
# right action is to drop that electrode, not to throw away 30 s of 32-channel EEG. Hence the channel
# triage helpers (channel_badness / channels_over_threshold / recompute_reject's epoch rule) below.
HD_CHANNEL_THRESHOLD = 12       # channels above which the high-density layout/behaviour kicks in
MONTAGE_ROW_IN       = 0.50     # inches per channel row (classic)
MONTAGE_ROW_MIN_IN   = 0.16     # floor: thinner rows are unreadable anyway
MONTAGE_MAX_H_IN     = 13.0     # cap on the montage figure height (inches)
N_1F_SUBSAMPLE_HD    = 250      # default 1/f epoch subsample above the threshold (0 = fit every epoch)
DEFAULT_TABLE_ROWS   = 8        # per-channel metric table: flagged channels + this many worst ones

AASM_STAGES = ['W', 'N1', 'N2', 'N3', 'R']

# ---- Custom (non-AASM) sleep stages (duplicated from tool 6, kept in sync) ----
BASE_STAGE_COLORS   = {'W': '#969696', 'N1': '#9e9ac8', 'N2': '#807dba', 'N3': '#6a51a3', 'R': '#c994c7'}
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


def parse_custom_field(text):
    """Parse the comma-separated 'Custom stages' field into a clean, de-duplicated list."""
    seen = []
    for tok in str(text).split(','):
        tok = tok.strip()
        if tok and tok not in seen:
            seen.append(tok)
    return seen


def custom_stage_style(custom_stages):
    """Stage -> y-position, stage -> colour, and axis ticks for AASM + custom stages.
    AASM keeps YASA heights (W top ... N3 bottom); custom stages stack below N3 (-1, -2, ...).
    Returns (stage_y, stage_colors, ytick_pos, ytick_labels)."""
    stage_y = {'W': 4, 'R': 3, 'N1': 2, 'N2': 1, 'N3': 0}
    stage_colors = dict(BASE_STAGE_COLORS)
    for i, cs in enumerate(custom_stages):
        stage_y[cs] = -1 - i
        stage_colors[cs] = CUSTOM_STAGE_PALETTE[i % len(CUSTOM_STAGE_PALETTE)]
    ordered = sorted(stage_y.items(), key=lambda kv: kv[1])
    ytick_pos = [y for _, y in ordered]
    ytick_labels = ['REM' if s == 'R' else s for s, _ in ordered]
    return stage_y, stage_colors, ytick_pos, ytick_labels


# ---- Default thresholds (tool-6 defaults; fallback when the params JSON is absent) ----
DEFAULT_THRESHOLDS = {
    'amplitude_ptp_uV': {'W': 300.0, 'N1': 250.0, 'N2': 200.0, 'N3': 200.0, 'R': 250.0},
    'flat_ptp_uV': 1.0,
    'gradient_uV_per_sample': 100.0,
    '1f_mae_max': 0.15,
    '1f_r2_min': 0.95,
}

# Frequency bands for the per-epoch band-power bars (Hz).
BANDS = [('delta', 0.5, 4), ('theta', 4, 8), ('alpha', 8, 12), ('sigma', 12, 16), ('beta', 16, 30)]


# ---------------------------------------------------------------------------
# Discovery + loading
# ---------------------------------------------------------------------------
def find_participants(deriv_root):
    """List participants under a derivatives root: every *_all-epo.fif found recursively.
    Returns a list of dicts {file_id, fif, folder} sorted by file_id."""
    root = Path(deriv_root)
    out = []
    if not root.exists():
        return out
    suffix = '_all-epo.fif'
    for fif in root.rglob('*' + suffix):
        fid = fif.name[:-len(suffix)]
        out.append({'file_id': fid, 'fif': fif, 'folder': fif.parent})
    out.sort(key=lambda d: d['file_id'])
    return out


def load_params(folder, file_id, custom_stages_fallback=()):
    """Read {file_id}_preprocessing_params.json (thresholds + preprocessing params).
    Missing / unreadable -> DEFAULT_THRESHOLDS and found=False (non-fatal)."""
    path = Path(folder) / f'{file_id}_preprocessing_params.json'
    thresholds = {
        'amplitude_ptp_uV': dict(DEFAULT_THRESHOLDS['amplitude_ptp_uV']),
        'flat_ptp_uV': DEFAULT_THRESHOLDS['flat_ptp_uV'],
        'gradient_uV_per_sample': DEFAULT_THRESHOLDS['gradient_uV_per_sample'],
        '1f_mae_max': DEFAULT_THRESHOLDS['1f_mae_max'],
        '1f_r2_min': DEFAULT_THRESHOLDS['1f_r2_min'],
    }
    info = {'found': False, 'custom_stages': list(custom_stages_fallback),
            'resample': None, 'filter': None, 'methods_run': None,
            'fit_range': (2.0, 45.0),   # 1/f fit window (Hz); tool-6 default when absent
            'epoch_length_s': 30,       # scoring epoch length (s); classic 30 s when absent
            'psd_smoothing': None}      # tool-6 PSD-smoothing block; None when absent (older sidecars)
    if not path.exists():
        return thresholds, info
    try:
        with open(path, encoding='utf-8') as f:
            data = json.load(f)
        rt = data.get('rejection_thresholds', {})
        if isinstance(rt.get('amplitude_ptp_uV'), dict):
            thresholds['amplitude_ptp_uV'].update({k: float(v) for k, v in rt['amplitude_ptp_uV'].items()})
        for k in ('flat_ptp_uV', 'gradient_uV_per_sample', '1f_mae_max', '1f_r2_min'):
            if k in rt:
                thresholds[k] = float(rt[k])
        fr = rt.get('1f_fit_range_hz')
        if isinstance(fr, (list, tuple)) and len(fr) == 2:
            info['fit_range'] = (float(fr[0]), float(fr[1]))
        info['found'] = True
        info['custom_stages'] = list(data.get('custom_stages', custom_stages_fallback))
        info['resample'] = data.get('resample')
        info['filter'] = data.get('filter')
        info['methods_run'] = data.get('methods_run')
        info['epoch_length_s'] = int(data.get('epoch_length_s', 30))
        info['psd_smoothing'] = data.get('psd_smoothing')   # honoured by compute_psds (tool-7 plots)
    except Exception:
        pass
    return thresholds, info


def load_participant(fif_path):
    """Read one {file_id}_all-epo.fif into a dict with signal + metadata.
    Returns dict: epochs, data_uV (n_ep, n_ch, n_t), sfreq, ch_names, meta (DataFrame),
    stages (array of str), reject_flag (bool array), methods_present (list)."""
    epochs = mne.read_epochs(str(fif_path), preload=True, verbose=False)
    meta = epochs.metadata.reset_index(drop=True).copy()
    data_uV = epochs.get_data() * 1e6                       # (n_ep, n_ch, n_t), µV
    stages = meta['stage'].astype(str).values
    reject_flag = meta['reject_flag'].astype(bool).values
    methods_present = [m for m in METHOD_ORDER
                       if 'flag_' + m in meta.columns]
    return {
        'epochs': epochs, 'data_uV': data_uV, 'sfreq': float(epochs.info['sfreq']),
        'ch_names': list(epochs.ch_names), 'meta': meta, 'stages': stages,
        'reject_flag': reject_flag, 'methods_present': methods_present,
    }


def load_context_epochs(folder, file_id):
    """Load the optional {file_id}_context-epo.fif companion written by tool 6 — the EOG/EMG/ECG
    "context" channels epoched identically to the EEG (same epoch length, same count), so they align 1:1
    with the EEG epochs by index. Returns {'data_uV' (n_ep, n_ch, n_t) µV, 'sfreq', 'labels'} or
    None if the companion is absent/unreadable (non-fatal: the per-epoch view just omits context)."""
    path = Path(folder) / f'{file_id}_context-epo.fif'
    if not path.exists():
        return None
    try:
        ep = mne.read_epochs(str(path), preload=True, verbose=False)
        return {'data_uV': ep.get_data() * 1e6, 'sfreq': float(ep.info['sfreq']),
                'labels': list(ep.ch_names)}
    except Exception:
        return None


def load_event_onsets(folder, file_id):
    """Load the optional {file_id}_event_onsets.tsv written by tool 6 — one row per scored event
    (columns type, onset_s, duration_s; onset in seconds from recording start). Returns a DataFrame or
    None if absent/unreadable (non-fatal: the montage just omits event markers). Lets tool 7 draw event
    onset lines on the per-epoch montage WITHOUT reloading the raw recording."""
    path = Path(folder) / f'{file_id}_event_onsets.tsv'
    if not path.exists():
        return None
    try:
        df = pd.read_csv(path, sep='\t')
        if not {'type', 'onset_s'}.issubset(df.columns):
            return None
        return df
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Spectral analysis (Welch PSD + specparam 1/f fit) — same config as tool 6
# ---------------------------------------------------------------------------
# --- PSD smoothing helpers (copied VERBATIM from 6_preprocessing_voila — keep in sync) ---------
def smooth_psd_median(psds_uV2, freqs, span_hz=3.0):
    """Running-median smooth of each (epoch, channel) linear PSD along frequency.
    Reproduces oscip.smooth_spectrum_median (MATLAB movmedian over a span given in Hz): a centred
    median filter that removes narrow spikes (residual line-noise harmonics, single-bin artefacts)
    BEFORE the LOWESS mean smoothing. The Hz span is converted to a point count
    (span_pts = round(span_hz / freq_res)); pandas rolling(center=True, min_periods=1) matches
    MATLAB's truncated-window edge behaviour. Non-fatal: a failing (epoch, channel) is left as-is."""
    freqs = np.asarray(freqs)
    if len(freqs) < 4:
        return psds_uV2
    freq_res = float(np.median(np.diff(freqs)))            # 0.25 Hz for a 4 s Welch window
    span_pts = max(3, int(round(span_hz / freq_res)))      # 3 Hz span -> ~12 points at 0.25 Hz
    out = np.array(psds_uV2, dtype=float, copy=True)
    n_ep, n_ch, _ = out.shape
    for ei in range(n_ep):
        for ci in range(n_ch):
            try:
                out[ei, ci] = pd.Series(out[ei, ci]).rolling(
                    window=span_pts, center=True, min_periods=1).median().to_numpy()
            except Exception:
                pass                                       # keep the raw spectrum for this (ep, ch)
    return out


def smooth_psd_lowess(psds_uV2, freqs, span_hz=2.0):
    """LOWESS-smooth each (epoch, channel) linear PSD along frequency.
    Reproduces oscip.smooth_spectrum from the Snipes MATLAB pipeline: MATLAB's non-robust
    'lowess' (local linear regression, tricube weights) over a span given in Hz, applied to the
    LINEAR power (matching oscip; specparam logs it internally afterwards). The Hz span is converted
    to a fraction of points because statsmodels expresses the window as a fraction (frac), whereas
    MATLAB uses a point count. it=0 = non-robust (matches 'lowess', not 'rlowess').
    Returns a smoothed copy with the same shape (n_ep, n_ch, n_freqs). Non-fatal: on failure a given
    (epoch, channel) spectrum is left unsmoothed so the 1/f fit still runs."""
    freqs = np.asarray(freqs)
    if len(freqs) < 4:
        return psds_uV2
    freq_res = float(np.median(np.diff(freqs)))            # 0.25 Hz for a 4 s Welch window
    span_pts = max(3, int(round(span_hz / freq_res)))      # 2 Hz span -> 8 points at 0.25 Hz
    frac = min(1.0, span_pts / len(freqs))                 # statsmodels wants a fraction of points
    out = np.array(psds_uV2, dtype=float, copy=True)
    n_ep, n_ch, _ = out.shape
    for ei in range(n_ep):
        for ci in range(n_ch):
            try:
                out[ei, ci] = lowess(out[ei, ci], freqs, frac=frac, it=0, return_sorted=False)
            except Exception:
                pass                                       # keep the raw spectrum for this (ep, ch)
    return out


def compute_psds(epochs, fmax=45.0, smoothing=None):
    """Welch PSD per epoch/channel, identical config to 6_preprocessing_voila.
    `fmax` is the 1/f fit upper bound (default 2–45 Hz range); the PSD spans up to it.
    `smoothing` is the tool-6 `psd_smoothing` sidecar block
    ({enabled, method, median_span_hz, lowess_span_hz}); when enabled the PSD is smoothed
    (median then LOWESS, same order/params as tool 6) BEFORE return, so tool 7's recomputed
    plots and per-channel 1/f attribution match what tool 6 flagged. Default None -> no smoothing
    (byte-identical to before). Returns (freqs, psds_uV2), psds shape (n_ep, n_ch, n_freqs), µV²/Hz."""
    sf = float(epochs.info['sfreq'])
    epoch_len_s = len(epochs.times) / sf                  # supports non-30 s epochs
    n_per_seg = int(min(4, epoch_len_s) * sf)              # 4 s Welch window, capped to the epoch length
    n_overlap = int(n_per_seg / 2)
    fmax_psd = min(fmax, sf / 2 - 0.5)
    psds_obj = epochs.compute_psd(method='welch', fmin=0.5, fmax=fmax_psd,
                                  n_fft=n_per_seg, n_overlap=n_overlap, n_per_seg=n_per_seg,
                                  window='hann', verbose=False)   # match tool 6 (MNE defaults to 'hamming')
    freqs, psds_uV2 = psds_obj.freqs, psds_obj.get_data() * 1e12
    if smoothing and smoothing.get('enabled'):
        try:
            if smoothing.get('method') == 'median+lowess' and smoothing.get('median_span_hz'):
                psds_uV2 = smooth_psd_median(psds_uV2, freqs, span_hz=float(smoothing['median_span_hz']))
            psds_uV2 = smooth_psd_lowess(psds_uV2, freqs, span_hz=float(smoothing.get('lowess_span_hz', 2.0)))
        except Exception:
            pass   # non-fatal: fall back to the raw PSD if smoothing fails
    return freqs, psds_uV2


def fit_1f(freqs, psd_uV2, fmin=2.0):
    """specparam aperiodic fit (fixed mode, freqs >= fmin) — same model as tool 6.
    `fmin` is the 1/f fit lower bound (default 2 Hz); the upper bound is set by `compute_psds`.
    Returns (mae, r2, offset, exponent, model) or (nan, nan, nan, nan, None) on failure."""
    if not HAS_SPECPARAM:
        return np.nan, np.nan, np.nan, np.nan, None
    fmask = freqs >= fmin
    try:
        sm = SpectralModel(peak_width_limits=[0.5, 20], aperiodic_mode='fixed',
                           min_peak_height=0.3, max_n_peaks=8, verbose=False)
        sm.fit(freqs[fmask], np.asarray(psd_uV2)[fmask])
        mae = float(sm.get_metrics('error', 'mae'))
        r2 = float(sm.get_metrics('gof', 'squared'))
        ap = sm.get_params('aperiodic')   # specparam 2.x: [offset, exponent] in 'fixed' mode
        off, exp = float(ap[0]), float(ap[1])
        return mae, r2, off, exp, sm
    except Exception:
        return np.nan, np.nan, np.nan, np.nan, None


def band_powers(freqs, psd_uV2):
    """Integrated band power (µV²) per BANDS entry for one channel PSD."""
    out = {}
    freqs = np.asarray(freqs)
    psd_uV2 = np.asarray(psd_uV2)
    for name, lo, hi in BANDS:
        m = (freqs >= lo) & (freqs < hi)
        out[name] = float(np.trapz(psd_uV2[m], freqs[m])) if m.sum() > 1 else 0.0
    return out


def _worst_1f_of_epoch(freqs, psds_ch, fmin=2.0):
    """Worst-channel 1/f fit of ONE epoch (max MAE, min R² across channels) -> (mae, r2).
    Module level so joblib can pickle it; the serial and parallel paths call this same function, so
    they are numerically identical."""
    worst_mae, best_r2 = np.nan, np.nan
    for ci in range(np.asarray(psds_ch).shape[0]):
        m, r, _o, _e, _sm = fit_1f(freqs, psds_ch[ci], fmin=fmin)
        if not np.isnan(m):
            worst_mae = m if np.isnan(worst_mae) else max(worst_mae, m)
            best_r2 = r if np.isnan(best_r2) else min(best_r2, r)
    return worst_mae, best_r2


def compute_epoch_metrics(data_uV, freqs, psds_uV2, progress=None, fmin=2.0,
                          fit_mask=None, do_1f=True, subsample_n=0, n_jobs=1, random_state=0):
    """Per-epoch rejection-metric VALUES (for the report plots — the flags themselves come from
    the .fif metadata). Peak-to-peak and gradient are per (epoch, channel); the per-epoch value
    is the worst channel (max), matching 'any channel flags the epoch'. 1/f MAE/R² are fitted per
    (epoch, channel); the per-epoch value is the worst channel (max MAE, min R²).

    `progress(done, total)` is called periodically so a UI can show advancement.
    `fit_mask` (optional bool (n_ep,)): fit 1/f ONLY on those epochs (others left NaN) — tool 7 passes
    the in-scope-stage mask so a stage sub-selection skips the slow fit on out-of-scope epochs.
    `do_1f=False` skips the 1/f fit entirely (tool 7 when no 1/f method is selected). ptp/gradient are
    always computed (cheap, vectorised).

    HIGH-DENSITY levers (cost is per epoch x channel: a 32-channel night is ~42 000 fits, ~8 min):
    `subsample_n` (int): fit only this many RANDOMLY CHOSEN in-scope epochs (0 = every one, the classic
        default). Only the METRIC DISTRIBUTIONS consume these values, and they are statistically
        identical on a subsample; the keep/reject flags always come from the .fif metadata, so a
        subsample can never change a decision. The caller must state it in the figure/UI.
    `n_jobs` (int): joblib parallelism over epochs (1 = serial, -1 = all cores). EXACT — same numbers,
        just faster. Falls back to the serial loop if joblib is unavailable or fails.

    Defaults (fit_mask=None, do_1f=True, subsample_n=0, n_jobs=1) are byte-identical to before.
    Returns dict of (n_ep,) arrays: ptp, gradient, mae, r2, (n_ep, n_ch) arrays ptp_ch, grad_ch, plus
    'n_fitted' / 'n_candidates' (how many epochs the 1/f fit actually ran on).
    """
    n_ep, n_ch, _ = data_uV.shape
    ptp_ch = np.ptp(data_uV, axis=-1)                              # (n_ep, n_ch)
    grad_ch = np.max(np.abs(np.diff(data_uV, axis=-1)), axis=-1)   # (n_ep, n_ch)
    mae = np.full(n_ep, np.nan)
    r2 = np.full(n_ep, np.nan)
    n_candidates = 0
    if do_1f and HAS_SPECPARAM and psds_uV2 is not None:
        idxs = list(range(n_ep)) if fit_mask is None \
            else list(np.where(np.asarray(fit_mask, dtype=bool))[0])
        n_candidates = len(idxs)
        if subsample_n and n_candidates > int(subsample_n):
            rng = np.random.default_rng(random_state)   # seeded: the same run gives the same figure
            idxs = sorted(int(i) for i in rng.choice(idxs, size=int(subsample_n), replace=False))
        total = len(idxs)
        par = None
        if n_jobs != 1 and total:
            try:
                from joblib import Parallel, delayed
                par = Parallel(n_jobs=n_jobs)
            except Exception:
                par = None                              # joblib missing -> serial, same numbers
        if par is not None:
            try:
                chunk = max(1, min(50, total))          # chunked so the progress bar still advances
                done = 0
                for start in range(0, total, chunk):
                    part = idxs[start:start + chunk]
                    res = par(delayed(_worst_1f_of_epoch)(freqs, psds_uV2[e], fmin) for e in part)
                    for e, (wm, br) in zip(part, res):
                        mae[e], r2[e] = wm, br
                    done += len(part)
                    if progress is not None:
                        progress(done, total)
            except Exception:
                par = None                              # any parallel failure -> redo serially
        if par is None:
            for done, ei in enumerate(idxs, 1):
                mae[ei], r2[ei] = _worst_1f_of_epoch(freqs, psds_uV2[ei], fmin=fmin)
                if progress is not None and (done % 25 == 0 or done == total):
                    progress(done, total)
    return {'ptp': ptp_ch.max(axis=1), 'gradient': grad_ch.max(axis=1),
            'mae': mae, 'r2': r2, 'ptp_ch': ptp_ch, 'grad_ch': grad_ch,
            'n_fitted': int(np.isfinite(mae).sum()), 'n_candidates': int(n_candidates)}


def ordered_present_stages(stages, custom_stages):
    """AASM order then custom stages, keeping only stages actually present."""
    present = set(np.asarray(stages).tolist())
    order = AASM_STAGES + [s for s in custom_stages if s not in AASM_STAGES]
    ordered = [s for s in order if s in present]
    ordered += [s for s in present if s not in ordered]   # any unexpected label, last
    return ordered


# ---------------------------------------------------------------------------
# Channel triage (high-density montages) — per-(epoch, channel) flags
# ---------------------------------------------------------------------------
def load_channel_flags(reports_folder, file_id, ch_names, n_epochs):
    """Read tool 6's `{file_id}_epoch_channel_rejection.tsv` into per-method (n_ep, n_ch) bool arrays.

    This is the ONLY thing tool 7 reads from `reports_preprocessing/`, and only when the user points at
    that folder: it is what makes per-channel badness exact (it carries the per-channel 1/f flags, which
    tool 7 would otherwise have to refit at ~8 min on a 32-channel night). Absent/unreadable -> None,
    and the caller falls back to `recompute_channel_flags` (time-domain only). Non-fatal by design.
    Returns {method: (n_ep, n_ch) bool} for the methods present, or None."""
    path = Path(reports_folder) / f'{file_id}_epoch_channel_rejection.tsv'
    if not path.exists():
        return None
    try:
        df = pd.read_csv(path, sep='\t')
        if not {'epoch_idx', 'channel'}.issubset(df.columns):
            return None
        out = {}
        for m in METHOD_ORDER:
            col = 'flag_' + m
            if col not in df.columns:
                continue                       # 'event' is epoch-level: never a per-channel column
            piv = df.pivot(index='epoch_idx', columns='channel', values=col)
            piv = piv.reindex(index=range(n_epochs), columns=list(ch_names))
            out[m] = piv.fillna(0).values.astype(bool)
        return out or None
    except Exception:
        return None


def recompute_channel_flags(P, thresholds):
    """Fallback for `load_channel_flags`: per-(epoch, channel) TIME-DOMAIN flags recomputed from the
    signal with the tool-6 formulas (amplitude / gradient / flat). Covers ~89 % of the flagged pairs on
    a real 32-channel night but NOT the 1/f ones, so the caller must say so in the UI.
    Returns {method: (n_ep, n_ch) bool} for amplitude / gradient / flat."""
    data = P['data_uV']
    stages = np.asarray(P['stages'])
    ptp = np.ptp(data, axis=-1)                                   # (n_ep, n_ch)
    grad = np.max(np.abs(np.diff(data, axis=-1)), axis=-1)        # (n_ep, n_ch)
    amp_thr = np.array([thresholds['amplitude_ptp_uV'].get(s, 250.0) for s in stages])[:, None]
    return {'amplitude': ptp > amp_thr,
            'gradient': grad > thresholds['gradient_uV_per_sample'],
            'flat': ptp < thresholds['flat_ptp_uV']}


def channel_badness(pair_flags, methods_sel, in_scope, n_ch):
    """Per-channel badness = % of IN-SCOPE epochs where the channel is flagged by a selected method
    (the same definition as tool 7bis's channel-first step, so both tools rank channels identically).
    Returns (overall (n_ch,) float %, {method: (n_ch,) float %})."""
    insc = np.asarray(in_scope, dtype=bool)
    n_in = int(insc.sum())
    per_method, any_flag = {}, np.zeros((int(insc.sum()), n_ch), dtype=bool) if n_in else None
    for m in METHOD_ORDER:
        if m not in pair_flags or m not in methods_sel:
            continue
        sub = np.asarray(pair_flags[m], dtype=bool)[insc]
        per_method[m] = 100.0 * sub.mean(axis=0) if n_in else np.zeros(n_ch)
        if n_in:
            any_flag |= sub
    overall = 100.0 * any_flag.mean(axis=0) if n_in else np.zeros(n_ch)
    return overall, per_method


def channels_over_threshold(badness_pct, ch_names, threshold_pct):
    """Channel names whose badness exceeds `threshold_pct` (tool 7bis's channel-first rule). Tool 7
    applies it as a SUGGESTION: it unticks those channels, and the user can re-tick any of them by hand.
    Careful with 0: the test is a strict `>`, so `threshold_pct=0` returns every channel flagged even
    once. Tool 7 spells 0 as "keep every channel" and therefore does not call this at all in that case."""
    return [c for c, b in zip(ch_names, np.asarray(badness_pct, dtype=float)) if b > float(threshold_pct)]


def recompute_reject(meta, methods_sel, stages_sel, event_types_sel=None,
                     pair_flags=None, ch_names=None, dropped_channels=(),
                     epoch_rule='any', epoch_rule_value=20.0):
    """Recompute the per-epoch reject decision from a USER-SELECTED subset of flagging methods, sleep
    stages and (for the 'event' method) event types — the tool-7 analogue of tool-7bis's build_pair_matrix,
    but read from the .fif metadata. An epoch is rejected when its stage is in `stages_sel` (in scope) AND at
    least one SELECTED method flags it; the 'event' method contributes the OR of the selected `evt_<type>`
    columns (falling back to `flag_event` when no per-type columns exist — pre-feature data). With every
    present method + stage + type selected this reproduces tool 6's stored `reject_flag`.

    HIGH-DENSITY channel layer (all optional; when `pair_flags` is None none of it runs and the result
    is exactly the pre-high-density one):
      `pair_flags`   {method: (n_ep, n_ch) bool} from load_channel_flags / recompute_channel_flags
      `ch_names`     channel names matching the pair_flags columns
      `dropped_channels`  channels the user dropped in the triage step: they no longer flag ANY epoch
      `epoch_rule`   'any' -> one flagged kept channel is enough (the classic rule)
                     'pct' -> more than `epoch_rule_value` % of the kept channels (tool 7bis's rule)
    Why this matters: with 3 channels 'any' is sound, with 32 it rejects 72.6 % of a real night, a third
    of it for a single electrode. Dropping the bad channel and/or switching to 'pct' is what makes the
    manual review tractable — see the measured figures in the HD constants block at the top.

    Returns (base_reject, reject_method, in_scope):
      base_reject   (n_ep,) bool  — in-scope AND flagged by a selected method
      reject_method (n_ep,) str   — the selected method that flagged it, 'multiple' for >=2, '' otherwise
                                    ('event' names the event contribution)
      in_scope      (n_ep,) bool  — stage in stages_sel (the epochs the report/navigator restrict to)"""
    stages = meta['stage'].astype(str).values
    n = len(stages)
    stages_sel = set(str(s) for s in stages_sel)
    in_scope = np.array([s in stages_sel for s in stages], dtype=bool)
    event_types_sel = list(event_types_sel) if event_types_sel else []
    # per-method boolean hits (only methods actually present as flag_<m> columns)
    hits = {}
    for m in methods_sel:
        if m == 'event':
            evt_cols = ['evt_' + t for t in event_types_sel if 'evt_' + t in meta.columns]
            if evt_cols:
                col = np.zeros(n, dtype=bool)
                for c in evt_cols:
                    col |= meta[c].astype(bool).values
            elif 'flag_event' in meta.columns:
                col = meta['flag_event'].astype(bool).values     # fallback: all event types
            else:
                col = np.zeros(n, dtype=bool)
            hits['event'] = col
        elif 'flag_' + m in meta.columns:
            hits[m] = meta['flag_' + m].astype(bool).values
    any_hit = np.zeros(n, dtype=bool)
    for col in hits.values():
        any_hit |= col
    # --- Channel-aware path (high density) -------------------------------------------------------
    # Recompute the per-epoch hits from the (epoch x channel) matrix restricted to the KEPT channels,
    # then apply the epoch rule. Skipped entirely when pair_flags is None (the default), so the sparse
    # PSG path stays byte-identical.
    if pair_flags is not None and ch_names is not None:
        dropped = {str(c) for c in (dropped_channels or ())}
        kept = [i for i, c in enumerate(ch_names) if str(c) not in dropped]
        n_kept = len(kept)
        pair_hits = {m: np.asarray(pair_flags[m], dtype=bool)[:, kept]
                     for m in methods_sel if m != 'event' and m in pair_flags}
        # Selected methods with NO per-channel column keep their epoch-level flag, broadcast across
        # every kept channel (tool 7bis's convention). That is always 'event', and also the two 1/f
        # methods when the flags were recomputed from the signal (time-domain only) instead of read
        # from the reports TSV — without this their contribution would silently vanish.
        epoch_only = {m: hits[m] for m in methods_sel
                      if m in hits and (m == 'event' or m not in pair_flags)}
        comb = np.zeros((n, n_kept), dtype=bool)          # OR of the selected per-channel methods
        for M in pair_hits.values():
            comb |= M
        if n_kept:
            for col in epoch_only.values():
                comb |= col[:, None]
        n_flag = comb.sum(axis=1)
        if epoch_rule == 'pct':
            any_hit = (n_flag / max(n_kept, 1)) > (float(epoch_rule_value) / 100.0)
        else:                                             # 'any' — the classic rule
            any_hit = n_flag > 0
        # attribution: the selected methods that flagged at least one KEPT channel in that epoch
        hits = {m: M.any(axis=1) for m, M in pair_hits.items()}
        hits.update(epoch_only)
    base_reject = in_scope & any_hit
    # per-epoch attribution (stable method order)
    reject_method = np.full(n, '', dtype=object)
    order = [m for m in METHOD_ORDER if m in hits]
    for ei in np.where(base_reject)[0]:
        got = [m for m in order if hits[m][ei]]
        reject_method[ei] = got[0] if len(got) == 1 else ('multiple' if len(got) >= 2 else '')
    return base_reject, reject_method.astype(str), in_scope


# ---------------------------------------------------------------------------
# Section 2 — per-stage report figures (shared by notebook + batch)
# ---------------------------------------------------------------------------
def _draw_psd_overlay(ax, freqs, psds_uV2, stages, reject_flag, reject_method, stage, custom_stages=()):
    """Draw one stage's per-epoch mean-across-channel PSD overlay into `ax`.
    The clean median + IQR band are drawn ON TOP (foreground) so a clean reference stands out; the
    individual clean traces (kept, per user preference) and the method-coloured rejected traces sit
    BEHIND them, attenuated. Returns nothing."""
    sel = (np.asarray(stages) == stage)
    if sel.sum() == 0:
        ax.text(0.5, 0.5, f'No {stage} epochs', ha='center', va='center')
        ax.set_axis_off()
        return
    psd_ep = np.asarray(psds_uV2)[sel].mean(axis=1)      # (n_stage_ep, n_freqs), µV²/Hz
    rej = np.asarray(reject_flag)[sel]
    meth = np.asarray(reject_method)[sel]
    # rejected epochs coloured by method — attenuated, BEHIND the clean median/IQR
    seen = set()
    for row, mth in zip(psd_ep[rej], meth[rej]):
        col = HEATMAP_COLORS[MULTIPLE_CODE] if mth == 'multiple' else METHOD_COLOR.get(mth, '#c0392b')
        lab = None
        if mth and mth not in seen:
            seen.add(mth)
            lab = mth
        ax.semilogy(freqs, row, color=col, lw=0.7, alpha=0.45, zorder=2, label=lab)
    # clean epochs: faint individual traces, then median + IQR ON TOP (foreground)
    clean = psd_ep[~rej]
    for row in clean:
        ax.semilogy(freqs, row, color='0.8', lw=0.4, alpha=0.30, zorder=1)
    if clean.shape[0] >= 3:
        med = np.median(clean, axis=0)
        q1, q3 = np.percentile(clean, [25, 75], axis=0)
        ax.fill_between(freqs, q1, q3, color='0.55', alpha=0.40, zorder=4, label='clean IQR')
        ax.semilogy(freqs, med, color='0.10', lw=2.0, zorder=5, label='clean median')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('PSD (µV²/Hz)')
    ax.set_title(f'{stage} — {int(sel.sum())} epochs, {int(rej.sum())} rejected', fontsize=9)
    ax.legend(fontsize=6, ncol=2, loc='upper right')


def plot_psd_overlay(freqs, psds_uV2, stages, reject_flag, reject_method, stage, custom_stages=(),
                     title_prefix=''):
    """Single-stage PSD overlay as its own Figure (thin wrapper over `_draw_psd_overlay`)."""
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    _draw_psd_overlay(ax, freqs, psds_uV2, stages, reject_flag, reject_method, stage, custom_stages)
    if title_prefix:
        ax.set_title(title_prefix + ax.get_title(), fontsize=10)
    fig.tight_layout()
    return fig


def plot_psd_overlays_grid(freqs, psds_uV2, stages, reject_flag, reject_method, stage_order,
                           custom_stages=()):
    """All per-stage PSD overlays in ONE figure on a 2-column grid (shorter notebook + report).
    Clean median/IQR in the foreground; clean + rejected traces attenuated behind. Returns a Figure."""
    stage_order = list(stage_order)
    n = len(stage_order)
    if n == 0:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, 'No stages in scope', ha='center', va='center'); ax.set_axis_off()
        return fig
    ncol = 2 if n > 1 else 1
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.4 * ncol, 3.4 * nrow), squeeze=False)
    for k, st in enumerate(stage_order):
        r, c = divmod(k, ncol)
        _draw_psd_overlay(axes[r][c], freqs, psds_uV2, stages, reject_flag, reject_method, st, custom_stages)
    for k in range(n, nrow * ncol):                       # blank any unused grid cell
        r, c = divmod(k, ncol)
        axes[r][c].set_axis_off()
    fig.tight_layout()
    return fig


def plot_metric_distributions(metrics, stages, reject_flag, thresholds, stage_order, title_prefix='',
                              reject_method=None, fit_note=''):
    """Box/strip of the four rejection metrics (p-p, gradient, 1/f MAE, 1/f R²), clean vs rejected,
    per stage, with the threshold line drawn. When `reject_method` is given, each REJECTED point is
    coloured by the method that flagged it (inside the red box); default None keeps the uniform red
    (batch twin unchanged). Returns a Figure."""
    specs = [('ptp', 'Peak-to-peak (µV)', 'amp'),
             ('gradient', 'Max |gradient| (µV/sample)', thresholds['gradient_uV_per_sample']),
             ('mae', '1/f MAE', thresholds['1f_mae_max']),
             ('r2', '1/f R²', thresholds['1f_r2_min'])]
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.2))
    stages = np.asarray(stages)
    rej = np.asarray(reject_flag)
    for ax, (key, ylab, thr) in zip(axes, specs):
        vals = np.asarray(metrics[key], dtype=float)
        positions, ticklabels = [], []
        for i, st in enumerate(stage_order):
            base_x = i * 3
            for j, (mask_lbl, colour, sub) in enumerate(
                    [('clean', '#2c7fb8', ~rej), ('rejected', '#c0392b', rej)]):
                m = (stages == st) & sub
                vv = vals[m]
                keep = ~np.isnan(vv)
                v = vv[keep]
                x = base_x + j
                positions.append(x)
                if v.size:
                    jit = np.full(v.size, x) + np.random.uniform(-0.15, 0.15, v.size)
                    if j == 1 and reject_method is not None:
                        # colour each rejected point by the method that flagged it
                        rmv = np.asarray(reject_method)[m][keep]
                        pcols = [HEATMAP_COLORS[MULTIPLE_CODE] if mm == 'multiple'
                                 else METHOD_COLOR.get(mm, '#c0392b') for mm in rmv]
                        ax.scatter(jit, v, s=8, c=pcols, alpha=0.6, zorder=3)
                    else:
                        ax.scatter(jit, v, s=4, color=colour, alpha=0.35, zorder=2)
                    ax.boxplot(v, positions=[x], widths=0.6, showfliers=False,
                               patch_artist=True,
                               boxprops=dict(facecolor='none', color=colour),
                               medianprops=dict(color=colour),
                               whiskerprops=dict(color=colour), capprops=dict(color=colour))
            ticklabels.append(st)
        # threshold line(s)
        if key == 'ptp':
            for i, st in enumerate(stage_order):
                t = thresholds['amplitude_ptp_uV'].get(st, 250.0)
                ax.plot([i * 3 - 0.5, i * 3 + 1.5], [t, t], color='k', ls='--', lw=1.0, zorder=3)
        else:
            ax.axhline(thr, color='k', ls='--', lw=1.0, zorder=3)
        ax.set_xticks([i * 3 + 0.5 for i in range(len(stage_order))])
        ax.set_xticklabels(ticklabels, fontsize=8)
        ax.set_ylabel(ylab, fontsize=9)
        ax.grid(True, axis='y', alpha=0.2)
    axes[0].set_title(f'{title_prefix}Metric distributions — clean (blue) vs rejected (red); '
                      f'dashed = threshold{fit_note}', fontsize=10, loc='left')
    fig.tight_layout()
    return fig


def plot_metric_scatter(metrics, reject_flag, reject_method, thresholds, title_prefix=''):
    """Scatter p-p vs max-gradient, coloured clean/grey or by reject method. Returns a Figure.

    Used by the DATABASE-POOLED section of `7_reject_manually_batch.py` only — it was dropped from the
    per-participant report (`build_participant_report_figs`), where the metric distributions already
    show both metrics with their thresholds and it never drove a decision."""
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    ptp = np.asarray(metrics['ptp'], dtype=float)
    grad = np.asarray(metrics['gradient'], dtype=float)
    rej = np.asarray(reject_flag)
    meth = np.asarray(reject_method)
    ax.scatter(ptp[~rej], grad[~rej], s=6, color='0.7', alpha=0.5, label='clean', zorder=1)
    for m in METHOD_ORDER + ['multiple']:
        sub = rej & (meth == m)
        if sub.sum():
            col = HEATMAP_COLORS[MULTIPLE_CODE] if m == 'multiple' else METHOD_COLOR.get(m, '#c0392b')
            ax.scatter(ptp[sub], grad[sub], s=10, color=col, alpha=0.7, label=m, zorder=2)
    ax.axhline(thresholds['gradient_uV_per_sample'], color='k', ls='--', lw=0.8, alpha=0.6)
    ax.set_xlabel('Peak-to-peak (µV)')
    ax.set_ylabel('Max |gradient| (µV/sample)')
    ax.set_title(f'{title_prefix}p-p vs gradient (method-coloured)', fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    return fig


def plot_channel_flag_heatmap(P, pair_flags, methods_sel, in_scope, badness_pct=None,
                              dropped_channels=(), custom_stages=(), flags_source=''):
    """Channels × epochs flagged-pair heatmap (coloured by flagging method) with a hypnogram strip on
    top and a per-channel badness bar on the right. The counterpart of tool 6's and tool 7bis's
    heatmaps, in tool 7's report — it is what makes "CPz is bad on 59 % of the night" visible BEFORE
    reviewing a single epoch, which is the decisive question on a dense montage.

    The matrix height scales with the channel count (same formula as 7bis) so a 3-channel PSG montage
    stays a short strip and a 64-channel one stays under ~12 in. `badness_pct` (n_ch,) is the bar; when
    None it is derived from `pair_flags` over the in-scope epochs. Dropped channels get a red bold
    label. `flags_source` is echoed in the title ('reports TSV' vs 'recomputed, time-domain only').
    Returns a Figure."""
    ch_names = list(P['ch_names'])
    n_ch = len(ch_names)
    stages = np.asarray(P['stages'])
    n_ep = len(stages)
    insc = np.asarray(in_scope, dtype=bool)
    used = [m for m in METHOD_ORDER if m in pair_flags and m in methods_sel]

    # Per-cell method code: 0 none, 1..6 the single method, 7 multiple (same code space as tool 6).
    n_hits = np.zeros((n_ch, n_ep), dtype=int)
    first = np.zeros((n_ch, n_ep), dtype=int)
    for m in used:
        Mm = np.asarray(pair_flags[m], dtype=bool).T                # (n_ch, n_ep)
        n_hits += Mm
        first = np.where((first == 0) & Mm, METHOD_CODE[m], first)
    code = np.where(n_hits >= 2, MULTIPLE_CODE, first)
    if badness_pct is None:
        badness_pct, _ = channel_badness(pair_flags, methods_sel, insc, n_ch)
    badness_pct = np.asarray(badness_pct, dtype=float)

    h_mat = float(np.clip(0.22 * n_ch, 1.0, 12.0))
    fig_w = float(np.clip(n_ep / 6.0, 9.0, 18.0)) + 2.6
    h_hyp = 0.9
    fig_h = h_hyp + h_mat + 1.2
    fig, axes = plt.subplots(2, 2, figsize=(fig_w, fig_h),
                             gridspec_kw={'height_ratios': [h_hyp, h_mat], 'width_ratios': [4.2, 1.0],
                                          'hspace': 0.06, 'wspace': 0.03})
    # --- hypnogram strip (shares the epoch axis with the matrix) ---
    ax_hyp = axes[0][0]
    stage_y, stage_colors, ytick_pos, ytick_labels = custom_stage_style(custom_stages)
    floor = min(ytick_pos) - 1
    hypno_y = np.array([stage_y.get(s, floor) for s in stages])
    ax_hyp.step(np.arange(n_ep + 1), np.append(hypno_y, hypno_y[-1]), where='post',
                color='#555555', lw=1.0)
    for ei in range(n_ep):                                          # REM in red (YASA convention)
        if stages[ei] == 'R':
            ax_hyp.plot([ei, ei + 1], [stage_y['R'], stage_y['R']], color='#c0392b', lw=2.0,
                        solid_capstyle='butt')
    ax_hyp.set_yticks(ytick_pos); ax_hyp.set_yticklabels(ytick_labels, fontsize=6)
    ax_hyp.set_xlim(0, n_ep); ax_hyp.set_xticks([])
    for sp in ['top', 'right', 'bottom']:
        ax_hyp.spines[sp].set_visible(False)
    # Title on the HYPNOGRAM axis (the TOP one) — on the matrix axis it was drawn inside the 6%-of-height
    # gap between the two axes, i.e. straight over the hypnogram strip.
    src = f' — flags: {flags_source}' if flags_source else ''
    ax_hyp.set_title(f'Flagged (epoch × channel) pairs — {int((code > 0).sum())} of {n_ep * n_ch} '
                     f'({100.0 * (code > 0).mean():.1f}%){src}', fontsize=9, pad=6)
    axes[0][1].set_axis_off()

    # --- channels x epochs matrix ---
    ax_m = axes[1][0]
    cmap = mcolors.ListedColormap(HEATMAP_COLORS)
    norm = mcolors.BoundaryNorm(np.arange(-0.5, len(HEATMAP_COLORS) + 0.5, 1), cmap.N)
    ax_m.imshow(code, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest',
                extent=[0, n_ep, n_ch - 0.5, -0.5])
    # label size follows the actual row height (inches -> points), so labels never collide
    lbl_fs = float(np.clip((h_mat / max(n_ch, 1)) * 72.0 * 0.7, 4.5, 8.0))
    ax_m.set_yticks(range(n_ch))
    ax_m.set_yticklabels(ch_names, fontsize=lbl_fs)
    dropped = {str(c) for c in (dropped_channels or ())}
    for tick, cn in zip(ax_m.get_yticklabels(), ch_names):
        if cn in dropped:                                           # dropped channels: red bold label
            tick.set_color('#c0392b'); tick.set_fontweight('bold')
    ax_m.set_xlabel('Epoch index', fontsize=9)

    # --- per-channel badness bar (right), sharing the channel axis ---
    ax_b = axes[1][1]
    colours = ['#c0392b' if cn in dropped else '#4a3aa7' for cn in ch_names]
    ax_b.barh(np.arange(n_ch), badness_pct, color=colours, height=0.8)
    ax_b.set_ylim(n_ch - 0.5, -0.5)
    ax_b.set_yticks([])
    ax_b.set_xlabel('% in-scope epochs', fontsize=8)
    ax_b.tick_params(axis='x', labelsize=7)
    ax_b.grid(True, axis='x', alpha=0.25)
    ax_b.set_title(f'Badness (n={int(insc.sum())})', fontsize=8)
    for sp in ['top', 'right']:
        ax_b.spines[sp].set_visible(False)

    handles = [Rectangle((0, 0), 1, 1, color=HEATMAP_COLORS[METHOD_CODE[m]], label=METHOD_LABEL[m])
               for m in used]
    if (n_hits >= 2).any():
        handles.append(Rectangle((0, 0), 1, 1, color=HEATMAP_COLORS[MULTIPLE_CODE], label='multiple'))
    if handles:
        ax_m.legend(handles=handles, loc='upper left', bbox_to_anchor=(0, -0.12),
                    ncol=min(len(handles), 7), fontsize=7, frameon=False)
    # Keep at least 0.32 in of top margin for the title (only binding on a short, few-channel figure —
    # a dense one is tall enough that 0.93 already leaves more).
    fig.subplots_adjust(left=0.09, right=0.98, top=min(0.93, 1.0 - 0.32 / fig_h), bottom=0.13)
    return fig


def build_rejection_table_html(meta, methods_present, stage_order, title='', reject_flag=None):
    """Per-stage × per-method rejection table (% and raw count), from the .fif metadata flags.
    `reject_flag` (optional (n_ep,) bool) overrides the 'Any' column / 'All' total with a recomputed
    decision (tool 7 passes its selected-methods `base_reject`); default None uses the tool-6 `reject_flag`.
    The 'All' total is taken over the epochs whose stage is in `stage_order` (the in-scope stages), so a
    stage sub-selection stays self-consistent."""
    rows = []
    header = ['Stage', 'N epochs'] + [METHOD_LABEL[m] for m in methods_present] + ['Any']
    stages = meta['stage'].astype(str).values
    rf = meta['reject_flag'].astype(bool).values if reject_flag is None else np.asarray(reject_flag, dtype=bool)
    for st in stage_order:
        sel = (stages == st)
        n = int(sel.sum())
        if n == 0:
            continue
        cells = [st, str(n)]
        for m in methods_present:
            c = int(meta.loc[sel, 'flag_' + m].astype(bool).sum())
            cells.append(f'{100.0 * c / n:.1f}% ({c})')
        a = int(rf[sel].sum())
        cells.append(f'<b>{100.0 * a / n:.1f}% ({a})</b>')
        rows.append(cells)
    # total row (over the in-scope stages only)
    scope = np.isin(stages, list(stage_order))
    n = int(scope.sum())
    tot = ['<b>All</b>', f'<b>{n}</b>']
    for m in methods_present:
        c = int(meta.loc[scope, 'flag_' + m].astype(bool).sum())
        tot.append(f'<b>{100.0 * c / n:.1f}% ({c})</b>' if n else '<b>—</b>')
    a = int(rf[scope].sum())
    tot.append(f'<b>{100.0 * a / n:.1f}% ({a})</b>' if n else '<b>—</b>')
    rows.append(tot)
    th = ''.join(f'<th style="padding:3px 8px;border:1px solid #ccc;">{h}</th>' for h in header)
    body = ''
    for r in rows:
        body += '<tr>' + ''.join(
            f'<td style="padding:3px 8px;border:1px solid #ccc;text-align:center;">{c}</td>' for c in r) + '</tr>'
    return (f'<h4>{title}</h4>' if title else '') + \
           f'<table style="border-collapse:collapse;font-size:.85em;"><tr>{th}</tr>{body}</table>'


def build_participant_report_figs(P, metrics, freqs, psds_uV2, thresholds, custom_stages,
                                  reject_flag=None, reject_method=None, stage_order=None, methods=None,
                                  pair_flags=None, in_scope=None, badness_pct=None,
                                  dropped_channels=(), flags_source=''):
    """Assemble the Section-2 figures for one participant. Returns (list_of_(title, fig), table_html).
    `reject_flag` / `reject_method` / `stage_order` / `methods` let tool 7 pass its RECOMPUTED
    selected-methods decision + in-scope stage order; all default to the tool-6 metadata (batch twin
    unchanged). When `pair_flags` is given the channels × epochs flagging heatmap is prepended (it is
    the per-channel overview; see plot_channel_flag_heatmap) — omitted when None, so the batch twin's
    output is unchanged."""
    meta = P['meta']
    stages = P['stages']
    if reject_flag is None:
        reject_flag = P['reject_flag']
    reject_flag = np.asarray(reject_flag)
    if reject_method is None:
        reject_method = meta['reject_method'].astype(str).values
    else:
        reject_method = np.asarray(reject_method).astype(str)
    if stage_order is None:
        stage_order = ordered_present_stages(stages, custom_stages)
    if methods is None:
        methods = P['methods_present']
    # State it when the 1/f fit ran on a subsample — a reader must never mistake it for the full night.
    n_fit, n_cand = metrics.get('n_fitted'), metrics.get('n_candidates')
    fit_note = ''
    if n_fit is not None and n_cand and n_fit < n_cand:
        fit_note = f'  |  1/f fitted on {n_fit} of {n_cand} in-scope epochs (random subsample)'
    figs = []
    if pair_flags is not None:
        figs.append(('Flagged pairs (channels × epochs) + per-channel badness',
                     plot_channel_flag_heatmap(P, pair_flags, methods,
                                               np.ones(len(stages), dtype=bool) if in_scope is None
                                               else in_scope,
                                               badness_pct=badness_pct,
                                               dropped_channels=dropped_channels,
                                               custom_stages=custom_stages,
                                               flags_source=flags_source)))
    figs.append(('PSD overlays (clean median/IQR in front)',
                 plot_psd_overlays_grid(freqs, psds_uV2, stages, reject_flag, reject_method,
                                        stage_order, custom_stages=custom_stages)))
    figs.append(('Metric distributions',
                 plot_metric_distributions(metrics, stages, reject_flag, thresholds, stage_order,
                                           reject_method=reject_method, fit_note=fit_note)))
    # No p-p vs gradient scatter here: it only adds the JOINT distribution of two metrics the figure
    # above already shows separately (with thresholds and method colours), which is an exploration
    # view, not a decision one. It is kept POOLED OVER THE DATABASE in 7_reject_manually_batch.py,
    # where comparing the two metrics across participants does carry information.
    table_html = build_rejection_table_html(meta, methods, stage_order,
                                            title='Rejection by stage × method', reject_flag=reject_flag)
    return figs, table_html


# ---------------------------------------------------------------------------
# Section 3 — per-epoch figures (used by the interactive navigator)
# ---------------------------------------------------------------------------
def _epoch_channel_method(cur_ch_uV, stage, thresholds):
    """Recompute which method (if any) flags one channel of one epoch (tool-6 formulas + thresholds).
    Priority amplitude > gradient > flat (spectral flags are handled in the detail panel)."""
    ptp = float(np.ptp(cur_ch_uV))
    grad = float(np.max(np.abs(np.diff(cur_ch_uV))))
    m = None
    if ptp < thresholds['flat_ptp_uV']:
        m = 'flat'
    if grad > thresholds['gradient_uV_per_sample']:
        m = 'gradient'
    if ptp > thresholds['amplitude_ptp_uV'].get(stage, 250.0):
        m = 'amplitude'
    return m, ptp, grad


def select_montage_channels(P, ei, thresholds, mode='all', n_worst=16, neighbours=2, custom=None):
    """Channel INDICES to draw in the per-epoch montage of epoch `ei`, always returned in FILE ORDER
    (never re-sorted by badness, so the spatial reading of the montage is preserved).

    On a 32-64 channel montage the median rejected epoch is flagged on ONE channel, so drawing all of
    them to find it is backwards - hence the 'flagged' modes. Modes:
      'all'                -> every channel (the classic behaviour)
      'flagged'            -> only the channels a time-domain method flags in this epoch
      'flagged+neighbours' -> those, plus their +/-`neighbours` file-order neighbours
      'worst'              -> the `n_worst` channels with the largest peak-to-peak in this epoch
      'custom'             -> the explicit `custom` list (channel names or indices)
    Never returns an empty list: an epoch with no flagged channel falls back to the `n_worst` worst,
    so the montage always shows something."""
    ch_names = P['ch_names']
    n_ch = len(ch_names)
    cur = P['data_uV'][ei]
    if mode == 'all':
        return list(range(n_ch))
    if mode == 'custom':
        wanted = list(custom or [])
        idx = [ch_names.index(c) if isinstance(c, str) and c in ch_names else int(c)
               for c in wanted if (isinstance(c, str) and c in ch_names) or not isinstance(c, str)]
        idx = sorted({i for i in idx if 0 <= i < n_ch})
        return idx or list(range(n_ch))
    ptp = np.ptp(cur, axis=-1)
    if mode == 'worst':
        return sorted(np.argsort(ptp)[::-1][:max(1, int(n_worst))].tolist())
    stg = P['stages'][ei]
    flagged = [c for c in range(n_ch) if _epoch_channel_method(cur[c], stg, thresholds)[0] is not None]
    if not flagged:                       # nothing flagged (e.g. a 1/f-only epoch): show the worst ones
        return sorted(np.argsort(ptp)[::-1][:max(1, int(n_worst))].tolist())
    if mode == 'flagged+neighbours':
        keep = set()
        for c in flagged:
            keep.update(range(max(0, c - int(neighbours)), min(n_ch, c + int(neighbours) + 1)))
        return sorted(keep)
    return sorted(flagged)


def plot_epoch_montage(P, ei, thresholds, context=1, ctx=None, onsets=None, title_method=None,
                       scales=None, show_channels=None):
    """Stacked montage of epoch `ei` (± `context` epochs), MNE-raw-plot style. EEG channels on top
    (current-epoch trace coloured by its recomputed flagging method, steepest-gradient jump boxed), then
    the optional EOG/EMG/ECG context traces below — **each channel TYPE on its own FIXED amplitude scale**
    (`DISPLAY_SCALE_UV`, clinical conventions), with a left-margin **scale bar (µV)** per type. The scales
    do not adapt to the window, so signals are directly comparable from one epoch to the next. Scored-event
    onsets are drawn as labelled vertical lines; a right-side legend gives the colour code. Returns a Figure.

    `ctx`   : dict from load_context_epochs (EOG/EMG/ECG epoched 1:1 with the EEG); None -> EEG only.
    `onsets`: DataFrame [type, onset_s, duration_s] from load_event_onsets; None -> no event markers.
    `title_method`: reject-method label to show (tool 7 passes its recomputed value); None -> tool-6 meta.
    `scales`: optional {type: µV per row} overriding DISPLAY_SCALE_UV (e.g. {'EEG': 200}).
    `show_channels`: optional list of channel INDICES to draw (see select_montage_channels) - the
        high-density lever; None (default) draws every channel, i.e. the classic behaviour."""
    data, sf, ch_all = P['data_uV'], P['sfreq'], P['ch_names']
    stages, meta = P['stages'], P['meta']
    n_ep, n_ch_all, n_t = data.shape
    # Channel subset (high-density): a plain gather, so drawing every channel is numerically identical
    # to the former code. Kept in file order; an empty/invalid selection falls back to all channels.
    if show_channels is None:
        ch_idx = list(range(n_ch_all))
    else:
        ch_idx = sorted({int(i) for i in show_channels if 0 <= int(i) < n_ch_all})
        if not ch_idx:
            ch_idx = list(range(n_ch_all))
    ch = [ch_all[i] for i in ch_idx]
    n_ch = len(ch_idx)
    lo, hi = max(0, ei - context), min(n_ep - 1, ei + context)
    idxs = list(range(lo, hi + 1))
    seg = np.concatenate([data[k][ch_idx] for k in idxs], axis=1)   # (n_ch, T)
    t = np.arange(seg.shape[1]) / sf
    stg = stages[ei]; cur = data[ei][ch_idx]; cur_pos = idxs.index(ei)

    # context rows aligned to the same window, grouped by physiological type (EOG/EMG/ECG)
    ctx_rows = []                                               # (label, type, series, sfreq)
    if ctx is not None and len(ctx.get('labels', [])):
        cdata, csf, clabels = ctx['data_uV'], ctx['sfreq'], ctx['labels']
        cidxs = [k for k in idxs if k < cdata.shape[0]]
        if cidxs:
            cseg = np.concatenate([cdata[k] for k in cidxs], axis=1)
            for j, cl in enumerate(clabels):
                ctyp = 'EOG' if str(cl).upper().startswith('EOG') else str(cl).split('-')[0].upper()
                ctx_rows.append((cl, ctyp, cseg[j], csf))
    n_ctx = len(ctx_rows)
    n_total = n_ch + n_ctx
    row_h, fill = 1.0, 0.8

    # --- High-density geometry ---------------------------------------------------------------
    # Only the INCHES PER ROW shrink; `row_h` stays 1.0 in data units, so `g = fill*row_h/type_scale`
    # is untouched and the FIXED clinical µV scales keep their meaning (a 75 µV slow wave still fills
    # exactly half a row, overflow into the neighbouring row still happens as in a clinical viewer).
    # The clip is a no-op up to 26 rows, so the 3-6 channel PSG montages render byte-identically.
    row_in = float(np.clip(MONTAGE_MAX_H_IN / max(n_total, 1), MONTAGE_ROW_MIN_IN, MONTAGE_ROW_IN))
    lbl_fs = float(np.clip(8.0 * row_in / MONTAGE_ROW_IN, 5.0, 8.0))
    trace_lw = 0.8 if row_in >= 0.35 else 0.55
    # Width follows the montage DENSITY, not the displayed subset, so the figure does not resize when
    # the user switches between 'all' and 'flagged only'.
    fig_w = 12.0 if n_ch_all <= HD_CHANNEL_THRESHOLD else 14.0

    # FIXED per-type amplitude scale (µV peak-to-peak per row) — deliberately NOT derived from the
    # displayed window, so the same waveform keeps the same height across epochs and visual amplitude
    # criteria (75 µV slow waves, EMG tone, …) remain comparable. `scales` overrides the defaults.
    type_scale = dict(DISPLAY_SCALE_UV)
    if scales:
        type_scale.update({k: float(v) for k, v in scales.items() if v})
    for typ in {r[1] for r in ctx_rows}:
        type_scale.setdefault(typ, DISPLAY_SCALE_UV['EEG'])   # unexpected context type: fall back to EEG

    fig, ax = plt.subplots(figsize=(fig_w, row_in * n_total + 1.9))
    ax.axvspan(cur_pos * n_t / sf, (cur_pos + 1) * n_t / sf, color='#fff3cd', alpha=0.7, zorder=0)
    off_of = lambda r: (n_total - 1 - r) * row_h               # row 0 = top
    tt = np.arange(cur_pos * n_t, (cur_pos + 1) * n_t) / sf
    g_eeg = fill * row_h / type_scale['EEG']
    methods_here = set()
    for i, cn in enumerate(ch):
        off = off_of(i); base = float(np.median(seg[i]))
        ax.plot(t, (seg[i] - base) * g_eeg + off, color='0.6', lw=0.5, zorder=2)
        m, ptp, grad = _epoch_channel_method(cur[i], stg, thresholds)
        if m is not None:
            methods_here.add(m)
        col = '0.2' if m is None else METHOD_COLOR[m]
        ax.plot(tt, (cur[i] - base) * g_eeg + off, color=col, lw=trace_lw, zorder=3)
        if grad > thresholds['gradient_uV_per_sample']:
            j = int(np.argmax(np.abs(np.diff(cur[i]))))
            # Frame the steepest sample-to-sample jump with a hollow box (keeps the trace visible),
            # redraw the jump edge in red.
            y_lo, y_hi = sorted([(cur[i][j] - base) * g_eeg + off, (cur[i][j + 1] - base) * g_eeg + off])
            pad_x, pad_y = 5.0 / sf, 0.12 * row_h
            ax.add_patch(Rectangle((tt[j] - pad_x, y_lo - pad_y), (tt[j + 1] - tt[j]) + 2 * pad_x,
                                   (y_hi - y_lo) + 2 * pad_y, fill=False,
                                   edgecolor=METHOD_COLOR['gradient'], lw=1.2, zorder=4))
            ax.plot(tt[j:j + 2], (cur[i][j:j + 2] - base) * g_eeg + off, color='red', lw=1.4, zorder=5)
        ax.text(-0.008, off, cn, ha='right', va='center', fontsize=lbl_fs,
                color=('k' if m is None else METHOD_COLOR[m]), transform=ax.get_yaxis_transform())
    ctx_types_present = []
    for j, (cl, ctyp, series, csf) in enumerate(ctx_rows):
        off = off_of(n_ch + j); base = float(np.median(series))
        g = fill * row_h / type_scale[ctyp]
        tc = np.arange(len(series)) / csf
        col = CTX_COLOR.get(cl, '#8e44ad')
        ax.plot(tc, (series - base) * g + off, color=col, lw=0.6, zorder=2)
        ax.text(-0.008, off, cl, ha='right', va='center', fontsize=lbl_fs, color=col,
                transform=ax.get_yaxis_transform())
        if ctyp not in ctx_types_present:
            ctx_types_present.append(ctyp)
    if n_ctx:
        ax.axhline(off_of(n_ch) + 0.6 * row_h, color='0.4', lw=0.8, ls=':', zorder=1)   # EEG|context divider

    # per-type scale bars in the left margin (axes-fraction x via get_yaxis_transform)
    def _scalebar(rows_off, typ):
        if not rows_off:
            return
        yc = float(np.mean(rows_off)); half = fill * row_h / 2.0; xb = -0.085
        trans = ax.get_yaxis_transform()
        ax.plot([xb, xb], [yc - half, yc + half], color='k', lw=1.3, transform=trans, clip_on=False)
        for yy in (yc - half, yc + half):
            ax.plot([xb - 0.006, xb + 0.006], [yy, yy], color='k', lw=1.0, transform=trans, clip_on=False)
        ax.text(xb - 0.012, yc, f'{typ}\n{type_scale[typ]:.0f} µV', ha='right', va='center', fontsize=7,
                transform=trans, clip_on=False)
    _scalebar([off_of(i) for i in range(n_ch)], 'EEG')
    for typ in ctx_types_present:
        _scalebar([off_of(n_ch + j) for j, r in enumerate(ctx_rows) if r[1] == typ], typ)

    # scored-event onset markers (from the tool-6 _event_onsets.tsv sidecar): a full-height vertical line
    # with the event name written HORIZONTALLY just under the top border of the plot (axes-fraction y, so
    # it hugs the top edge regardless of the amplitude scaling).
    epoch_len = n_t / sf; win_start = lo * epoch_len
    drew_event = False
    if onsets is not None and len(onsets):
        for _, evr in onsets.iterrows():
            try:
                x = float(evr['onset_s']) - win_start
            except Exception:
                continue
            if 0 <= x <= t[-1]:
                ax.axvline(x, color=METHOD_COLOR['event'], lw=1.1, alpha=0.85, zorder=6)
                ax.text(x + 0.2, 0.995, str(evr['type']), transform=ax.get_xaxis_transform(),
                        va='top', ha='left', fontsize=7, color=METHOD_COLOR['event'], clip_on=True, zorder=7)
                drew_event = True

    for k in range(len(idxs) + 1):
        ax.axvline(k * n_t / sf, color='0.85', lw=0.6, zorder=1)
    ax.set_yticks([]); ax.set_xlabel('Time (s)'); ax.set_xlim(t[0], t[-1])
    ax.set_ylim(off_of(n_total - 1) - row_h, off_of(0) + 1.3 * row_h)
    rm = title_method if title_method is not None else meta['reject_method'].astype(str).values[ei]
    scale_txt = ', '.join(f'{k} {type_scale[k]:.0f} µV'
                          for k in ['EEG'] + ctx_types_present if k in type_scale)
    # State the subset explicitly: on a dense montage channels are hidden, and that must never be silent.
    sub_txt = f'showing {n_ch}/{n_ch_all} channels; ' if n_ch < n_ch_all else ''
    ax.set_title(f'Epoch {ei} — stage {stg} — reject_method={rm}   '
                 f'({sub_txt}context ±{context}; fixed scales: {scale_txt}; '
                 f'coloured = flagging method)', fontsize=10)

    # right-side legend (colour code): flagging methods used in this epoch + context types + event marker
    handles = [Line2D([0], [0], color='0.2', lw=2, label='clean (kept)')]
    handles += [Line2D([0], [0], color=METHOD_COLOR[m], lw=2, label=m)
                for m in METHOD_ORDER if m in methods_here]
    for typ in ctx_types_present:
        rep = next((CTX_COLOR.get(r[0], '#8e44ad') for r in ctx_rows if r[1] == typ), '#8e44ad')
        handles.append(Line2D([0], [0], color=rep, lw=2, label=f'{typ} (context)'))
    if drew_event:
        handles.append(Line2D([0], [0], color=METHOD_COLOR['event'], lw=1.5, label='event onset'))
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.005, 1.0), fontsize=7, framealpha=0.9)
    fig.subplots_adjust(left=0.17, right=0.84, top=0.93, bottom=0.08)
    return fig


def has_positions(P):
    """True when the epochs carry usable electrode positions (Curry `.cdt` keeps a real DigMontage
    through tool 6; Compumedics EDF has none). Gates the topomap panel — non-fatal when absent."""
    try:
        loc = np.array([c['loc'][:3] for c in P['epochs'].info['chs']], dtype=float)
    except Exception:
        return False
    return bool(loc.shape[0] and np.isfinite(loc).all() and (np.abs(loc).sum(axis=1) > 0).all())


def _draw_channel_topomap(ax, P, values, title, cmap='viridis', mask=None):
    """Small per-channel topomap of `values` (one per channel) for the current epoch. `mask` (bool per
    channel) rings the flagged electrodes in red so a bad sensor is identifiable without reading the
    colour scale. Non-fatal: any MNE/layout failure prints a note in the panel instead of raising."""
    try:
        info = P['epochs'].info
        vals = np.asarray(values, dtype=float)
        good = np.isfinite(vals)
        if not good.any():
            raise ValueError('no finite value')
        vals = np.where(good, vals, float(np.median(vals[good])))   # NaN 1/f fits -> neutral colour
        im, _ = mne.viz.plot_topomap(
            vals, info, axes=ax, show=False, cmap=cmap, sensors=True, contours=0, outlines='head',
            mask=None if mask is None else np.asarray(mask, dtype=bool),
            mask_params=dict(marker='o', markerfacecolor='none', markeredgecolor='#c0392b',
                             linewidth=0, markersize=9, markeredgewidth=2.0))
        cb = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.ax.tick_params(labelsize=6)
        ax.set_title(title, fontsize=8)
    except Exception as e:
        ax.set_axis_off()
        ax.text(0.5, 0.5, f'topomap n/a\n({e})', ha='center', va='center', fontsize=7, color='0.4')


def plot_epoch_detail(P, ei, freqs, psds_uV2, thresholds, spectro_ch_idx=0, fmin=2.0,
                      clean_mask=None, table_rows=None, topomap=None):
    """2×2 detail panel for epoch `ei`: per-channel PSD (top-left) with the per-channel metric table right
    beside it (top-right); mean band power + 50 Hz ratio and the epoch spectrogram on the bottom row
    (the spectrogram's channel selector sits under the figure in the navigator).

    PSD panel (same units/scale as the Section-2 overlay — semilogy, µV²/Hz, linear frequency axis, so
    non-experts read both plots the same way):
      * background = median + IQR of the CLEAN epochs OF THE SAME SLEEP STAGE (mean across channels),
        so it is immediately visible whether this epoch's spectrum sits inside the normal range;
      * each channel whose 1/f fit is flagged (MAE > threshold or R² < threshold) is drawn in its METHOD
        colour and gets its aperiodic fit dashed; unflagged channels stay grey (no 'worst channel'
        emphasis anymore — what matters is which channels actually breach a threshold).
    Table: value (threshold) per metric, the value shown in RED when it breaches its threshold.

    `clean_mask` (optional (n_ep,) bool): epochs to use as the clean reference; default = not rejected
    according to the tool-6 metadata.
    `table_rows` (optional int): HIGH-DENSITY lever — keep every channel that breaches a threshold plus
        this many worst-by-p-p ones in the metric table, instead of one row per channel (illegible past
        ~15 rows). None (default) = every channel, i.e. the classic behaviour.
    `topomap` (optional bool): draw the spatial panels (peak-to-peak + 1/f MAE). None (default) = auto:
        on when the epochs carry electrode positions AND the montage is dense. Curry keeps a real
        DigMontage through tool 6; position-less EDF simply never gets the panel."""
    from scipy.signal import welch as _welch
    data, sf, ch = P['data_uV'], P['sfreq'], P['ch_names']
    stg = P['stages'][ei]
    n_ch = len(ch)
    if topomap is None:
        topomap = (n_ch > HD_CHANNEL_THRESHOLD) and has_positions(P)
    fits = [fit_1f(freqs, psds_uV2[ei, c], fmin=fmin) for c in range(n_ch)]
    mae_thr = thresholds['1f_mae_max']
    r2_thr = thresholds['1f_r2_min']
    amp_thr = thresholds['amplitude_ptp_uV'].get(stg, 250.0)
    grad_thr = thresholds['gradient_uV_per_sample']
    flat_thr = thresholds['flat_ptp_uV']

    # which 1/f method (if any) flags each channel -> trace colour
    def _fit_flag(c):
        mae_c, r2_c = fits[c][0], fits[c][1]
        if not np.isnan(r2_c) and r2_c < r2_thr:
            return '1f_r2'
        if not np.isnan(mae_c) and mae_c > mae_thr:
            return '1f_error'
        return None

    # Layout: PSD (top-left) with the metric table right beside it (top-right); mean band power and the
    # epoch spectrogram on the bottom row. The spectrogram's channel selector sits under the figure in the UI.
    # A third column carries the two topomaps when positions are available (high-density montages).
    if topomap:
        fig = plt.figure(figsize=(16.5, 8))
        gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 1.0, 0.75], height_ratios=[1.15, 1.0],
                              wspace=0.25, hspace=0.35)
        ax_topo_ptp = fig.add_subplot(gs[0, 2])
        ax_topo_mae = fig.add_subplot(gs[1, 2])
    else:
        fig = plt.figure(figsize=(13, 8))
        gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0], height_ratios=[1.15, 1.0],
                              wspace=0.2, hspace=0.35)
        ax_topo_ptp = ax_topo_mae = None
    ax_psd = fig.add_subplot(gs[0, 0])
    ax_table = fig.add_subplot(gs[0, 1])
    ax_band = fig.add_subplot(gs[1, 0])
    ax_spec = fig.add_subplot(gs[1, 1])

    # (a) PSD — clean stage reference in the background, flagged channels highlighted
    ax = ax_psd
    clean = np.asarray(clean_mask, dtype=bool) if clean_mask is not None else ~P['reject_flag']
    ref = np.asarray(psds_uV2)[(np.asarray(P['stages']) == stg) & clean]
    if ref.shape[0] >= 3:
        ref_ep = ref.mean(axis=1)                       # (n_clean_ep, n_freqs), mean across channels
        med = np.median(ref_ep, axis=0)
        q1, q3 = np.percentile(ref_ep, [25, 75], axis=0)
        ax.fill_between(freqs, q1, q3, color='0.55', alpha=0.35, zorder=1,
                        label=f'clean {stg} IQR (n={ref_ep.shape[0]})')
        ax.semilogy(freqs, med, color='0.10', lw=1.8, zorder=2, label=f'clean {stg} median')
    n_flagged = 0
    labelled_grey = False
    hd = n_ch > HD_CHANNEL_THRESHOLD
    unflagged_psd = []            # high-density: collected and drawn once as a p5–p95 band
    max_flagged_labels = 12       # keep the legend readable when many channels breach at once
    for c in range(n_ch):
        fl = _fit_flag(c)
        if fl is None:
            if hd:
                # One grey line per channel is both unreadable and a legend flood at 32–64 channels:
                # the unflagged population is drawn as a band after the loop instead.
                unflagged_psd.append(psds_uV2[ei, c])
                continue
            # unflagged channels of THIS epoch: grey, labelled once so the legend explains them
            ax.semilogy(freqs, psds_uV2[ei, c], color='0.75', lw=0.8, zorder=4,
                        label=None if labelled_grey else 'this epoch — channels not flagged')
            labelled_grey = True
        else:
            n_flagged += 1
            col = METHOD_COLOR[fl]
            off, exp = fits[c][2], fits[c][3]
            lab = None
            if n_flagged <= max_flagged_labels:
                lab = f'{ch[c]} ({METHOD_LABEL[fl]}'
                lab += f', exp={exp:.2f})' if not np.isnan(exp) else ')'
            ax.semilogy(freqs, psds_uV2[ei, c], color=col, lw=1.4, zorder=5, label=lab)
            if not np.isnan(off):                        # aperiodic fit of the flagged channel only
                fm = freqs >= fmin
                ax.semilogy(freqs[fm], 10 ** (off - exp * np.log10(freqs[fm])),
                            color=col, lw=1.0, ls='--', alpha=0.8, zorder=6)
    if unflagged_psd:
        arr = np.asarray(unflagged_psd)
        lo_env, hi_env = np.percentile(arr, [5, 95], axis=0)
        ax.fill_between(freqs, lo_env, hi_env, color='0.78', alpha=0.55, zorder=3,
                        label=f'this epoch — channels not flagged (n={len(unflagged_psd)}, p5–p95)')
        ax.semilogy(freqs, np.median(arr, axis=0), color='0.55', lw=1.0, zorder=4)
    if n_flagged > max_flagged_labels:
        ax.plot([], [], ' ', label=f'… +{n_flagged - max_flagged_labels} more flagged')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('PSD (µV²/Hz)')
    ax.set_title(f'PSD vs clean {stg} reference — {n_flagged} channel(s) flagged by the 1/f fit '
                 f'(dashed = aperiodic fit)', fontsize=9)
    ax.legend(fontsize=6, ncol=3 if hd else 2, loc='upper right')

    # (b) per-channel metric table, right beside the PSD (red = breaches its threshold)
    ax = ax_table
    ax.set_axis_off()
    cur = data[ei]
    # Vectorised per-channel values (identical to the former per-channel float() calls) — also feed
    # the topomaps below.
    ptp_all = np.ptp(cur, axis=-1)
    grad_all = np.max(np.abs(np.diff(cur, axis=-1)), axis=-1)
    mae_all = np.array([fits[c][0] for c in range(n_ch)], dtype=float)
    r2_all = np.array([fits[c][1] for c in range(n_ch)], dtype=float)
    breached = {c for c in range(n_ch)
                if (ptp_all[c] > amp_thr or ptp_all[c] < flat_thr or grad_all[c] > grad_thr
                    or (not np.isnan(mae_all[c]) and mae_all[c] > mae_thr)
                    or (not np.isnan(r2_all[c]) and r2_all[c] < r2_thr))}
    # Which channels get a row: all of them (classic), or — on a dense montage — every breaching
    # channel plus the `table_rows` worst by peak-to-peak, kept in file order.
    if table_rows is None:
        show_rows = list(range(n_ch))
    else:
        worst = np.argsort(ptp_all)[::-1][:max(0, int(table_rows))]
        show_rows = sorted(breached | {int(w) for w in worst})
    header = ['ch', 'p-p (thr)', 'grad (thr)', 'MAE (thr)', 'R² (thr)']
    table = [header]
    bad_cells = []                                       # (row, col) to paint red
    for row, c in enumerate(show_rows, start=1):
        ptp, grad = float(ptp_all[c]), float(grad_all[c])
        mae_c, r2_c = fits[c][0], fits[c][1]
        table.append([ch[c],
                      f'{ptp:.0f} ({amp_thr:.0f})',
                      f'{grad:.0f} ({grad_thr:.0f})',
                      f'{mae_c:.3f} ({mae_thr:.2f})',
                      f'{r2_c:.3f} ({r2_thr:.2f})'])
        if ptp > amp_thr or ptp < flat_thr:              # amplitude OR flat both breach the p-p cell
            bad_cells.append((row, 1))
        if grad > grad_thr:
            bad_cells.append((row, 2))
        if not np.isnan(mae_c) and mae_c > mae_thr:
            bad_cells.append((row, 3))
        if not np.isnan(r2_c) and r2_c < r2_thr:
            bad_cells.append((row, 4))
    # top-anchored (aligns with the PSD); the channel column is narrower so the value columns breathe
    tb = ax.table(cellText=table, loc='upper center', cellLoc='center',
                  colWidths=[0.15, 0.2125, 0.2125, 0.2125, 0.2125])
    tb.auto_set_font_size(False)
    tb.set_fontsize(8 if len(show_rows) <= 15 else 7)
    tb.scale(1, 1.4 if len(show_rows) <= 15 else 1.05)   # tighter rows once the table gets long
    for col in range(len(header)):                       # header in bold
        tb[(0, col)].get_text().set_fontweight('bold')
    for (r, cc) in bad_cells:                            # breached values in red bold
        cell = tb[(r, cc)]
        cell.get_text().set_color('#c0392b')
        cell.get_text().set_fontweight('bold')
    # Subset note on a SECOND line: appended horizontally it would collide with the PSD title.
    subset_txt = ('' if len(show_rows) == n_ch
                  else f'\n{len(show_rows)}/{n_ch} channels shown ({len(breached)} flagged)')
    ax.set_title(f'Per-channel metrics (stage {stg}) — red = over threshold{subset_txt}', fontsize=9)

    # (c) mean band power + 50 Hz ratio (worst channel — line noise is per channel, so report the
    #     worst one by name instead of tying this panel to the spectrogram's channel selector)
    ax = ax_band
    bp = band_powers(freqs, psds_uV2[ei].mean(axis=0))
    ax.bar([b[0] for b in BANDS], [bp[b[0]] for b in BANDS], color='#807dba')
    ax.set_ylabel('Band power (µV²)')
    ax.set_yscale('log')
    ax.set_title('Mean band power (across channels)', fontsize=9)
    if sf / 2 > 52:
        worst_ratio, worst_ch = np.nan, ''
        for c in range(n_ch):
            f2, p2 = _welch(data[ei, c], fs=sf, nperseg=min(int(4 * sf), data.shape[-1]))
            def _bp(lo, hi):
                m = (f2 >= lo) & (f2 < hi)
                return float(np.trapz(p2[m], f2[m])) if m.sum() > 1 else np.nan
            nb = _bp(40, 47)
            r = _bp(48, 52) / nb if nb and nb > 0 else np.nan
            if not np.isnan(r) and (np.isnan(worst_ratio) or r > worst_ratio):
                worst_ratio, worst_ch = r, ch[c]
        txt = (f'50 Hz / 40–47 Hz = {worst_ratio:.2f}  (worst: {worst_ch})'
               if not np.isnan(worst_ratio) else '50 Hz ratio n/a')
        ax.text(0.98, 0.96, txt, transform=ax.transAxes, ha='right', va='top', fontsize=8)
    else:
        ax.text(0.98, 0.96, '50 Hz ratio n/a (Nyquist ≤ 52 Hz)',
                transform=ax.transAxes, ha='right', va='top', fontsize=8)

    # (d) epoch spectrogram of the selected channel — separates a brief transient artefact (vertical
    #     smear) from a sustained contamination (horizontal band, e.g. line noise).
    ax = ax_spec
    f3, t3, Sxx = sp_spectrogram(data[ei, spectro_ch_idx], fs=sf,
                                 nperseg=int(1.5 * sf), noverlap=int(0.75 * sf))
    fm = f3 <= 40
    ax.pcolormesh(t3, f3[fm], 10 * np.log10(Sxx[fm] + 1e-12), shading='auto', cmap='viridis')
    ax.set_ylabel('Frequency (Hz)')
    ax.set_xlabel('Time (s)')
    ax.set_title(f'Epoch spectrogram — {ch[spectro_ch_idx]}', fontsize=9)

    # (e) spatial view — only a high-density montage with real electrode positions can show this, and
    #     it answers the question that dominates at 32-64 channels: ONE bad electrode (a hot spot) or
    #     a whole-head artefact (a global shift)? Red rings = channels breaching any threshold.
    if topomap:
        mask = np.array([c in breached for c in range(n_ch)], dtype=bool)
        _draw_channel_topomap(ax_topo_ptp, P, ptp_all, 'Peak-to-peak (µV)', cmap='viridis', mask=mask)
        _draw_channel_topomap(ax_topo_mae, P, mae_all, '1/f MAE', cmap='magma', mask=mask)
    return fig


def plot_review_strip(P, final_reject, overridden, custom_stages=(), in_scope=None):
    """Section-4 review strip: hypnogram step-line on top; below, a per-epoch bar coloured by the
    final keep/reject decision, with overridden epochs outlined. Returns a Figure.
    `in_scope` (optional (n_ep,) bool): epochs whose stage is out of the selected scope are drawn GREY
    (excluded — not written to the clean-epo), not green; default None treats every epoch as in scope."""
    stages = P['stages']
    n_ep = len(stages)
    stage_y, stage_colors, ytick_pos, ytick_labels = custom_stage_style(custom_stages)
    floor = min(ytick_pos) - 1
    hypno_y = np.array([stage_y.get(s, floor) for s in stages])
    fig, axes = plt.subplots(2, 1, figsize=(12, 3.2),
                             gridspec_kw={'height_ratios': [1.4, 1.0], 'hspace': 0.15})
    ax = axes[0]
    xs = np.arange(n_ep + 1)
    ax.step(xs, np.append(hypno_y, hypno_y[-1]), where='post', color='#555555', lw=1.0)
    for ei in range(n_ep):
        if stages[ei] == 'R':
            ax.plot([ei, ei + 1], [stage_y['R'], stage_y['R']], color='#c0392b', lw=2.0,
                    solid_capstyle='butt')
    ax.set_yticks(ytick_pos)
    ax.set_yticklabels(ytick_labels, fontsize=7)
    ax.set_xlim(0, n_ep)
    ax.set_xticks([])
    for sp in ['top', 'right', 'bottom']:
        ax.spines[sp].set_visible(False)
    ax = axes[1]
    fr = np.asarray(final_reject, dtype=bool)
    insc = np.ones(n_ep, dtype=bool) if in_scope is None else np.asarray(in_scope, dtype=bool)
    idx = np.arange(n_ep)
    excl = ~insc                       # out-of-scope stages: excluded from the clean-epo -> grey
    keep = insc & ~fr
    rej = insc & fr
    ax.bar(idx[excl], np.ones(int(excl.sum())), width=1.0, color='#bdbdbd', align='edge')
    ax.bar(idx[keep], np.ones(int(keep.sum())), width=1.0, color='#2ecc71', align='edge')
    ax.bar(idx[rej], np.ones(int(rej.sum())), width=1.0, color='#c0392b', align='edge')
    for ei in sorted(overridden):
        ax.plot([ei + 0.5], [1.15], marker='v', color='k', ms=4)
    ax.set_ylim(0, 1.35)
    ax.set_yticks([])
    ax.set_xlim(0, n_ep)
    ax.set_xlabel('Epoch index', fontsize=9)
    ax.set_title(f'Final decision — green = keep ({int(keep.sum())}), red = reject ({int(rej.sum())}), '
                 f'grey = excluded/out-of-scope ({int(excl.sum())}), ▼ = overridden ({len(overridden)})',
                 fontsize=9)
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    # subplots_adjust (not tight_layout): the spine-stripped strip axes are not tight_layout-compatible.
    fig.subplots_adjust(left=0.05, right=0.99, top=0.90, bottom=0.16, hspace=0.15)
    return fig


def build_manual_decision_row(file_id, P, base_reject, final_reject, in_scope, overridden,
                              stages_sel, methods_sel, event_types_sel, thresholds, epoch_length_s=30,
                              dropped_channels=(), epoch_rule='any', epoch_rule_value=20.0,
                              flags_source=''):
    """One-row durable record of a tool-7 manual review — the analogue of tool 7bis's
    `{file_id}_autoreject_decision.tsv`, and the source rebuilt into the global summary. Counts are over
    the IN-SCOPE epochs (the selected stages, i.e. exactly what the clean-epo contains).
    The channel-triage provenance (`dropped_channels`, `epoch_rule`, …) is written as ADDITIVE columns
    whose values are constant on the classic path (no channel dropped, rule 'any').
    Returns a one-row DataFrame."""
    base = np.asarray(base_reject, dtype=bool)
    fin = np.asarray(final_reject, dtype=bool)
    insc = np.asarray(in_scope, dtype=bool)
    stages = np.asarray(P['stages'])
    n_in = int(insc.sum())
    n_rej = int((insc & fin).sum())
    amp = thresholds.get('amplitude_ptp_uV', {})
    row = {
        'file_id': file_id,
        'n_epochs': int(len(stages)),
        'n_in_scope': n_in,
        'n_out_of_scope': int((~insc).sum()),
        'n_rejected': n_rej,
        'n_kept': int((insc & ~fin).sum()),
        'pct_rejected': round(100.0 * n_rej / n_in, 2) if n_in else np.nan,
        'n_flagged_by_selection': int((insc & base).sum()),
        'n_overrides': len(overridden),
        'n_rescued': int((insc & base & ~fin).sum()),
        'n_added': int((insc & ~base & fin).sum()),
        'stages_used': '+'.join(stages_sel),
        'methods_used': '+'.join(methods_sel),
        'event_types_used': '+'.join(event_types_sel) if event_types_sel else '',
        'n_channels': len(P['ch_names']),
        'channels': '+'.join(P['ch_names']),
        'n_channels_dropped': len(dropped_channels or ()),
        'dropped_channels': '+'.join(dropped_channels or ()),
        'epoch_rule': epoch_rule,
        'epoch_rule_value': epoch_rule_value if epoch_rule != 'any' else '',
        'channel_flags_source': flags_source,
        'sfreq': float(P['sfreq']),
        'epoch_length_s': epoch_length_s,
        'thr_amplitude_ptp_uV': '|'.join(f'{k}:{v:g}' for k, v in sorted(amp.items())),
        'thr_flat_ptp_uV': thresholds.get('flat_ptp_uV'),
        'thr_gradient_uV_per_sample': thresholds.get('gradient_uV_per_sample'),
        'thr_1f_mae_max': thresholds.get('1f_mae_max'),
        'thr_1f_r2_min': thresholds.get('1f_r2_min'),
        'reviewed_at': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
    }
    # per-stage rejected counts over the in-scope stages (one column per selected stage)
    for st in stages_sel:
        sel = insc & (stages == st)
        n = int(sel.sum())
        row[f'n_{st}'] = n
        row[f'n_rejected_{st}'] = int((sel & fin).sum())
    return pd.DataFrame([row])


def manual_decision_html(file_id, decision_row, stage_table_html=''):
    """Human-readable recap of a manual review for the report: headline counts + the parameters used,
    followed by the per-stage x per-method table. Returns an HTML string."""
    r = decision_row.iloc[0]
    def _cell(label, value):
        return (f'<tr><td style="padding:3px 10px;border:1px solid #ccc;">{label}</td>'
                f'<td style="padding:3px 10px;border:1px solid #ccc;text-align:right;"><b>{value}</b></td></tr>')
    pct = '—' if pd.isna(r['pct_rejected']) else f'{r["pct_rejected"]:.1f}%'
    body = (_cell('Epochs in file', r['n_epochs'])
            + _cell('In scope (selected stages)', r['n_in_scope'])
            + _cell('Excluded (out-of-scope stages)', r['n_out_of_scope'])
            + _cell('Rejected', f'{r["n_rejected"]}  ({pct} of in-scope)')
            + _cell('Kept &rarr; clean-epo', r['n_kept'])
            + _cell('Flagged by the selection', r['n_flagged_by_selection'])
            + _cell('Manual overrides', f'{r["n_overrides"]} '
                                        f'(rescued {r["n_rescued"]}, newly rejected {r["n_added"]})')
            + _cell('Stages used', r['stages_used'])
            + _cell('Methods used', r['methods_used'])
            + _cell('Event types used', r['event_types_used'] or '—')
            + _cell('Channels', f'{r["n_channels"]} ({r["channels"]})')
            + _cell('Channels dropped', (f'{r["n_channels_dropped"]} ({r["dropped_channels"]})'
                                         if r.get('n_channels_dropped') else '—'))
            + _cell('Epoch rule', (f'{r["epoch_rule"]} {r["epoch_rule_value"]}'
                                   if r.get('epoch_rule', 'any') != 'any'
                                   else 'any flagged channel'))
            + _cell('Epoch length', f'{r["epoch_length_s"]} s')
            + _cell('Reviewed at', r['reviewed_at']))
    return (f'<h3>{file_id} — manual epoch rejection</h3>'
            f'<table style="border-collapse:collapse;font-size:.9em;">{body}</table>'
            + (stage_table_html or ''))


def save_report_html(out_path, title, figs, table_html):
    """Write an mne.Report HTML with the Section-2 figures + rejection table."""
    report = mne.Report(title=title, verbose=False)
    report.add_html(html=table_html, title='Rejection summary', section='Overview')
    for name, fig in figs:
        report.add_figure(fig=fig, title=name, section='Per-stage', image_format='PNG')
        plt.close(fig)
    report.save(str(out_path), overwrite=True, open_browser=False, verbose=False)
