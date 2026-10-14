"""Spectral features of tool 9 (9_spectral_features_voila), shared with tool 0.

Moved verbatim from tool 9 (its setup constants, its computation and plotting cells, and the computation lines
of its run cell): tool 9 and tool 0's step 8 call the same code, so they give the same features by
construction. Edit here, never in a notebook.

PSD in uV^2/Hz (Welch or multitaper with normalization='full'), the specparam aperiodic fit on a smoothed copy
(the band powers always use the raw PSD), band integrals as rectangular sums, the peak model at specparam's
defaults (deliberately not tool 7's: see SPEC section 9).
"""
import os
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import mne

import participant_selection_lib as PSL     # the excluded-participant banner of the report (same folder)
from edf_io_lib import (   # one copy, in edf_io_lib (same folder)
    BASE_STAGE_COLORS, CUSTOM_STAGE_PALETTE, load_custom_stages, parse_custom_field,
)

# numpy 2 renamed trapz -> trapezoid; keep both working
_trapz = getattr(np, 'trapezoid', np.trapz)

# --- optional dependencies (all probed, never fatal at import time) -----------
try:
    from specparam import SpectralModel
    HAS_SPECPARAM = True
except Exception:
    HAS_SPECPARAM = False

try:
    from statsmodels.nonparametric.smoothers_lowess import lowess
    HAS_LOWESS = True
except Exception:
    HAS_LOWESS = False

try:
    from joblib import Parallel, delayed
    HAS_JOBLIB = True
except Exception:
    HAS_JOBLIB = False

# --- sleep stages ------------------------------------------------------------
AASM_STAGES = ['W', 'N1', 'N2', 'N3', 'R']

# --- frequency bands ---------------------------------------------------------
DEFAULT_BANDS = [('delta', 0.5, 4.0), ('theta', 4.0, 8.0), ('alpha', 8.0, 12.0),
                 ('sigma', 12.0, 16.0), ('beta', 16.0, 30.0), ('gamma', 30.0, 45.0)]

# --- exportable band-power measures ------------------------------------------
# (key, label, one-line description shown in the UI, default state, needs the 1/f fit)
MEASURES = [
    ('power_mean_uV2_Hz', 'Mean power density',
     'Mean of the PSD over the band (uV^2/Hz) - independent of the band width.',
     True, False),
    ('power_abs_uV2', 'Absolute power',
     'Integral of the PSD over the band (uV^2) - scales with the band width.',
     True, False),
    ('power_db', 'Power in dB',
     '10*log10 of the mean power density (dB/Hz) - the usual scale for statistics.',
     True, False),
    ('power_rel', 'Relative power',
     'Band absolute power divided by the absolute power of the whole PSD range - '
     'dimensionless, all bands sum to ~1.',
     True, False),
    ('power_ap_removed_db', 'Power above the 1/f fit (dB)',
     'Mean over the band of (PSD in dB - aperiodic fit in dB): how far the spectrum sits '
     'above the 1/f floor. Removes offset differences between subjects and stages.',
     True, True),
    ('power_ratio_over_aperiodic', 'Fold-change over the 1/f fit',
     'Mean over the band of (PSD / aperiodic fit) in linear space - a fold-change '
     '(2 = twice the 1/f floor). Always positive; not a re-expression of the dB measure '
     '(the mean of a ratio is not the ratio of the means).',
     False, True),
    ('power_rel_periodic', 'Relative periodic power',
     'Residual power max(PSD - aperiodic fit, 0) integrated over the band, divided by the same '
     'over the whole range. Like the two measures above it rests on the specparam fit, so it is '
     'only as good as that fit; the clipping at 0 makes it the most fragile of the three when the '
     'fit is poor - check the fit_ok rate in the report before using it.',
     False, True),
]

MEASURE_NEEDS_FIT = {k: needs for k, _, _, _, needs in MEASURES}

# Measures whose epoch -> stage average is reported in BOTH averaging spaces
# (see the note in Section 2): the dB ones. Everything else gets a plain mean.
DB_MEASURES = ['power_db', 'power_ap_removed_db']

# --- defaults ----------------------------------------------------------------
DEF_PSD = dict(method='welch', fmin=0.5, fmax=45.0, win_s=4.0, overlap_pct=50.0,
               window='hann', bandwidth=2.0)

# Peak model kept close to the specparam defaults (peak_width_limits=(0.5, 12), min_peak_height=0)
# rather than to tool 7's wider settings: a 20 Hz-wide "peak" can absorb part of the aperiodic
# slope, and a 0.3 threshold discards small oscillations whose power then leaks into the fit.
DEF_FIT = dict(fit_fmin=2.0, fit_fmax=45.0, aperiodic_mode='fixed',
               peak_width_min=0.5, peak_width_max=12.0, min_peak_height=0.1,
               max_n_peaks=8, peak_threshold=2.0, r2_min=0.90, mae_max=0.15)

DEF_SMOOTH = dict(enabled=True, lowess=True, median=True,
                  median_span_hz=3.0, lowess_span_hz=2.0)

DEF_MIN_EPOCHS = 20

# Below this many (epoch x channel) fits the serial loop beats the parallel one (worker
# start-up costs ~3 s on Windows, one fit ~3 ms) - measured, see fit_spectra.
PARALLEL_MIN_FITS = 5000

# A knee this large means specparam found no real bend and the parameter ran away.
KNEE_DEGENERATE = 1e4

def stage_style(custom_stages):
    """Stage -> y position and colour for the hypnogram strips (AASM + custom below N3).
    Duplicated from qc_rejected_epochs_lib.custom_stage_style (keep in sync)."""
    stage_y = {'W': 4, 'R': 3, 'N1': 2, 'N2': 1, 'N3': 0}
    colors = dict(BASE_STAGE_COLORS)
    for i, cs in enumerate(custom_stages):
        stage_y[cs] = -1 - i
        colors[cs] = CUSTOM_STAGE_PALETTE[i % len(CUSTOM_STAGE_PALETTE)]
    ordered = sorted(stage_y.items(), key=lambda kv: kv[1])
    return stage_y, colors, [y for _, y in ordered], ['REM' if s == 'R' else s for s, _ in ordered]

def norm(p):
    """Normalised path/id string, used on BOTH sides of every comparison (Windows safety)."""
    return os.path.normcase(str(p))


# =============================================================================
# Core computation - discovery, PSD, specparam fit, band measures, aggregation
# (pure functions: no widget is touched here)
# =============================================================================

# ---------------------------------------------------------------------------
# Discovery and provenance
# ---------------------------------------------------------------------------
def find_participants(root):
    """List the participants of the chosen preprocessed-epochs folder: tool 7's `*_all-epo.fif` (every epoch,
    the deselected channels included as bad). Which epochs / channels count is read from the decision
    tables of tool 8 / 8bis (see `load_decision`), never from a separate cleaned copy of the epochs."""
    root = Path(root)
    if not root.exists():
        return []
    parts = []
    for f in sorted(root.rglob('*_all-epo.fif'), key=lambda p: p.name):
        parts.append({'file_id': f.name[:-len('_all-epo.fif')], 'fif': f, 'folder': f.parent,
                      'subtree': f.parent.relative_to(root)})
    parts.sort(key=lambda d: d['file_id'])
    return parts


def detect_source(decision_root):
    """Label the origin of the keep/reject decision, stored in every output as the `source` column:
    manual (tool 8) / auto (tool 8bis) / raw (no decision: every epoch kept, or the tool-7 flags dropped
    when 'Exclude the epochs tool 7 flagged' is ticked -> raw_tool7flags, set at run time)."""
    if decision_root is None:
        return 'raw'
    name = norm(Path(decision_root).name)
    if 'manual' in name:
        return 'manual'
    if 'auto' in name:
        return 'auto'
    return 'decision'


def load_decision(p, meta, ch_names, bads, source):
    """(in_scope, rejected, dropped channels, note) for the epochs of `meta` (row order).

    With a decision folder: joined on epoch_idx from {id}_epoch_decision.tsv; the dropped channels come
    from {id}_channel_decision.tsv plus the channels tool 7 kept as bad. Without one: every epoch in
    scope, none rejected (or the tool-7 flags when source == 'raw_tool7flags'), the bad channels dropped.
    An epoch missing from the decision table is treated as out of scope (and reported)."""
    epoch_table = channel_table = None
    if source not in ('raw', 'raw_tool7flags'):
        if p.get('epoch_decision') is None:
            raise ValueError('no _epoch_decision.tsv for this participant in the decision folder')
        epoch_table = pd.read_csv(p['epoch_decision'], sep='\t')
        if p.get('channel_decision') is not None:
            channel_table = pd.read_csv(p['channel_decision'], sep='\t')
    return decision_from_tables(meta, ch_names, bads, source, epoch_table, channel_table)


def decision_from_tables(meta, ch_names, bads, source, epoch_table=None, channel_table=None):
    """load_decision on tables already in memory (tool 0 passes its step 7's): `epoch_table` /
    `channel_table` = the {id}_epoch_decision.tsv / _channel_decision.tsv tables (None without a decision)."""
    n = len(meta)
    epoch_idx = meta['epoch_idx'].astype(int).values
    dropped = set(bads)
    note = ''
    if source in ('raw', 'raw_tool7flags'):
        in_scope = np.ones(n, dtype=bool)
        if source == 'raw_tool7flags' and 'reject_flag' in meta.columns:
            rejected = meta['reject_flag'].astype(bool).values
        else:
            rejected = np.zeros(n, dtype=bool)
    else:
        ed = epoch_table.drop_duplicates('epoch_idx').set_index('epoch_idx')
        known = np.isin(epoch_idx, ed.index.values)
        in_scope = ed['in_scope'].reindex(epoch_idx).fillna(False).astype(bool).values
        rejected = ed['rejected'].reindex(epoch_idx).fillna(False).astype(bool).values & in_scope
        if not known.all():
            note = (f'{int((~known).sum())} epoch(s) missing from the decision table: treated as out of '
                    f'scope (decision older than the tool-7 run?)')
        if channel_table is not None:
            cd = channel_table
            dropped |= set(cd.loc[cd['dropped'].astype(bool), 'channel'].astype(str))
    return in_scope, rejected, [c for c in ch_names if c in dropped], note

def find_sidecar(file_id, fif_folder, clean_root, data_root):
    """Locate tool 7's {file_id}_preprocessing_params.json - it lives in `preproc_epo/`, not next
    to the clean `.fif`. Searched beside the fif, then under the clean folder's parent
    (i.e. `derivatives/`, which holds `preproc_epo*/`), then under the study folder.
    Returns a Path or None (absence is never fatal - the sidecar is provenance only)."""
    name = f'{file_id}_preprocessing_params.json'
    cand = Path(fif_folder) / name
    if cand.exists():
        return cand
    for root in [Path(clean_root).parent if clean_root else None,
                 Path(data_root) if data_root else None]:
        if root is None or not root.exists():
            continue
        try:
            hits = sorted(root.rglob(name))
        except Exception:
            hits = []
        if hits:
            return hits[0]
    return None


def load_sidecar(path):
    """Read tool 7's params sidecar - provenance only (filters, epoch length, PSD smoothing).
    Missing / unreadable -> found=False and every field None (non-fatal)."""
    info = {'found': False, 'epoch_length_s': None, 'filter': None, 'resample': None,
            'notch': None, 'psd_smoothing': None, 'methods_run': None,
            'custom_stages': [], 'raw': {}}
    if path is None or not Path(path).exists():
        return info
    try:
        with open(path, encoding='utf-8') as f:
            data = json.load(f)
        info['found'] = True
        info['raw'] = data
        info['epoch_length_s'] = data.get('epoch_length_s')
        info['filter'] = data.get('filter')
        info['resample'] = data.get('resample')
        info['notch'] = data.get('notch')
        info['psd_smoothing'] = data.get('psd_smoothing')
        info['methods_run'] = data.get('methods_run')
        info['custom_stages'] = list(data.get('custom_stages', []))
    except Exception:
        pass
    return info


def fit_settings_diff(cfg, side):
    """1/f settings that differ between THIS run and the tool-7 run that flagged the epochs.

    Tool 7 already rejects epochs on a 1/f criterion, so a low `fit_ok` rate here is worth
    explaining: most of the time it comes from the fit being recomputed with different settings
    (a different fit range, different thresholds, or smoothing that was added to tool 7 later)."""
    if not side.get('found'):
        return []
    rt = (side.get('raw') or {}).get('rejection_thresholds', {}) or {}
    fit = cfg['fit']
    diffs = []
    fr = rt.get('1f_fit_range_hz')
    if isinstance(fr, (list, tuple)) and len(fr) == 2:
        if (float(fr[0]), float(fr[1])) != (fit['fit_fmin'], fit['fit_fmax']):
            diffs.append(f"fit range {float(fr[0]):g}-{float(fr[1]):g} Hz in tool 7, "
                         f"{fit['fit_fmin']:g}-{fit['fit_fmax']:g} Hz here")
    if '1f_r2_min' in rt and float(rt['1f_r2_min']) != fit['r2_min']:
        diffs.append(f"R2 bound {float(rt['1f_r2_min']):g} in tool 7, {fit['r2_min']:g} here")
    if '1f_mae_max' in rt and float(rt['1f_mae_max']) != fit['mae_max']:
        diffs.append(f"MAE bound {float(rt['1f_mae_max']):g} in tool 7, {fit['mae_max']:g} here")
    t6_smooth = bool((side.get('psd_smoothing') or {}).get('enabled'))
    if t6_smooth != bool(cfg['smooth']['enabled']):
        diffs.append(f"PSD smoothing {'on' if t6_smooth else 'off'} in tool 7, "
                     f"{'on' if cfg['smooth']['enabled'] else 'off'} here")
    return diffs


def load_subject_info(path):
    """Read the optional participant-info table (.csv/.tsv). Returns (DataFrame, error)."""
    if not path:
        return None, ''
    try:
        sep = '\t' if norm(path).endswith('.tsv') else None
        df = pd.read_csv(path, sep=sep, engine='python')
        return df, ''
    except Exception as exc:
        return None, str(exc)


def assign_thirds(epoch_idx, stages):
    """Split the sleep period into three equal spans of epoch index (draft feature).

    The sleep period runs from the first to the last epoch scored as anything but W; it is
    divided into three equal spans of the ORIGINAL epoch index, so the thirds follow real
    time even though the rejected epochs are missing from the file. Epochs outside are
    labelled `pre` / `post`. NREM-REM cycle detection is deliberately not attempted."""
    epoch_idx = np.asarray(epoch_idx)
    stages = np.asarray(stages)
    out = np.array(['none'] * len(stages), dtype=object)
    sleep = np.where(stages != 'W')[0]
    if len(sleep) == 0:
        return out
    lo, hi = float(epoch_idx[sleep[0]]), float(epoch_idx[sleep[-1]])
    span = max((hi - lo + 1) / 3.0, 1e-9)
    for i, e in enumerate(epoch_idx):
        if e < lo:
            out[i] = 'pre'
        elif e > hi:
            out[i] = 'post'
        else:
            out[i] = 'T%d' % (min(2, int((e - lo) // span)) + 1)
    return out


# ---------------------------------------------------------------------------
# PSD
# ---------------------------------------------------------------------------
def compute_psd_array(epochs, method='welch', fmin=0.5, fmax=45.0, win_s=4.0,
                      overlap_pct=50.0, window='hann', bandwidth=2.0):
    """Power spectral density per (epoch, channel), in uV^2/Hz.

    Welch: sliding window of `win_s` seconds with `overlap_pct` overlap (4 s -> 0.25 Hz
    resolution, the sleep-research convention; the rule of thumb is a resolution fine enough
    for ~2 cycles of the lowest frequency of interest). Multitaper: `bandwidth` Hz smoothing.
    fmax is capped at the Nyquist frequency and the window at the epoch length; both are
    reported back in `notes` so the run log can warn. Returns (freqs, psds_uV2, notes)."""
    notes = []
    sf = float(epochs.info['sfreq'])
    epoch_len_s = len(epochs.times) / sf
    fmax_eff = min(float(fmax), sf / 2 - 0.5)
    if fmax_eff < fmax:
        notes.append(f'fmax capped to {fmax_eff:.1f} Hz (Nyquist at sfreq={sf:g} Hz)')
    if method == 'multitaper':
        # normalization='full' is REQUIRED: MNE defaults to 'length', which is not a power
        # density - it is larger by a factor of sfreq (x256 = +24 dB at 256 Hz) and, worse,
        # the factor follows the sampling rate, so recordings at 256 and 512 Hz would not be
        # comparable. With 'full' the multitaper PSD matches Welch to within ~0.6 dB.
        spec = epochs.compute_psd(method='multitaper', fmin=fmin, fmax=fmax_eff,
                                  bandwidth=bandwidth, normalization='full', verbose=False)
    else:
        win_eff = min(float(win_s), epoch_len_s)
        if win_eff < win_s:
            notes.append(f'Welch window clipped to the epoch length ({win_eff:g} s)')
        n_per_seg = max(int(round(win_eff * sf)), 8)
        n_overlap = int(n_per_seg * float(overlap_pct) / 100.0)
        spec = epochs.compute_psd(method='welch', fmin=fmin, fmax=fmax_eff,
                                  n_fft=n_per_seg, n_overlap=n_overlap, n_per_seg=n_per_seg,
                                  window=window, verbose=False)
    return spec.freqs, spec.get_data() * 1e12, notes      # V^2/Hz -> uV^2/Hz


# --- PSD smoothing (duplicated VERBATIM from 7_preprocessing_voila / -------------
# --- qc_rejected_epochs_lib - keep in sync). Applied to the FIT COPY only. -------
def smooth_psd_median(psds_uV2, freqs, span_hz=3.0):
    """Running-median smooth of each (epoch, channel) linear PSD along frequency.
    Reproduces oscip.smooth_spectrum_median (MATLAB movmedian over a span given in Hz):
    a centred median filter that removes narrow spikes (residual line-noise harmonics,
    single-bin artefacts) BEFORE the LOWESS mean smoothing."""
    freqs = np.asarray(freqs)
    if len(freqs) < 4:
        return psds_uV2
    freq_res = float(np.median(np.diff(freqs)))
    span_pts = max(3, int(round(span_hz / freq_res)))
    out = np.array(psds_uV2, dtype=float, copy=True)
    n_ep, n_ch, _ = out.shape
    for ei in range(n_ep):
        for ci in range(n_ch):
            try:
                out[ei, ci] = pd.Series(out[ei, ci]).rolling(
                    window=span_pts, center=True, min_periods=1).median().to_numpy()
            except Exception:
                pass
    return out


def smooth_psd_lowess(psds_uV2, freqs, span_hz=2.0):
    """LOWESS-smooth each (epoch, channel) linear PSD along frequency.
    Reproduces oscip.smooth_spectrum (MATLAB non-robust 'lowess', local linear regression)
    over a span given in Hz, applied to the LINEAR power. Non-fatal per spectrum."""
    freqs = np.asarray(freqs)
    if len(freqs) < 4 or not HAS_LOWESS:
        return psds_uV2
    freq_res = float(np.median(np.diff(freqs)))
    span_pts = max(3, int(round(span_hz / freq_res)))
    frac = min(1.0, span_pts / len(freqs))
    out = np.array(psds_uV2, dtype=float, copy=True)
    n_ep, n_ch, _ = out.shape
    for ei in range(n_ep):
        for ci in range(n_ch):
            try:
                out[ei, ci] = lowess(out[ei, ci], freqs, frac=frac, it=0, return_sorted=False)
            except Exception:
                pass
    return out


def smooth_for_fit(psds_uV2, freqs, smooth_cfg):
    """Apply the tool-7 smoothing (median first, then LOWESS - the oscip order) to the copy that
    feeds specparam. Each filter can be used alone. The band powers are ALWAYS computed on the
    unsmoothed PSD."""
    if not smooth_cfg.get('enabled'):
        return psds_uV2
    out = psds_uV2
    if smooth_cfg.get('median'):
        out = smooth_psd_median(out, freqs, span_hz=float(smooth_cfg.get('median_span_hz', 3.0)))
    if smooth_cfg.get('lowess', True):
        out = smooth_psd_lowess(out, freqs, span_hz=float(smooth_cfg.get('lowess_span_hz', 2.0)))
    return out


# ---------------------------------------------------------------------------
# Aperiodic / periodic decomposition (specparam)
# ---------------------------------------------------------------------------
def aperiodic_log10(freqs, mode, params):
    """Aperiodic component in log10(power), evaluated on ANY frequency axis.

    Evaluating the model analytically (rather than reading specparam's fitted curve) lets the
    aperiodic component be extrapolated below the fit's lower bound, so the delta band can be
    corrected too - at the cost of trusting the model outside the fitted range."""
    freqs = np.asarray(freqs, dtype=float)
    if mode == 'knee':
        offset, knee, exponent = params[0], params[1], params[2]
        return offset - np.log10(np.maximum(knee + freqs ** exponent, 1e-30))
    offset, exponent = params[0], params[1]
    return offset - exponent * np.log10(freqs)


def _fit_one_epoch(freqs_fit, psd_ep_fit, cfg, ei=0, peaks=None):
    """specparam fit of every channel of ONE epoch.
    psd_ep_fit: (n_ch, n_freqs_fit). Returns the (n_ch, 6) parameter array and appends
    (epoch, channel, peak_idx, n_peaks, cf, pw, bw) rows to `peaks`."""
    n_ch = psd_ep_fit.shape[0]
    out = np.full((n_ch, 6), np.nan)          # offset, knee, exponent, r2, mae, n_peaks
    if peaks is None:
        peaks = []
    if not HAS_SPECPARAM:
        return out, peaks
    for ci in range(n_ch):
        psd = psd_ep_fit[ci]
        if not np.all(np.isfinite(psd)) or np.any(psd <= 0):
            continue                          # dead / clipped spectrum: leave NaN
        try:
            sm = SpectralModel(peak_width_limits=[cfg['peak_width_min'], cfg['peak_width_max']],
                               aperiodic_mode=cfg['aperiodic_mode'],
                               min_peak_height=cfg['min_peak_height'],
                               max_n_peaks=int(cfg['max_n_peaks']),
                               peak_threshold=cfg['peak_threshold'], verbose=False)
            sm.fit(freqs_fit, psd)
            ap = np.asarray(sm.get_params('aperiodic'), dtype=float).ravel()
            if cfg['aperiodic_mode'] == 'knee':
                offset, knee, exponent = ap[0], ap[1], ap[2]
            else:
                offset, knee, exponent = ap[0], np.nan, ap[1]
            pk = np.asarray(sm.get_params('peak'), dtype=float)
            pk = pk.reshape(-1, 3) if pk.size else np.zeros((0, 3))
            out[ci] = [offset, knee, exponent,
                       float(sm.get_metrics('gof', 'squared')),
                       float(sm.get_metrics('error', 'mae')), len(pk)]
            for pi in range(len(pk)):
                peaks.append((ei, ci, pi + 1, len(pk), pk[pi, 0], pk[pi, 1], pk[pi, 2]))
        except Exception:
            continue                          # failed fit stays NaN (counted as fit_ok=False)
    return out, peaks


def _fit_chunk(freqs_fit, psd_chunk, cfg, ei0):
    """Fit a contiguous block of epochs - the unit of work handed to each parallel worker.
    Chunking (rather than one task per epoch) is what makes the parallel path actually
    faster: a single epoch is only a few fits, far too little to pay for the dispatch."""
    n_ep = psd_chunk.shape[0]
    arr = np.full((n_ep, psd_chunk.shape[1], 6), np.nan)
    peaks = []
    for k in range(n_ep):
        arr[k], _ = _fit_one_epoch(freqs_fit, psd_chunk[k], cfg, ei=ei0 + k, peaks=peaks)
    return ei0, arr, peaks


def fit_spectra(freqs, psds_uV2, cfg, n_jobs=1, progress=None):
    """specparam fit of every (epoch, channel) spectrum.

    Returns (res, peak_rows) where res holds (n_ep, n_ch) arrays for offset / knee /
    exponent / r2 / mae / n_peaks / fit_ok, and peak_rows is a list of
    (epoch_i, channel_i, peak_idx, n_peaks, center_freq, power, bandwidth).
    The parallel path is exact - it calls the same per-epoch function as the serial one -
    and any failure falls back to the serial loop."""
    n_ep, n_ch, _ = psds_uV2.shape
    fmask = (np.asarray(freqs) >= cfg['fit_fmin']) & (np.asarray(freqs) <= cfg['fit_fmax'])
    freqs_fit = np.asarray(freqs)[fmask]
    psds_fit = psds_uV2[:, :, fmask]
    params = np.full((n_ep, n_ch, 6), np.nan)
    peak_rows = []

    # Spawning workers costs ~3 s on Windows, so the parallel path only pays off on a big
    # job (measured 2.8 ms/fit: below ~5000 fits the serial loop is simply faster).
    done = False
    if n_jobs != 1 and HAS_JOBLIB and n_ep > 8 and n_ep * n_ch >= PARALLEL_MIN_FITS:
        try:
            n_workers = os.cpu_count() or 4 if n_jobs < 0 else n_jobs
            chunk = max(8, int(np.ceil(n_ep / max(1, n_workers))))
            starts = list(range(0, n_ep, chunk))
            results = Parallel(n_jobs=n_jobs)(
                delayed(_fit_chunk)(freqs_fit, psds_fit[s:s + chunk], cfg, s) for s in starts)
            for ei0, arr, peaks in results:
                params[ei0:ei0 + len(arr)] = arr
                peak_rows.extend(peaks)
            done = True
        except Exception:
            params[:] = np.nan
            peak_rows = []
            done = False                      # fall back to the serial loop
    if not done:
        for ei in range(n_ep):
            params[ei], _ = _fit_one_epoch(freqs_fit, psds_fit[ei], cfg, ei=ei, peaks=peak_rows)
            if progress is not None and (ei % 25 == 0 or ei == n_ep - 1):
                progress(ei + 1, n_ep)

    res = {'offset': params[:, :, 0], 'knee': params[:, :, 1], 'exponent': params[:, :, 2],
           'r2': params[:, :, 3], 'mae': params[:, :, 4], 'n_peaks': params[:, :, 5]}
    res['fit_ok'] = (np.isfinite(res['r2']) & np.isfinite(res['mae'])
                     & (res['r2'] >= cfg['r2_min']) & (res['mae'] <= cfg['mae_max']))
    return res, peak_rows


def aperiodic_psd_log10(freqs, res, mode):
    """Aperiodic model evaluated on the full PSD axis -> (n_ep, n_ch, n_freqs) in log10 power."""
    n_ep, n_ch = res['offset'].shape
    out = np.full((n_ep, n_ch, len(freqs)), np.nan)
    for ei in range(n_ep):
        for ci in range(n_ch):
            off, exp_, knee = res['offset'][ei, ci], res['exponent'][ei, ci], res['knee'][ei, ci]
            if not (np.isfinite(off) and np.isfinite(exp_)):
                continue
            if mode == 'knee' and not np.isfinite(knee):
                continue
            out[ei, ci] = aperiodic_log10(freqs, mode, (off, knee, exp_) if mode == 'knee'
                                          else (off, exp_))
    return out


# ---------------------------------------------------------------------------
# Band measures
# ---------------------------------------------------------------------------
def compute_band_measures(freqs, psds_uV2, ap_log10, bands, measures_on):
    """Band power for every selected measure.

    Returns {band_name: {measure: (n_ep, n_ch) array}}. Non-positive PSD bins (dead or
    clipped channels) become NaN before any log, so a bad channel yields NaN rather than -inf.
    Bin selection is half-open [fmin, fmax) so adjacent bands never share a bin, and the
    integral is a rectangular sum (each bin stands for a `df`-wide slab). With a 0.25 Hz
    resolution that gives each band its exact width (delta 0.5-4 Hz -> 14 bins x 0.25 = 3.5 Hz)
    and, when the bands tile the PSD range without gaps, relative powers that sum to 1 -
    which a trapezoidal rule over the band's bins alone would not do."""
    freqs = np.asarray(freqs, dtype=float)
    df = float(np.median(np.diff(freqs))) if len(freqs) > 1 else 1.0
    psd = np.where(psds_uV2 > 0, psds_uV2, np.nan)
    psd_db = 10.0 * np.log10(psd)
    ap_lin = 10.0 ** ap_log10 if ap_log10 is not None else None
    ap_db = 10.0 * ap_log10 if ap_log10 is not None else None
    ratio = psd / ap_lin if ap_lin is not None else None
    residual = np.maximum(psd - ap_lin, 0.0) if ap_lin is not None else None

    total_abs = np.nansum(psd, axis=2) * df
    total_res = np.nansum(residual, axis=2) * df if residual is not None else None

    out = {}
    for name, lo, hi in bands:
        m = (freqs >= lo) & (freqs < hi)
        if m.sum() < 1:
            continue                          # band outside the PSD range: skipped upstream
        vals = {}
        if 'power_mean_uV2_Hz' in measures_on or 'power_db' in measures_on:
            mean_lin = np.nanmean(psd[:, :, m], axis=2)
            if 'power_mean_uV2_Hz' in measures_on:
                vals['power_mean_uV2_Hz'] = mean_lin
            if 'power_db' in measures_on:
                vals['power_db'] = 10.0 * np.log10(np.where(mean_lin > 0, mean_lin, np.nan))
        if 'power_abs_uV2' in measures_on or 'power_rel' in measures_on:
            abs_pow = np.nansum(psd[:, :, m], axis=2) * df
            if 'power_abs_uV2' in measures_on:
                vals['power_abs_uV2'] = abs_pow
            if 'power_rel' in measures_on:
                vals['power_rel'] = abs_pow / np.where(total_abs > 0, total_abs, np.nan)
        if 'power_ap_removed_db' in measures_on and ap_db is not None:
            vals['power_ap_removed_db'] = np.nanmean(psd_db[:, :, m] - ap_db[:, :, m], axis=2)
        if 'power_ratio_over_aperiodic' in measures_on and ratio is not None:
            vals['power_ratio_over_aperiodic'] = np.nanmean(ratio[:, :, m], axis=2)
        if 'power_rel_periodic' in measures_on and residual is not None:
            res_band = np.nansum(residual[:, :, m], axis=2) * df
            vals['power_rel_periodic'] = res_band / np.where(total_res > 0, total_res, np.nan)
        out[name] = vals
    return out


# ---------------------------------------------------------------------------
# Long-format tables (epoch level)
# ---------------------------------------------------------------------------
def build_bandpower_epoch(file_id, source, meta, ch_names, band_vals, measures_on, thirds=None):
    """Long table: one row per (epoch, channel, band)."""
    rows = []
    epoch_idx = meta['epoch_idx'].values
    stages = meta['stage'].astype(str).values
    for band, vals in band_vals.items():
        for ci, ch in enumerate(ch_names):
            for ei in range(len(stages)):
                row = {'file_id': file_id, 'source': source, 'epoch_idx': int(epoch_idx[ei]),
                       'stage': stages[ei], 'channel': ch, 'band': band}
                if thirds is not None:
                    row['third'] = thirds[ei]
                for meas in measures_on:
                    if meas in vals:
                        row[meas] = float(vals[meas][ei, ci])
                rows.append(row)
    return pd.DataFrame(rows)


def build_aperiodic_epoch(file_id, source, meta, ch_names, res, thirds=None):
    """Long table: one row per (epoch, channel) with the aperiodic parameters and fit quality."""
    rows = []
    epoch_idx = meta['epoch_idx'].values
    stages = meta['stage'].astype(str).values
    for ci, ch in enumerate(ch_names):
        for ei in range(len(stages)):
            row = {'file_id': file_id, 'source': source, 'epoch_idx': int(epoch_idx[ei]),
                   'stage': stages[ei], 'channel': ch}
            if thirds is not None:
                row['third'] = thirds[ei]
            row.update({'offset': float(res['offset'][ei, ci]),
                        'knee': float(res['knee'][ei, ci]),
                        'exponent': float(res['exponent'][ei, ci]),
                        'r2': float(res['r2'][ei, ci]),
                        'mae': float(res['mae'][ei, ci]),
                        'n_peaks': float(res['n_peaks'][ei, ci]),
                        'fit_ok': bool(res['fit_ok'][ei, ci])})
            rows.append(row)
    return pd.DataFrame(rows)


def build_peaks_table(file_id, source, meta, ch_names, peak_rows):
    """Long table: one row per specparam peak (centre frequency, power, bandwidth)."""
    epoch_idx = meta['epoch_idx'].values
    stages = meta['stage'].astype(str).values
    rows = []
    for ei, ci, pi, npk, cf, pw, bw in peak_rows:
        rows.append({'file_id': file_id, 'source': source, 'epoch_idx': int(epoch_idx[ei]),
                     'stage': stages[ei], 'channel': ch_names[ci], 'peak_idx': int(pi),
                     'n_peaks': int(npk), 'center_freq_hz': float(cf),
                     'peak_power': float(pw), 'bandwidth_hz': float(bw)})
    return pd.DataFrame(rows, columns=['file_id', 'source', 'epoch_idx', 'stage', 'channel',
                                       'peak_idx', 'n_peaks', 'center_freq_hz',
                                       'peak_power', 'bandwidth_hz'])


# ---------------------------------------------------------------------------
# Aggregation (epoch -> stage)
# ---------------------------------------------------------------------------
def stage_masks(stages, stages_sel):
    """Ordered {stage: boolean epoch mask} for the selected stages.
    Groupings such as NREM or whole-night are deliberately left to the user: they are one
    groupby away in the exported tables, and computing them here doubled the table size."""
    stages = np.asarray(stages)
    return {st: (stages == st) for st in stages_sel}


def _agg_stats(values, min_epochs):
    """mean / median / sd over the epochs of one cell; NaN when too few epochs."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) < max(1, min_epochs):
        return np.nan, np.nan, np.nan
    return float(np.mean(v)), float(np.median(v)), float(np.std(v, ddof=1)) if len(v) > 1 else 0.0


def aggregate_bandpower(file_id, source, band_vals, ch_names, masks, measures_on,
                        min_epochs, spaces=('log',), fit_ok=None, exclude_bad_fit=False):
    """Aggregate the per-epoch band measures to one row per (stage, channel, band).

    `spaces` picks how the dB measures are averaged over the epochs of a stage: 'log' = mean of
    the per-epoch dB (what a per-epoch statistical model sees, robust to one loud epoch), 'lin' =
    dB of the mean in linear space (preserves total power, but a few loud epochs dominate). The
    two are exported side by side when both are asked for, suffixed `_from_log` / `_from_lin`.
    The linear measures get a plain mean/median/sd. Cells with fewer than `min_epochs` epochs are
    written as NaN, with n_epochs kept, so a shortage is visible instead of silent."""
    rows = []
    for label, mask in masks.items():
        for ci, ch in enumerate(ch_names):
            use = mask.copy()
            if exclude_bad_fit and fit_ok is not None:
                use = use & fit_ok[:, ci]
            n_ep = int(use.sum())
            n_ok = int((mask & fit_ok[:, ci]).sum()) if fit_ok is not None else np.nan
            for band, vals in band_vals.items():
                row = {'file_id': file_id, 'source': source, 'stage': label, 'channel': ch,
                       'band': band, 'n_epochs': n_ep, 'n_fit_ok': n_ok}
                for meas in measures_on:
                    if meas not in vals:
                        continue
                    v = vals[meas][use, ci]
                    mean, med, sd = _agg_stats(v, min_epochs)
                    if meas not in DB_MEASURES:
                        row[f'{meas}_mean'] = mean
                        row[f'{meas}_median'] = med
                        row[f'{meas}_sd'] = sd
                        continue
                    if 'log' in spaces:
                        row[f'{meas}_from_log_mean'] = mean
                        row[f'{meas}_from_log_median'] = med
                        row[f'{meas}_from_log_sd'] = sd
                    if 'lin' in spaces:
                        lin = 10.0 ** (np.asarray(v, dtype=float) / 10.0)
                        lm, lmd, _ = _agg_stats(lin, min_epochs)
                        row[f'{meas}_from_lin_mean'] = (10.0 * np.log10(lm)
                                                        if np.isfinite(lm) and lm > 0 else np.nan)
                        row[f'{meas}_from_lin_median'] = (10.0 * np.log10(lmd)
                                                          if np.isfinite(lmd) and lmd > 0 else np.nan)
                rows.append(row)
    return pd.DataFrame(rows)


def aggregate_aperiodic(file_id, source, res, ch_names, masks, min_epochs,
                        exclude_bad_fit=False):
    """One row per (stage, channel) with the aperiodic parameters and the fit quality."""
    rows = []
    for label, mask in masks.items():
        for ci, ch in enumerate(ch_names):
            use = mask.copy()
            if exclude_bad_fit:
                use = use & res['fit_ok'][:, ci]
            row = {'file_id': file_id, 'source': source, 'stage': label, 'channel': ch,
                   'n_epochs': int(use.sum()),
                   'n_fit_ok': int((mask & res['fit_ok'][:, ci]).sum())}
            for key in ('offset', 'exponent', 'knee', 'r2', 'mae', 'n_peaks'):
                mean, med, sd = _agg_stats(res[key][use, ci], min_epochs)
                row[f'{key}_mean'], row[f'{key}_median'], row[f'{key}_sd'] = mean, med, sd
            rows.append(row)
    return pd.DataFrame(rows)


def aggregate_psd(file_id, source, freqs, psds_uV2, ap_log10, ch_names, masks, min_epochs,
                  spaces=('log',)):
    """Mean PSD per (stage, channel, frequency) - lets the figures be redrawn without recomputing
    anything. The dB columns follow the selected averaging space(s); the linear mean is always
    kept, since it is the raw material every other column derives from.

    The `_sd` columns are the standard deviation ACROSS EPOCHS, so the +/- SEM band of the report
    figures can be reproduced from this table alone: sem = sd / sqrt(n_epochs)."""
    psd = np.where(psds_uV2 > 0, psds_uV2, np.nan)
    psd_db = 10.0 * np.log10(psd)
    ap_removed = psd_db - 10.0 * ap_log10 if ap_log10 is not None else None
    rows = []
    for label, mask in masks.items():
        n_ep = int(mask.sum())
        if n_ep < max(1, min_epochs):
            continue
        with np.errstate(invalid='ignore'):
            lin = np.nanmean(psd[mask], axis=0)                 # (n_ch, n_freqs)
            lin_sd = np.nanstd(psd[mask], axis=0, ddof=1)
            log = np.nanmean(psd_db[mask], axis=0)
            log_sd = np.nanstd(psd_db[mask], axis=0, ddof=1)
            apr = np.nanmean(ap_removed[mask], axis=0) if ap_removed is not None else None
            apr_sd = np.nanstd(ap_removed[mask], axis=0, ddof=1) if ap_removed is not None else None
        for ci, ch in enumerate(ch_names):
            for fi, fr in enumerate(freqs):
                row = {'file_id': file_id, 'source': source, 'stage': label, 'channel': ch,
                       'freq_hz': float(fr), 'n_epochs': n_ep,
                       'psd_uV2_Hz_from_lin': float(lin[ci, fi]),
                       'psd_uV2_Hz_sd': float(lin_sd[ci, fi])}
                if 'log' in spaces:
                    row['psd_db_from_log'] = float(log[ci, fi])
                    row['psd_db_from_log_sd'] = float(log_sd[ci, fi])
                if 'lin' in spaces:
                    row['psd_db_from_lin'] = (float(10.0 * np.log10(lin[ci, fi]))
                                              if np.isfinite(lin[ci, fi]) and lin[ci, fi] > 0
                                              else np.nan)
                if apr is not None:
                    row['psd_ap_removed_db'] = float(apr[ci, fi])
                    row['psd_ap_removed_db_sd'] = float(apr_sd[ci, fi])
                rows.append(row)
    return pd.DataFrame(rows)


# =============================================================================
# Figures and HTML tables (participant level and database level)
# =============================================================================
def _fig_width(n_ch):
    return min(5.0 * max(n_ch, 1), 20.0)


def _axes_row(n_ch, height=4.2, sharey=True):
    fig, axes = plt.subplots(1, max(n_ch, 1), figsize=(_fig_width(n_ch), height), sharey=sharey)
    return fig, np.atleast_1d(axes)


def stage_curves(values, masks, min_epochs):
    """Mean + SEM across epochs of a (n_ep, n_ch, n_freqs) array, per stage."""
    out = {}
    for label, m in masks.items():
        n = int(m.sum())
        if n < max(1, min_epochs):
            continue
        sub = values[m]
        with np.errstate(invalid='ignore'):
            mean = np.nanmean(sub, axis=0)
            sem = np.nanstd(sub, axis=0, ddof=1) / np.sqrt(max(n, 1))
        out[label] = (mean, sem, n)
    return out


# --- 1. hypnogram of the epochs that survived the rejection ------------------
def plot_clean_hypnogram(meta, custom_stages, file_id, n_total=None):
    """Stage of every surviving epoch against its ORIGINAL epoch index, so the epochs removed
    by tools 7/8/8bis show up as gaps, plus the clean-epoch count per stage."""
    stage_y, colors, ytick_pos, ytick_lab = stage_style(custom_stages)
    epoch_idx = meta['epoch_idx'].values
    stages = meta['stage'].astype(str).values
    fig, axes = plt.subplots(1, 2, figsize=(15, 3.2), gridspec_kw={'width_ratios': [3, 1]})

    ax = axes[0]
    y = np.array([stage_y.get(s, np.nan) for s in stages], dtype=float)
    ax.scatter(epoch_idx, y, s=6, c=[colors.get(s, '#888888') for s in stages])
    ax.set_yticks(ytick_pos)
    ax.set_yticklabels(ytick_lab)
    ax.set_xlabel('Epoch index in the recording')
    kept = len(stages)
    total = n_total if n_total else (int(epoch_idx.max()) + 1 if len(epoch_idx) else 0)
    ax.set_title(f'{file_id} - clean epochs: {kept} kept out of {total} '
                 f'({100.0 * kept / max(total, 1):.0f} %); gaps = removed epochs')
    ax.grid(alpha=0.2)

    ax = axes[1]
    order = [s for s in AASM_STAGES + list(custom_stages) if s in set(stages)]
    counts = [int((stages == s).sum()) for s in order]
    ax.bar(range(len(order)), counts, color=[colors.get(s, '#888888') for s in order])
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(['REM' if s == 'R' else s for s in order])
    ax.set_ylabel('clean epochs')
    ax.set_title('Clean epochs per stage')
    for i, c in enumerate(counts):
        ax.text(i, c, str(c), ha='center', va='bottom', fontsize=8)
    fig.tight_layout()
    return fig


# --- 2 & 3. mean PSD per stage (raw and aperiodic-removed) -------------------
def _curve_from_log(values_db, mask):
    """Average per-epoch dB values (log space); returns (mean, lo, hi, n) in dB."""
    sub = values_db[mask]
    n = int(mask.sum())
    with np.errstate(invalid='ignore'):
        mean = np.nanmean(sub, axis=0)
        sem = np.nanstd(sub, axis=0, ddof=1) / np.sqrt(max(n, 1))
    return mean, mean - sem, mean + sem, n


def _curve_from_lin(values_lin, mask):
    """Average linear power, then convert to dB; returns (mean, lo, hi, n) in dB."""
    sub = values_lin[mask]
    n = int(mask.sum())
    with np.errstate(invalid='ignore'):
        mean = np.nanmean(sub, axis=0)
        sem = np.nanstd(sub, axis=0, ddof=1) / np.sqrt(max(n, 1))

    def to_db(x):
        return 10.0 * np.log10(np.where(x > 0, x, np.nan))

    return to_db(mean), to_db(mean - sem), to_db(mean + sem), n


def plot_psd_per_stage(freqs, masks, ch_names, custom_stages, min_epochs, title, ylabel,
                       values_db=None, values_lin=None):
    """One subplot per channel; one curve per stage with a +/- SEM band.

    `values_db` holds per-epoch values already in dB (averaged in LOG space); `values_lin` holds
    per-epoch LINEAR power (averaged, then converted to dB). Whichever is available is drawn
    solid; when both are, the second is dashed behind so the gap between the two averaging
    spaces is visible on the figure itself."""
    _, colors, _, _ = stage_style(custom_stages)
    primary = ('log', values_db) if values_db is not None else ('lin', values_lin)
    secondary = ('lin', values_lin) if (values_db is not None and values_lin is not None) else None
    fig, axes = _axes_row(len(ch_names))
    for ci, ch in enumerate(ch_names):
        ax = axes[ci]
        for stage, mask in masks.items():
            if int(mask.sum()) < max(1, min_epochs):
                continue
            col = colors.get(stage, '#888888')
            fn = _curve_from_log if primary[0] == 'log' else _curve_from_lin
            mean, lo, hi, n = fn(primary[1], mask)
            ax.plot(freqs, mean[ci], color=col, label=f'{stage} (n={n})')
            ax.fill_between(freqs, lo[ci], hi[ci], color=col, alpha=0.25)
            if secondary is not None:
                alt, _, _, _ = _curve_from_lin(secondary[1], mask)
                ax.plot(freqs, alt[ci], color=col, ls='--', lw=0.8, alpha=0.7)
        ax.set_title(ch)
        ax.set_xlabel('Frequency (Hz)')
        ax.set_xlim(freqs[0], freqs[-1])
        ax.grid(alpha=0.2)
        if ci == 0:
            ax.set_ylabel(ylabel)
            ax.legend(fontsize=8)
    fig.suptitle(title + (' (solid: log-space average, dashed: linear-space average)'
                          if secondary is not None else ''))
    fig.tight_layout()
    return fig


# --- 4. log-log PSD + aperiodic fit ------------------------------------------
def plot_loglog_with_fit(freqs, psds_uV2, ap_log10, masks, ch_names, custom_stages,
                         min_epochs, fit_range, res=None):
    """Mean PSD per stage in log-log with the mean aperiodic fit dashed on top - the visual
    check of the fit, and of whether a `knee` is needed (a real bend in the straight line)."""
    _, colors, _, _ = stage_style(custom_stages)
    psd_db = 10.0 * np.log10(np.where(psds_uV2 > 0, psds_uV2, np.nan))
    curves = stage_curves(psd_db, masks, min_epochs)
    ap_curves = stage_curves(10.0 * ap_log10, masks, min_epochs) if ap_log10 is not None else {}
    fig, axes = _axes_row(len(ch_names))
    for ci, ch in enumerate(ch_names):
        ax = axes[ci]
        for stage, (mean, _, n) in curves.items():
            col = colors.get(stage, '#888888')
            lab = f'{stage} (n={n})'
            if res is not None:
                exps = res['exponent'][masks[stage], ci]
                if np.isfinite(exps).any():
                    lab += f'  exp={np.nanmean(exps):.2f}'
            ax.plot(freqs, mean[ci], color=col, label=lab)
            if stage in ap_curves:
                ax.plot(freqs, ap_curves[stage][0][ci], color=col, ls='--', lw=1.0, alpha=0.8)
        ax.axvspan(fit_range[0], fit_range[1], color='#000000', alpha=0.04)
        ax.set_xscale('log')
        ax.set_xticks([0.5, 1, 2, 4, 8, 16, 30, 45])
        ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        ax.set_xlim(freqs[0], freqs[-1])
        ax.set_title(ch)
        ax.set_xlabel('Frequency (Hz)')
        ax.grid(alpha=0.2, which='both')
        if ci == 0:
            ax.set_ylabel('PSD (dB/Hz)')
            ax.legend(fontsize=8)
    fig.suptitle('log-log PSD per stage with the aperiodic fit (dashed); '
                 f'shaded = fit range {fit_range[0]:g}-{fit_range[1]:g} Hz')
    fig.tight_layout()
    return fig


# --- 5. aperiodic parameters across the night --------------------------------
def plot_aperiodic_over_night(meta, res, ch_names, custom_stages):
    """Exponent and offset of every epoch against its position in the night, coloured by
    stage - the quickest way to see a drift, an outlier stretch or a stage effect."""
    _, colors, _, _ = stage_style(custom_stages)
    epoch_idx = meta['epoch_idx'].values
    stages = meta['stage'].astype(str).values
    n_ch = len(ch_names)
    fig, axes = plt.subplots(2, max(n_ch, 1), figsize=(_fig_width(n_ch), 6.0),
                             sharex=True, squeeze=False)
    for ci, ch in enumerate(ch_names):
        for row, key, lab in [(0, 'exponent', 'Exponent'), (1, 'offset', 'Offset')]:
            ax = axes[row][ci]
            vals = res[key][:, ci]
            for stage in dict.fromkeys(stages):
                m = stages == stage
                ax.scatter(epoch_idx[m], vals[m], s=5, alpha=0.6,
                           color=colors.get(stage, '#888888'),
                           label=stage if (ci == 0 and row == 0) else None)
            ax.grid(alpha=0.2)
            if ci == 0:
                ax.set_ylabel(lab)
            if row == 0:
                ax.set_title(ch)
            else:
                ax.set_xlabel('Epoch index')
        axes[0][ci].set_ylim(*np.nanpercentile(res['exponent'], [0.5, 99.5])
                             if np.isfinite(res['exponent']).any() else (0, 1))
    if n_ch:
        axes[0][0].legend(fontsize=8, markerscale=2)
    fig.suptitle('Aperiodic parameters across the night (one point per epoch)')
    fig.tight_layout()
    return fig


# --- 6. band power per stage -------------------------------------------------
def plot_band_power(band_vals, ch_names, masks, custom_stages, measure, measure_label):
    """Boxplot of the per-epoch band power, grouped by band and coloured by stage."""
    _, colors, _, _ = stage_style(custom_stages)
    bands = list(band_vals)
    stages = [s for s in masks if masks[s].sum() > 0]
    fig, axes = _axes_row(len(ch_names), height=4.6, sharey=True)
    width = 0.8 / max(len(stages), 1)
    for ci, ch in enumerate(ch_names):
        ax = axes[ci]
        for si, stage in enumerate(stages):
            data, pos = [], []
            for bi, band in enumerate(bands):
                v = band_vals[band].get(measure)
                if v is None:
                    continue
                vals = v[masks[stage], ci]
                vals = vals[np.isfinite(vals)]
                if len(vals) == 0:
                    continue
                data.append(vals)
                pos.append(bi + (si - (len(stages) - 1) / 2) * width)
            if not data:
                continue
            bp = ax.boxplot(data, positions=pos, widths=width * 0.85, patch_artist=True,
                            showfliers=False, medianprops=dict(color='black', lw=0.8))
            for patch in bp['boxes']:
                patch.set_facecolor(colors.get(stage, '#888888'))
                patch.set_alpha(0.75)
                patch.set_linewidth(0.5)
        ax.set_xticks(range(len(bands)))
        ax.set_xticklabels(bands, rotation=30, ha='right')
        ax.set_title(ch)
        ax.grid(alpha=0.2, axis='y')
        if ci == 0:
            ax.set_ylabel(measure_label)
    handles = [matplotlib.patches.Patch(facecolor=colors.get(s, '#888888'), alpha=0.75, label=s)
               for s in stages]
    if len(axes):
        axes[0].legend(handles=handles, fontsize=8, ncol=2)
    fig.suptitle(f'Band power per stage - {measure_label}')
    fig.tight_layout()
    return fig


# --- HTML tables -------------------------------------------------------------
def params_html(cfg):
    """Two-column 'parameter -> value' table of everything the run used."""
    rows = ''.join(f'<tr><td style="padding:2px 12px 2px 0"><b>{k}</b></td>'
                   f'<td>{v}</td></tr>' for k, v in cfg.items())
    return f'<table style="font-size:90%">{rows}</table>'


def coverage_html(agg, ch_names):
    """n_epochs (and n_fit_ok) per stage x channel - shows where `min epochs` produced NaN."""
    if agg is None or len(agg) == 0:
        return '<p><i>(no data)</i></p>'
    piv = agg.drop_duplicates(['stage', 'channel'])[['stage', 'channel', 'n_epochs', 'n_fit_ok']]
    wide = piv.pivot(index='stage', columns='channel', values='n_epochs').reindex(columns=ch_names)
    return wide.reset_index().to_html(index=False, border=0, na_rep='')


def guidance_html(cfg, data_dir, reports_dir, has_subject_info=False):
    """"Where do I look?" - the closing section of the database report.

    The run writes seven kinds of table; without a map it is easy to open the wrong one. This
    names the two files that answer most questions, then lists the rest with what each is for."""
    space = 'from_log' if 'log' in cfg['spaces'] else 'from_lin'
    band_col = (f'power_db_{space}_mean' if 'power_db' in cfg['measures']
                else (f'{cfg["measures"][0]}_mean' if cfg['measures'] else 'power_db_mean'))
    rows = [
        ('<b>Band power per participant</b>, ready for statistics',
         '<code>global_spectral_stage.tsv</code>', 'one row per participant x stage x channel x band; '
         f'the main column is <code>{band_col}</code>'),
        ('<b>Aperiodic exponent and offset</b> per participant',
         '<code>global_aperiodic_stage.tsv</code>',
         '<code>exponent_mean</code> and <code>offset_mean</code>, one row per participant x stage '
         'x channel'),
        ('Everything above in one file', '<code>spectral_features_database.xlsx</code>',
         'the same tables as Excel sheets, plus <code>coverage</code> and <code>parameters</code>'),
        ('Redraw a spectrum without recomputing', '<code>global_psd_stage.tsv</code>',
         'the mean PSD per stage x channel x frequency bin'),
        ('Look inside one participant', '<code>{file_id}_spectral_report.html</code>',
         'the per-participant figures and the parameters actually used'),
        ('Work at the epoch level (time course, outliers)',
         '<code>{file_id}_bandpower_epoch.tsv</code> / <code>{file_id}_aperiodic_epoch.tsv</code>',
         'one row per epoch x channel (x band); this is where a night-long time course lives'),
        ('Individual oscillatory peaks', '<code>{file_id}_periodic_peaks.tsv</code>',
         'centre frequency, power and bandwidth of every peak specparam fitted'),
    ]
    redraw = (
        '<p><small>Every figure of these reports can be redrawn from the tables, so you are never '
        'forced to re-run the tool to change a grouping:</small></p>'
        '<table style="font-size:90%">'
        '<tr><td style="padding:2px 14px 2px 0">mean PSD per stage (with its +/- SEM band)</td>'
        '<td><code>global_psd_stage.tsv</code>: <code>sem = psd_*_sd / sqrt(n_epochs)</code></td></tr>'
        '<tr><td style="padding:2px 14px 2px 0">PSD with the aperiodic component removed</td>'
        '<td><code>global_psd_stage.tsv</code>: <code>psd_ap_removed_db</code></td></tr>'
        '<tr><td style="padding:2px 14px 2px 0">log-log PSD + the aperiodic fit</td>'
        '<td>the same table, plus <code>offset_mean</code> / <code>exponent_mean</code> from '
        '<code>global_aperiodic_stage.tsv</code>: '
        '<code>ap_dB(f) = 10*(offset - exponent*log10(f))</code></td></tr>'
        '<tr><td style="padding:2px 14px 2px 0">exponent / offset across the night</td>'
        '<td><code>{file_id}_aperiodic_epoch.tsv</code></td></tr>'
        '<tr><td style="padding:2px 14px 2px 0">band-power boxplots, clean-epoch hypnogram</td>'
        '<td><code>{file_id}_bandpower_epoch.tsv</code> (it carries <code>epoch_idx</code> and '
        '<code>stage</code>)</td></tr></table>')
    joined = ('The participant-info table you selected is <b>already joined</b> into the global '
              'tables, so its columns can be used directly as grouping keys.'
              if has_subject_info else
              'No participant-info table was selected; to group by population variables, merge your '
              'own table on <code>file_id</code> (or re-run with the <i>Participant info</i> '
              'chooser filled in).')
    regroup = (
        f'<p><small>{joined}</small></p>'
        '<pre style="font-size:85%;background:#f6f6f6;padding:8px">'
        "import pandas as pd\n"
        "df = pd.read_csv('global_spectral_stage.tsv', sep='\\t', dtype={'file_id': str})\n"
        "# df = df.merge(pd.read_csv('my_population.csv'), on='file_id')   # if not joined already\n"
        "sel = df[(df.stage == 'N2') &amp; (df.band == 'sigma') &amp; (df.n_epochs >= 20)]\n"
        "sel.groupby(['group', 'channel'])['power_db_from_log_mean'].agg(['mean', 'sem', 'count'])"
        '</pre>')
    body = ''.join(f'<tr><td style="padding:3px 14px 3px 0;vertical-align:top">{q}</td>'
                   f'<td style="padding:3px 14px 3px 0;vertical-align:top">{f}</td>'
                   f'<td style="padding:3px 0;vertical-align:top"><small>{d}</small></td></tr>'
                   for q, f, d in rows)
    caveats = []
    if 'log' in cfg['spaces'] and 'lin' in cfg['spaces']:
        caveats.append('Both averaging spaces were exported: <code>_from_log</code> (mean of the '
                       'per-epoch dB) and <code>_from_lin</code> (dB of the linear mean). Pick one '
                       'and stay with it - mixing them across analyses is not comparable.')
    caveats.append('A cell computed on fewer than '
                   f'{cfg["min_epochs"]} epochs is <b>NaN</b>, but its <code>n_epochs</code> is '
                   'kept: always check that column before reading a value.')
    caveats.append('Rows with <code>n_epochs = 0</code> are channels a participant does not have '
                   '(dropped by tool 8/8bis); they are there so the table stays rectangular.')
    if any(m in cfg['measures'] for m in ('power_ap_removed_db', 'power_ratio_over_aperiodic',
                                          'power_rel_periodic')):
        caveats.append('The 1/f-corrected measures rest on the specparam fit, which is extrapolated '
                       'below its lower bound - read <b>delta</b> with care and check the '
                       '<code>n_fit_ok</code> column.')
    return (f'<p>Data tables: <code>{data_dir}</code><br>'
            f'Reports and database tables: <code>{reports_dir}</code></p>'
            f'<table style="font-size:95%">{body}</table>'
            '<p style="margin-top:14px"><b>Before you read a number</b></p><ul>'
            + ''.join(f'<li><small>{c}</small></li>' for c in caveats) + '</ul>'
            + '<p style="margin-top:14px"><b>Redrawing a figure yourself</b></p>' + redraw
            + '<p style="margin-top:14px"><b>Grouping participants differently</b></p>' + regroup)


def save_report_html(path, title, items):
    """Write an mne.Report from an ORDERED list of ('fig'|'html', title, payload) items.

    Ordered rather than "all figures then all text" so a figure can be placed after the section it
    belongs with. `None` payloads are skipped, which keeps the call sites free of conditionals."""
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
    """Close every figure of an item list (matplotlib keeps them open otherwise)."""
    for kind, _, payload in items:
        if kind == 'fig' and payload is not None:
            plt.close(payload)


# --- database-level figures --------------------------------------------------
def pick_psd_column(psd_df, prefer=('psd_db_from_log', 'psd_db_from_lin')):
    """First available dB column of the global PSD table (which one exists depends on the
    averaging space the user selected)."""
    if psd_df is None or len(psd_df) == 0:
        return None
    for col in prefer:
        if col in psd_df.columns and psd_df[col].notna().any():
            return col
    return None


def plot_group_psd(psd_df, custom_stages, value_col=None,
                   ylabel='PSD (dB/Hz)', title='Group mean PSD per stage'):
    """Mean across participants (+/- SEM between subjects) of the per-stage mean PSD."""
    if psd_df is None or len(psd_df) == 0:
        return None
    value_col = value_col or pick_psd_column(psd_df)
    if value_col is None or value_col not in psd_df.columns:
        return None
    _, colors, _, _ = stage_style(custom_stages)
    ch_names = sorted(psd_df['channel'].unique())
    fig, axes = _axes_row(len(ch_names), height=4.4)
    for ci, ch in enumerate(ch_names):
        ax = axes[ci]
        sub_ch = psd_df[psd_df['channel'] == ch]
        for stage in [s for s in AASM_STAGES + list(custom_stages)
                      if s in set(sub_ch['stage'])]:
            sub = sub_ch[sub_ch['stage'] == stage]
            piv = sub.pivot_table(index='freq_hz', columns='file_id', values=value_col)
            if piv.empty:
                continue
            freqs = piv.index.values
            mean = np.nanmean(piv.values, axis=1)
            n = np.sum(np.isfinite(piv.values), axis=1)
            sem = np.nanstd(piv.values, axis=1, ddof=1) / np.sqrt(np.maximum(n, 1))
            col = colors.get(stage, '#888888')
            ax.plot(freqs, mean, color=col, label=f'{stage} (n={int(np.nanmax(n))})')
            ax.fill_between(freqs, mean - sem, mean + sem, color=col, alpha=0.25)
        ax.set_title(ch)
        ax.set_xlabel('Frequency (Hz)')
        ax.grid(alpha=0.2)
        if ci == 0:
            ax.set_ylabel(ylabel)
            ax.legend(fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_group_band(stage_df, column, label, custom_stages):
    """Distribution across participants of one aggregated band measure, per stage and band."""
    if stage_df is None or len(stage_df) == 0 or column not in stage_df.columns:
        return None
    _, colors, _, _ = stage_style(custom_stages)
    ch_names = sorted(stage_df['channel'].unique())
    bands = list(dict.fromkeys(stage_df['band']))
    stages = [s for s in AASM_STAGES + list(custom_stages) if s in set(stage_df['stage'])]
    fig, axes = _axes_row(len(ch_names), height=4.6)
    width = 0.8 / max(len(stages), 1)
    for ci, ch in enumerate(ch_names):
        ax = axes[ci]
        sub_ch = stage_df[stage_df['channel'] == ch]
        for si, stage in enumerate(stages):
            data, pos = [], []
            for bi, band in enumerate(bands):
                v = sub_ch[(sub_ch['stage'] == stage) & (sub_ch['band'] == band)][column].values
                v = v[np.isfinite(v)]
                if len(v) == 0:
                    continue
                data.append(v)
                pos.append(bi + (si - (len(stages) - 1) / 2) * width)
            if not data:
                continue
            bp = ax.boxplot(data, positions=pos, widths=width * 0.85, patch_artist=True,
                            showfliers=False, medianprops=dict(color='black', lw=0.8))
            for patch in bp['boxes']:
                patch.set_facecolor(colors.get(stage, '#888888'))
                patch.set_alpha(0.75)
                patch.set_linewidth(0.5)
        ax.set_xticks(range(len(bands)))
        ax.set_xticklabels(bands, rotation=30, ha='right')
        ax.set_title(ch)
        ax.grid(alpha=0.2, axis='y')
        if ci == 0:
            ax.set_ylabel(label)
    handles = [matplotlib.patches.Patch(facecolor=colors.get(s, '#888888'), alpha=0.75, label=s)
               for s in stages]
    if len(axes):
        axes[0].legend(handles=handles, fontsize=8, ncol=2)
    fig.suptitle(f'{label} across participants (one point per participant)')
    fig.tight_layout()
    return fig


def plot_group_aperiodic(ap_df, custom_stages):
    """Exponent and offset per stage across participants (one point per participant)."""
    if ap_df is None or len(ap_df) == 0:
        return None
    _, colors, _, _ = stage_style(custom_stages)
    ch_names = sorted(ap_df['channel'].unique())
    stages = [s for s in AASM_STAGES + list(custom_stages) if s in set(ap_df['stage'])]
    fig, axes = plt.subplots(2, max(len(ch_names), 1), figsize=(_fig_width(len(ch_names)), 6.4),
                             squeeze=False)
    for ci, ch in enumerate(ch_names):
        sub_ch = ap_df[ap_df['channel'] == ch]
        for row, col_name, lab in [(0, 'exponent_mean', 'Exponent'), (1, 'offset_mean', 'Offset')]:
            ax = axes[row][ci]
            data, keep = [], []
            for stage in stages:
                v = sub_ch[sub_ch['stage'] == stage][col_name].values
                v = v[np.isfinite(v)]
                if len(v):
                    data.append(v)
                    keep.append(stage)
            if data:
                bp = ax.boxplot(data, patch_artist=True, showfliers=False,
                                medianprops=dict(color='black', lw=0.8))
                for patch, stage in zip(bp['boxes'], keep):
                    patch.set_facecolor(colors.get(stage, '#888888'))
                    patch.set_alpha(0.75)
                ax.set_xticks(range(1, len(keep) + 1))
                ax.set_xticklabels(keep)
            ax.grid(alpha=0.2, axis='y')
            if ci == 0:
                ax.set_ylabel(lab)
            if row == 0:
                ax.set_title(ch)
    fig.suptitle('Aperiodic parameters per stage across participants')
    fig.tight_layout()
    return fig


def plot_group_counts(stage_df):
    """Clean epochs per participant, split by stage - the sample size behind every value."""
    if stage_df is None or len(stage_df) == 0:
        return None
    sub = stage_df[stage_df['stage'].isin(AASM_STAGES)]
    if len(sub) == 0:
        return None
    counts = (sub.drop_duplicates(['file_id', 'stage', 'channel'])
              .groupby(['file_id', 'stage'])['n_epochs'].max().unstack(fill_value=0))
    counts = counts.reindex(columns=[s for s in AASM_STAGES if s in counts.columns])
    fig, ax = plt.subplots(figsize=(max(6.0, 0.45 * len(counts) + 3), 4.0))
    bottom = np.zeros(len(counts))
    for stage in counts.columns:
        vals = counts[stage].values.astype(float)
        ax.bar(range(len(counts)), vals, bottom=bottom,
               color=BASE_STAGE_COLORS.get(stage, '#888888'), label=stage)
        bottom += vals
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=90, fontsize=7)
    ax.set_ylabel('clean epochs')
    ax.set_title('Clean epochs per participant and stage')
    ax.legend(fontsize=8, ncol=5)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------------------------------------
# Lines of tool 9's run cell, as functions (tool 0's step 8 calls them too)
# ---------------------------------------------------------------------------------------------------------
def stage_order(stages_present, custom_stages):
    """Tool 9's stage checkboxes: the AASM stages, then the custom ones, then any other label present."""
    order = [s for s in AASM_STAGES + list(custom_stages) if s in stages_present]
    order += [s for s in stages_present if s not in order]
    return order


def prepare_epochs(epochs):
    """Tool 9's load step on epochs already read (tool 7's .fif): the metadata (with an epoch_idx), the
    channel names and the bad channels, then the bad channels cleared so MNE's PSD keeps them.
    Returns (meta, ch_names, bads)."""
    meta = epochs.metadata.reset_index(drop=True).copy()
    if 'stage' not in meta.columns:
        raise ValueError('the .fif metadata has no `stage` column')
    if 'epoch_idx' not in meta.columns:
        meta['epoch_idx'] = np.arange(len(meta))
    ch_names = list(epochs.ch_names)
    bads = list(epochs.info['bads'])
    epochs.info['bads'] = []     # MNE's compute_psd skips bad channels: their features are wanted too
    return meta, ch_names, bads


def fit_aperiodic(freqs, psds, cfg, progress=None, warn=print):
    """Tool 9's aperiodic step: the fit copy smoothed (median then LOWESS), the specparam fit of every
    (epoch, channel), the aperiodic model on the PSD axis, and a warning when the knee is degenerate.
    Returns (res, peak_rows, ap_log10)."""
    psd_fit = smooth_for_fit(psds, freqs, cfg['smooth'])
    res, peak_rows = fit_spectra(freqs, psd_fit, cfg['fit'], n_jobs=cfg['n_jobs'],
                                 progress=progress)
    del psd_fit
    ap_log10 = aperiodic_psd_log10(freqs, res, cfg['fit']['aperiodic_mode'])
    if cfg['fit']['aperiodic_mode'] == 'knee':
        med_knee = np.nanmedian(np.abs(res['knee']))
        if np.isfinite(med_knee) and med_knee > KNEE_DEGENERATE:
            warn(f'  &#9888; the knee parameter is degenerate (median |knee| = {med_knee:.3g}): '
                 f'this spectrum has no real bend: use the <code>fixed</code> mode')
    return res, peak_rows, ap_log10


def usable_bands(freqs, bands, warn=print):
    """The bands with at least two frequency bins inside the PSD range (a band partly outside is kept but
    warned about). Raises ValueError when none is left."""
    # A band needs at least two frequency bins to mean anything; one that is only partly inside the
    # PSD range is kept but flagged, since its value is not comparable to a fully covered band.
    df_hz = float(np.median(np.diff(freqs)))
    usable = []
    for n, lo, hi in bands:
        n_bins = int(((freqs >= lo) & (freqs < hi)).sum())
        if n_bins < 2:
            warn(f'  &#9888; band "{n}" ({lo:g}-{hi:g} Hz) has {n_bins} frequency bin(s) inside the '
                 f'PSD range - skipped')
            continue
        span = min(hi, float(freqs[-1]) + df_hz) - max(lo, float(freqs[0]))
        if span < 0.99 * (hi - lo):
            warn(f'  &#9888; band "{n}" ({lo:g}-{hi:g} Hz) is only {100 * span / (hi - lo):.0f} % '
                 f'inside the PSD range - its power is not comparable to a fully covered band')
        usable.append((n, lo, hi))
    if not usable:
        raise ValueError('no frequency band fits inside the PSD range')
    return usable


def participant_tables(fid, source, cfg, meta, ch_names, freqs, psds, res, peak_rows, ap_log10, band_vals,
                       in_scope, rejected, dropped, warn=print):
    """Tool 9's tables of one participant (the first part of its finish_participant): the decision columns,
    the stage aggregates (kept epochs only, every channel), the per-epoch tables and the night-thirds table.
    Returns a dict of DataFrames (None when not computed), plus kept / stages_sel / not_analysed."""
    stages = meta['stage'].astype(str).values
    kept = in_scope & ~rejected
    # night thirds over the WHOLE night (the sleep period does not depend on the rejection decision)
    thirds = assign_thirds(meta['epoch_idx'].values, stages) if cfg['thirds'] else None
    stages_sel = [s for s in cfg['stages'] if s in set(stages[kept])]
    # A stage asked for here but with no kept epoch has no row in the stage tables. Say why, so it is
    # not mistaken for a bug: out of scope in the tool 8/8bis decision (stage unticked there), every
    # epoch rejected, or simply absent from the hypnogram.
    not_analysed = []
    for st in cfg['stages']:
        if st in stages_sel:
            continue
        n_st = int((stages == st).sum())
        n_st_oos = int(((stages == st) & ~in_scope).sum())
        n_st_rej = int(((stages == st) & in_scope & rejected).sum())
        if n_st == 0:
            why = 'no epoch of this stage in the recording'
        elif n_st_oos == n_st:
            why = f'all {n_st} epochs out of scope (stage not selected in the tool 8/8bis decision)'
        elif n_st_rej == n_st:
            why = f'all {n_st} epochs rejected'
        else:
            why = f'{n_st} epochs: {n_st_rej} rejected, {n_st_oos} out of scope'
        not_analysed.append(f'{st} ({why})')
    not_analysed = '; '.join(not_analysed)
    masks_el = {st: m & kept for st, m in stage_masks(stages, stages_sel).items()}
    fit_ok = res['fit_ok'] if res is not None else None
    n_rej = {st: int(((stages == st) & in_scope & rejected).sum()) for st in set(stages)}
    n_oos = {st: int(((stages == st) & ~in_scope).sum()) for st in set(stages)}
    pos_e = {int(e): i for i, e in enumerate(meta['epoch_idx'])}

    def decision_cols_stage(df):
        if df is None or len(df) == 0:
            return df
        df['channel_dropped'] = df['channel'].isin(dropped)
        df['n_epochs_rejected'] = df['stage'].map(n_rej).fillna(0).astype(int)
        df['n_epochs_out_of_scope'] = df['stage'].map(n_oos).fillna(0).astype(int)
        return df

    def decision_cols_epoch(df):
        if df is None or len(df) == 0:
            return df
        idx = df['epoch_idx'].map(pos_e).values
        df['in_scope'] = in_scope[idx]
        df['rejected'] = rejected[idx]
        df['channel_dropped'] = df['channel'].isin(dropped)
        return df

    agg = decision_cols_stage(aggregate_bandpower(fid, source, band_vals, ch_names, masks_el, cfg['measures'],
                                                  cfg['min_epochs'], spaces=cfg['spaces'], fit_ok=fit_ok,
                                                  exclude_bad_fit=cfg['exclude_bad_fit']))
    psd_stage = decision_cols_stage(aggregate_psd(fid, source, freqs, psds, ap_log10, ch_names, masks_el,
                                                  cfg['min_epochs'], spaces=cfg['spaces']))
    ap_stage = (decision_cols_stage(aggregate_aperiodic(fid, source, res, ch_names, masks_el, cfg['min_epochs'],
                                                        exclude_bad_fit=cfg['exclude_bad_fit']))
                if res is not None else None)
    if not stages_sel:
        warn('  &#9888; no kept epoch in the selected stages: the stage tables are empty')
    bandpower_epoch = decision_cols_epoch(build_bandpower_epoch(fid, source, meta, ch_names, band_vals,
                                                                cfg['measures'], thirds))
    aperiodic_epoch = peaks = None
    if res is not None:
        aperiodic_epoch = decision_cols_epoch(build_aperiodic_epoch(fid, source, meta, ch_names, res, thirds))
        peaks = decision_cols_epoch(build_peaks_table(fid, source, meta, ch_names, peak_rows))
    third_agg = None
    if cfg['thirds']:
        masks_third = {f'{st}|{th}': (stages == st) & (thirds == th) & kept
                       for st in stages_sel for th in ['T1', 'T2', 'T3']}
        third_agg = aggregate_bandpower(fid, source, band_vals, ch_names, masks_third,
                                        cfg['measures'], cfg['min_epochs'], spaces=cfg['spaces'],
                                        fit_ok=fit_ok, exclude_bad_fit=cfg['exclude_bad_fit'])
        if len(third_agg):
            third_agg[['stage', 'third']] = third_agg['stage'].str.split('|', expand=True)
            # the counts of THIS stage x third cell (not of the whole stage), so a third resting on a small
            # share of its epochs (e.g. an artefact-ridden end of night) shows: pct_rejected =
            # rejected / (kept + rejected) of the in-scope epochs of the cell
            in_cell = {(st, th): (stages == st) & (thirds == th) for st in stages_sel for th in ['T1', 'T2', 'T3']}
            keys = list(zip(third_agg['stage'], third_agg['third']))
            n_kept_c = [int((in_cell[k] & kept).sum()) for k in keys]
            n_rej_c = [int((in_cell[k] & in_scope & rejected).sum()) for k in keys]
            third_agg['channel_dropped'] = third_agg['channel'].isin(dropped)
            third_agg['n_epochs_rejected'] = n_rej_c
            third_agg['n_epochs_out_of_scope'] = [int((in_cell[k] & ~in_scope).sum()) for k in keys]
            third_agg['pct_rejected'] = [round(100.0 * r / (k + r), 2) if (k + r) else float('nan')
                                         for k, r in zip(n_kept_c, n_rej_c)]
    return {'spectral_stage': agg, 'psd_stage': psd_stage, 'aperiodic_stage': ap_stage,
            'bandpower_epoch': bandpower_epoch, 'aperiodic_epoch': aperiodic_epoch, 'periodic_peaks': peaks,
            'spectral_stage_third': third_agg, 'kept': kept, 'stages_sel': stages_sel,
            'not_analysed': not_analysed}


def participant_report_items(p, cfg, source, custom, tables, meta, ch_names, freqs, psds, res, peak_rows,
                             ap_log10, band_vals, usable, in_scope, rejected, dropped, warn=print):
    """Tool 9's per-participant report (the second part of its finish_participant): the figures, the
    parameters and the coverage, on the kept epochs of the kept channels. Returns the report items.
    `p`: tool 9's participant record (file_id, fif, sfreq, epoch_length_s, sidecar_info, epoch_decision,
    excluded, exclude_reason); `tables`: participant_tables' result."""
    fid = p['file_id']
    kept, stages_sel, not_analysed = tables['kept'], tables['stages_sel'], tables['not_analysed']
    agg = tables['spectral_stage']
    good = [ci for ci, ch in enumerate(ch_names) if ch not in dropped]
    kidx = np.flatnonzero(kept)
    ch_k = [ch_names[ci] for ci in good]
    meta_k = meta.iloc[kidx].reset_index(drop=True)
    psds_k = psds[kidx][:, good]
    res_k = {k: np.asarray(v)[kidx][:, good] for k, v in res.items()} if res is not None else None
    ap_k = ap_log10[kidx][:, good] if ap_log10 is not None else None
    band_k = {b: {m: v[kidx][:, good] for m, v in vals.items()} for b, vals in band_vals.items()}
    masks_k = stage_masks(meta_k['stage'].astype(str).values, stages_sel)
    psd_lin = np.where(psds_k > 0, psds_k, np.nan)
    psd_db = 10.0 * np.log10(psd_lin)
    use_log, use_lin = 'log' in cfg['spaces'], 'lin' in cfg['spaces']
    items = []
    if p.get('excluded'):
        items.append(('html', 'Excluded participant', PSL.excluded_banner_html(p['exclude_reason'])))
    items += [('fig', 'Clean epochs', plot_clean_hypnogram(meta_k, custom, fid, n_total=len(meta))),
              ('fig', 'Mean PSD per stage',
               plot_psd_per_stage(freqs, masks_k, ch_k, custom, cfg['min_epochs'],
                                  'Mean PSD per stage', 'PSD (dB/Hz)',
                                  values_db=psd_db if use_log else None,
                                  values_lin=psd_lin if use_lin else None))]
    if ap_k is not None:
        # In linear space the "1/f removed" spectrum is the PSD / aperiodic-fit ratio, whose dB is
        # exactly the log-space difference averaged the other way round.
        items.append(('fig', 'PSD after removing the aperiodic component',
                      plot_psd_per_stage(freqs, masks_k, ch_k, custom, cfg['min_epochs'],
                                         'Mean PSD per stage after removing the aperiodic component',
                                         'Power above the 1/f fit (dB)',
                                         values_db=(psd_db - 10.0 * ap_k) if use_log else None,
                                         values_lin=(psd_lin / 10.0 ** ap_k) if use_lin else None)))
        items.append(('fig', 'log-log PSD with the aperiodic fit',
                      plot_loglog_with_fit(freqs, psds_k, ap_k, masks_k, ch_k, custom,
                                           cfg['min_epochs'],
                                           (cfg['fit']['fit_fmin'], cfg['fit']['fit_fmax']), res_k)))
        items.append(('fig', 'Aperiodic parameters across the night',
                      plot_aperiodic_over_night(meta_k, res_k, ch_k, custom)))
    band_measure = ('power_db' if 'power_db' in cfg['measures'] else cfg['measures'][0])
    band_label = dict((k, lbl) for k, lbl, _, _, _ in MEASURES)[band_measure]
    items.append(('fig', f'Band power per stage ({band_label})',
                  plot_band_power(band_k, ch_k, masks_k, custom, band_measure, band_label)))

    side = p['sidecar_info']
    params = {
        'Source': f"{source} ({p['fif'].name}"
                  + (f", decisions {p['epoch_decision'].name})" if p.get('epoch_decision') else ')'),
        'Epochs analysed': (f'{len(meta_k)} kept of {len(meta)} in the recording '
                            f'({int((in_scope & rejected).sum())} rejected, {int((~in_scope).sum())} out of '
                            f'scope; every epoch is in the per-epoch tables)'),
        'Channels': ', '.join(ch_k) + (f' (dropped, computed but not shown: {", ".join(dropped)})'
                                       if dropped else ''),
        'Sampling rate': f"{p['sfreq']:g} Hz",
        'Epoch length': f"{p['epoch_length_s']:g} s",
        'PSD': (f"{cfg['psd']['method']}, {cfg['psd']['fmin']:g}-{cfg['psd']['fmax']:g} Hz, "
                + (f"{cfg['psd']['win_s']:g} s {cfg['psd']['window']} window, "
                   f"{cfg['psd']['overlap_pct']:g} % overlap"
                   if cfg['psd']['method'] == 'welch'
                   else f"bandwidth {cfg['psd']['bandwidth']:g} Hz")
                + f", resolution {freqs[1] - freqs[0]:.3g} Hz"),
        'Bands': ', '.join(f'{n} {lo:g}-{hi:g}' for n, lo, hi in usable),
        'Measures': ', '.join(cfg['measures']),
        'Epoch averaging': ' + '.join({'log': 'log space (mean of the per-epoch dB)',
                                       'lin': 'linear space, then converted to dB'}[s]
                                      for s in cfg['spaces']),
        'Min epochs per cell': cfg['min_epochs'],
    }
    if not_analysed:
        # placed right after 'Epochs analysed' (absent when every requested stage was analysed)
        with_stages = {}
        for k, v in params.items():
            with_stages[k] = v
            if k == 'Epochs analysed':
                with_stages['Stages requested but not analysed'] = not_analysed
        params = with_stages
    if cfg['do_fit']:
        # Tool 7 already rejects epochs on a 1/f criterion, so "why do some spectra still fail?"
        # is the first question this line raises. Answer it in place, and name the settings that
        # differ from the tool-7 run when its sidecar is available.
        sm = cfg['smooth']
        applied = ([f"median {sm['median_span_hz']:g} Hz"] if sm['enabled'] and sm['median'] else []) \
            + ([f"LOWESS {sm['lowess_span_hz']:g} Hz"] if sm['enabled'] and sm['lowess'] else [])
        smooth_desc = ((' then '.join(applied) if applied else 'none')
                       + ' - the band powers use the raw PSD')
        diffs = fit_settings_diff(cfg, side)
        fit_note = ('Some failures are expected even though tool 7 already flagged on a 1/f '
                    'criterion: tool 7 flags each (epoch, channel) pair while tools 8/8bis decide '
                    'per epoch (an epoch flagged on one channel out of four can legitimately be '
                    'kept), and the fit is recomputed here with this notebook&#39;s own PSD and '
                    'specparam settings.')
        if diffs:
            fit_note += ' Settings that differ from the tool-7 run: ' + '; '.join(diffs) + '.'
        elif side['found']:
            fit_note += ' The 1/f settings match the tool-7 run.'
        params.update({
            '1/f fit': (f"specparam, {cfg['fit']['fit_fmin']:g}-{cfg['fit']['fit_fmax']:g} Hz, "
                        f"mode {cfg['fit']['aperiodic_mode']}, peak width "
                        f"{cfg['fit']['peak_width_min']:g}-{cfg['fit']['peak_width_max']:g} Hz, "
                        f"max {cfg['fit']['max_n_peaks']} peaks, min height "
                        f"{cfg['fit']['min_peak_height']:g}, threshold {cfg['fit']['peak_threshold']:g} SD"),
            'Fit quality bounds': (f"R2 >= {cfg['fit']['r2_min']:g} and MAE <= {cfg['fit']['mae_max']:g}"
                                   f" - {np.nanmean(res_k['fit_ok']) * 100:.0f} % of the kept spectra pass"
                                   + (', poor fits excluded from the aggregates'
                                      if cfg['exclude_bad_fit'] else ', all fits kept')
                                   + '<br><small style="color:#555">' + fit_note + '</small>'),
            'Smoothing (fit copy only)': smooth_desc,
        })
        ok = float(np.nanmean(res_k['fit_ok'])) * 100.0 if res_k['fit_ok'].size else float('nan')
        warn(f'  1/f fit: {ok:.0f} % of the kept spectra within the quality bounds, '
            f'{len(peak_rows)} peaks detected (all epochs)')
    params['Tool-7 provenance'] = (
        f"filter {side['filter']}, resample {side['resample']}, notch {side['notch']}, "
        f"PSD smoothing {side['psd_smoothing']}" if side['found']
        else '<span style="color:#ef6c00">no sidecar found - the preprocessing settings are unknown</span>')

    # The numbers themselves live in the TSV tables - the report only carries what helps read
    # them: the settings used and how many epochs stand behind each cell.
    agg_k = agg[~agg['channel_dropped']] if len(agg) else agg
    n_nan = int(agg_k.get('n_epochs', pd.Series(dtype=int)).lt(cfg['min_epochs']).sum()) if len(agg_k) else 0
    items += [
        ('html', 'Parameters used', params_html(params)),
        ('html', 'Coverage - clean epochs per stage and channel',
         coverage_html(agg_k, ch_k)
         + f'<p><small>{n_nan} (stage x channel x band) cells fall below the '
           f'{cfg["min_epochs"]}-epoch minimum and are written as NaN. The values themselves are in '
           f'<code>{fid}_spectral_stage.tsv</code>.</small></p>'),
    ]
    return items
