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

import mne
from scipy.signal import spectrogram as sp_spectrogram

try:
    from specparam import SpectralModel
    HAS_SPECPARAM = True
except Exception:
    HAS_SPECPARAM = False

mne.set_log_level('ERROR')

# ---- Rejection-method registry (single source of truth, kept in sync with tool 6) ----
METHOD_ORDER  = ['amplitude', 'flat', 'gradient', '1f_error', '1f_r2', 'event']
METHOD_CODE   = {m: i + 1 for i, m in enumerate(METHOD_ORDER)}   # 1..6  (0 = none)
MULTIPLE_CODE = len(METHOD_ORDER) + 1                            # 7 = multiple
METHOD_LABEL  = {'amplitude': 'Amplitude', 'flat': 'Flat', 'gradient': 'Gradient',
                 '1f_error': '1/f error', '1f_r2': '1/f R²', 'event': 'Event'}
# Heatmap / per-method colours (index = method code; identical to tool 6's plot_rejection_heatmap).
HEATMAP_COLORS = ['#1c0a3b', '#c0392b', '#2980b9', '#e67e22', '#f1c40f', '#27ae60', '#e84393', '#7b0000']
HEATMAP_LABELS = ['none', 'amplitude', 'flat', 'gradient', '1/f error', '1/f R2', 'event', 'multiple']
# Per-method overlay colour keyed by METHOD_ORDER label.
METHOD_COLOR   = {m: HEATMAP_COLORS[METHOD_CODE[m]] for m in METHOD_ORDER}

# Context-channel (EOG/EMG/ECG) trace colours for the per-epoch montage (role labels written by tool 6).
CTX_COLOR = {'EOG-L': '#2c7fb8', 'EOG-R': '#2c7fb8', 'EMG': '#27ae60', 'ECG': '#c0392b'}

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
            'epoch_length_s': 30}       # scoring epoch length (s); classic 30 s when absent
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


# ---------------------------------------------------------------------------
# Spectral analysis (Welch PSD + specparam 1/f fit) — same config as tool 6
# ---------------------------------------------------------------------------
def compute_psds(epochs, fmax=45.0):
    """Welch PSD per epoch/channel, identical config to 6_preprocessing_voila.
    `fmax` is the 1/f fit upper bound (default 2–45 Hz range); the PSD spans up to it.
    Returns (freqs, psds_uV2) with psds shape (n_ep, n_ch, n_freqs) in µV²/Hz."""
    sf = float(epochs.info['sfreq'])
    epoch_len_s = len(epochs.times) / sf                  # supports non-30 s epochs
    n_per_seg = int(min(4, epoch_len_s) * sf)              # 4 s Welch window, capped to the epoch length
    n_overlap = int(n_per_seg / 2)
    fmax_psd = min(fmax, sf / 2 - 0.5)
    psds_obj = epochs.compute_psd(method='welch', fmin=0.5, fmax=fmax_psd,
                                  n_fft=n_per_seg, n_overlap=n_overlap, n_per_seg=n_per_seg,
                                  verbose=False)
    return psds_obj.freqs, psds_obj.get_data() * 1e12


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


def compute_epoch_metrics(data_uV, freqs, psds_uV2, progress=None, fmin=2.0):
    """Per-epoch rejection-metric VALUES (for the report plots — the flags themselves come from
    the .fif metadata). Peak-to-peak and gradient are per (epoch, channel); the per-epoch value
    is the worst channel (max), matching 'any channel flags the epoch'. 1/f MAE/R² are fitted per
    (epoch, channel); the per-epoch value is the worst channel (max MAE, min R²).

    `progress(done, total)` is called periodically so a UI can show advancement.
    Returns dict of (n_ep,) arrays: ptp, gradient, mae, r2, and (n_ep, n_ch) arrays ptp_ch, grad_ch.
    """
    n_ep, n_ch, _ = data_uV.shape
    ptp_ch = np.ptp(data_uV, axis=-1)                              # (n_ep, n_ch)
    grad_ch = np.max(np.abs(np.diff(data_uV, axis=-1)), axis=-1)   # (n_ep, n_ch)
    mae = np.full(n_ep, np.nan)
    r2 = np.full(n_ep, np.nan)
    if HAS_SPECPARAM and psds_uV2 is not None:
        total = n_ep
        for ei in range(n_ep):
            worst_mae, best_r2 = np.nan, np.nan
            for ci in range(n_ch):
                m, r, _o, _e, _sm = fit_1f(freqs, psds_uV2[ei, ci], fmin=fmin)
                if not np.isnan(m):
                    worst_mae = m if np.isnan(worst_mae) else max(worst_mae, m)
                    best_r2 = r if np.isnan(best_r2) else min(best_r2, r)
            mae[ei], r2[ei] = worst_mae, best_r2
            if progress is not None and (ei % 25 == 0 or ei == total - 1):
                progress(ei + 1, total)
    return {'ptp': ptp_ch.max(axis=1), 'gradient': grad_ch.max(axis=1),
            'mae': mae, 'r2': r2, 'ptp_ch': ptp_ch, 'grad_ch': grad_ch}


def ordered_present_stages(stages, custom_stages):
    """AASM order then custom stages, keeping only stages actually present."""
    present = set(np.asarray(stages).tolist())
    order = AASM_STAGES + [s for s in custom_stages if s not in AASM_STAGES]
    ordered = [s for s in order if s in present]
    ordered += [s for s in present if s not in ordered]   # any unexpected label, last
    return ordered


# ---------------------------------------------------------------------------
# Section 2 — per-stage report figures (shared by notebook + batch)
# ---------------------------------------------------------------------------
def plot_psd_overlay(freqs, psds_uV2, stages, reject_flag, reject_method, stage, custom_stages=(),
                     title_prefix=''):
    """Per-epoch mean-across-channel PSD for one stage: clean epochs in light grey, rejected epochs
    coloured by their reject method, plus median + IQR band of the clean epochs. Returns a Figure."""
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    sel = (np.asarray(stages) == stage)
    if sel.sum() == 0:
        ax.text(0.5, 0.5, f'No {stage} epochs', ha='center', va='center')
        ax.set_axis_off()
        return fig
    psd_ep = np.asarray(psds_uV2)[sel].mean(axis=1)      # (n_stage_ep, n_freqs), µV²/Hz
    rej = np.asarray(reject_flag)[sel]
    meth = np.asarray(reject_method)[sel]
    # clean epochs (grey) + IQR band
    clean = psd_ep[~rej]
    for row in clean:
        ax.semilogy(freqs, row, color='0.8', lw=0.4, alpha=0.5, zorder=1)
    if clean.shape[0] >= 3:
        med = np.median(clean, axis=0)
        q1, q3 = np.percentile(clean, [25, 75], axis=0)
        ax.fill_between(freqs, q1, q3, color='0.5', alpha=0.25, zorder=2, label='clean IQR')
        ax.semilogy(freqs, med, color='0.25', lw=1.6, zorder=3, label='clean median')
    # rejected epochs coloured by method
    seen = set()
    for row, mth in zip(psd_ep[rej], meth[rej]):
        base = mth if mth in METHOD_COLOR else ('multiple' if mth == 'multiple' else 'amplitude')
        col = HEATMAP_COLORS[MULTIPLE_CODE] if mth == 'multiple' else METHOD_COLOR.get(mth, '#c0392b')
        lab = None
        if mth not in seen:
            seen.add(mth)
            lab = mth
        ax.semilogy(freqs, row, color=col, lw=0.7, alpha=0.7, zorder=4, label=lab)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('PSD (µV²/Hz)')
    n_rej = int(rej.sum())
    ax.set_title(f'{title_prefix}{stage} — PSD overlay  '
                 f'({sel.sum()} epochs, {n_rej} rejected)', fontsize=10)
    ax.legend(fontsize=7, ncol=2, loc='upper right')
    fig.tight_layout()
    return fig


def plot_metric_distributions(metrics, stages, reject_flag, thresholds, stage_order, title_prefix=''):
    """Box/strip of the four rejection metrics (p-p, gradient, 1/f MAE, 1/f R²), clean vs rejected,
    per stage, with the threshold line drawn. Returns a Figure."""
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
                v = vals[(stages == st) & sub]
                v = v[~np.isnan(v)]
                x = base_x + j
                positions.append(x)
                if v.size:
                    ax.scatter(np.full(v.size, x) + np.random.uniform(-0.15, 0.15, v.size),
                               v, s=4, color=colour, alpha=0.35, zorder=2)
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
                      f'dashed = threshold', fontsize=10, loc='left')
    fig.tight_layout()
    return fig


def plot_metric_scatter(metrics, reject_flag, reject_method, thresholds, title_prefix=''):
    """Scatter p-p vs max-gradient, coloured clean/grey or by reject method. Returns a Figure."""
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


def build_rejection_table_html(meta, methods_present, stage_order, title=''):
    """Per-stage × per-method rejection table (% and raw count), from the .fif metadata flags."""
    rows = []
    header = ['Stage', 'N epochs'] + [METHOD_LABEL[m] for m in methods_present] + ['Any']
    stages = meta['stage'].astype(str).values
    for st in stage_order:
        sel = (stages == st)
        n = int(sel.sum())
        if n == 0:
            continue
        cells = [st, str(n)]
        for m in methods_present:
            c = int(meta.loc[sel, 'flag_' + m].astype(bool).sum())
            cells.append(f'{100.0 * c / n:.1f}% ({c})')
        a = int(meta.loc[sel, 'reject_flag'].astype(bool).sum())
        cells.append(f'<b>{100.0 * a / n:.1f}% ({a})</b>')
        rows.append(cells)
    # total row
    n = len(meta)
    tot = ['<b>All</b>', f'<b>{n}</b>']
    for m in methods_present:
        c = int(meta['flag_' + m].astype(bool).sum())
        tot.append(f'<b>{100.0 * c / n:.1f}% ({c})</b>')
    a = int(meta['reject_flag'].astype(bool).sum())
    tot.append(f'<b>{100.0 * a / n:.1f}% ({a})</b>')
    rows.append(tot)
    th = ''.join(f'<th style="padding:3px 8px;border:1px solid #ccc;">{h}</th>' for h in header)
    body = ''
    for r in rows:
        body += '<tr>' + ''.join(
            f'<td style="padding:3px 8px;border:1px solid #ccc;text-align:center;">{c}</td>' for c in r) + '</tr>'
    return (f'<h4>{title}</h4>' if title else '') + \
           f'<table style="border-collapse:collapse;font-size:.85em;"><tr>{th}</tr>{body}</table>'


def build_participant_report_figs(P, metrics, freqs, psds_uV2, thresholds, custom_stages):
    """Assemble the Section-2 figures for one participant. Returns (list_of_(title, fig), table_html)."""
    meta = P['meta']
    stages = P['stages']
    reject_flag = P['reject_flag']
    reject_method = meta['reject_method'].astype(str).values
    stage_order = ordered_present_stages(stages, custom_stages)
    figs = []
    for st in stage_order:
        figs.append((f'PSD overlay — {st}',
                     plot_psd_overlay(freqs, psds_uV2, stages, reject_flag, reject_method, st,
                                      custom_stages=custom_stages)))
    figs.append(('Metric distributions',
                 plot_metric_distributions(metrics, stages, reject_flag, thresholds, stage_order)))
    figs.append(('Metric scatter',
                 plot_metric_scatter(metrics, reject_flag, reject_method, thresholds)))
    table_html = build_rejection_table_html(meta, P['methods_present'], stage_order,
                                            title='Rejection by stage × method')
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


def plot_epoch_montage(P, ei, thresholds, context=1, ctx=None):
    """Stacked EEG montage of epoch `ei` (± `context` epochs) with the current epoch shaded, the
    current-epoch trace of each channel coloured by its recomputed flagging method, the gradient-max
    sample marked, and channel labels coloured by method. Returns a Figure.

    `ctx` (optional) is the dict from `load_context_epochs` — EOG/EMG/ECG traces epoched 1:1 with the
    EEG. When given, those channels are stacked below the EEG montage (own scale, dotted divider) to
    help read eye movements / muscle bursts. `ctx=None` keeps the plot EEG-only."""
    data, sf, ch = P['data_uV'], P['sfreq'], P['ch_names']
    stages, meta = P['stages'], P['meta']
    n_ep, n_ch, n_t = data.shape
    lo, hi = max(0, ei - context), min(n_ep - 1, ei + context)
    idxs = list(range(lo, hi + 1))
    seg = np.concatenate([data[k] for k in idxs], axis=1)       # (n_ch, n_t*len)
    t = np.arange(seg.shape[1]) / sf
    spacing = max(float(np.percentile(np.ptp(seg, axis=1), 90)), 50.0) * 1.15
    n_ctx = len(ctx['labels']) if ctx is not None else 0
    fig, ax = plt.subplots(figsize=(11, 1.0 * n_ch + 1.4 + 0.9 * n_ctx))
    cur_pos = idxs.index(ei)
    ax.axvspan(cur_pos * n_t / sf, (cur_pos + 1) * n_t / sf, color='#fff3cd', alpha=0.7, zorder=0)
    stg = stages[ei]
    cur = data[ei]
    for i, cn in enumerate(ch):
        off = (n_ch - 1 - i) * spacing
        ax.plot(t, seg[i] + off, color='0.55', lw=0.5, zorder=2)
        m, ptp, grad = _epoch_channel_method(cur[i], stg, thresholds)
        col = '0.2' if m is None else METHOD_COLOR[m]
        tt = np.arange(cur_pos * n_t, (cur_pos + 1) * n_t) / sf
        ax.plot(tt, cur[i] + off, color=col, lw=0.8, zorder=3)
        if grad > thresholds['gradient_uV_per_sample']:
            j = int(np.argmax(np.abs(np.diff(cur[i]))))
            ax.plot(tt[j], cur[i][j] + off, 'o', color=METHOD_COLOR['gradient'], ms=5, zorder=4)
        ax.text(-0.008, off, cn, ha='right', va='center', fontsize=8,
                color=('k' if m is None else METHOD_COLOR[m]),
                transform=ax.get_yaxis_transform())
    # --- optional EOG/EMG/ECG context traces, stacked below the EEG montage (own scale) ---
    if ctx is not None and n_ctx > 0:
        cdata, csf, clabels = ctx['data_uV'], ctx['sfreq'], ctx['labels']
        n_cep = cdata.shape[0]
        cidxs = [k for k in idxs if k < n_cep]                    # align by epoch index
        if cidxs:
            cseg = np.concatenate([cdata[k] for k in cidxs], axis=1)   # (n_cctx, n_ct*len)
            tctx = np.arange(cseg.shape[1]) / csf
            cspacing = max(float(np.percentile(np.ptp(cseg, axis=1), 90)), 50.0) * 1.15
            ax.axhline(-0.6 * spacing, color='0.4', lw=0.8, ls=':', zorder=1)   # EEG | context divider
            for j, cl in enumerate(clabels):
                coff = -0.6 * spacing - (j + 1) * cspacing
                col = CTX_COLOR.get(cl, '#8e44ad')
                ax.plot(tctx, cseg[j] + coff, color=col, lw=0.6, zorder=2)
                ax.text(-0.008, coff, cl, ha='right', va='center', fontsize=8, color=col,
                        transform=ax.get_yaxis_transform())
    for k in range(len(idxs) + 1):
        ax.axvline(k * n_t / sf, color='0.85', lw=0.6, zorder=1)
    ax.set_yticks([])
    ax.set_xlabel('Time (s)')
    ax.set_xlim(t[0], t[-1])
    rm = meta['reject_method'].astype(str).values[ei]
    ax.set_title(f'Epoch {ei} — stage {stg} — reject_method={rm}   '
                 f'(context ±{context}; {spacing:.0f} µV/row; coloured = flagging method)', fontsize=10)
    fig.tight_layout()
    return fig


def plot_epoch_detail(P, ei, freqs, psds_uV2, thresholds, spectro_ch_idx=0, fmin=2.0):
    """2×2 detail for epoch `ei`: PSD + aperiodic fit (worst-R² channel), mean band power + 50 Hz
    ratio, epoch spectrogram, per-channel metric table (value vs threshold). Returns a Figure."""
    from scipy.signal import welch as _welch
    data, sf, ch = P['data_uV'], P['sfreq'], P['ch_names']
    stg = P['stages'][ei]
    n_ch = len(ch)
    fits = [fit_1f(freqs, psds_uV2[ei, c], fmin=fmin) for c in range(n_ch)]
    r2s = np.array([f[1] for f in fits])
    wc = int(np.nanargmin(r2s)) if not np.all(np.isnan(r2s)) else 0
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # (a) PSD + aperiodic fit
    ax = axes[0, 0]
    for c in range(n_ch):
        ax.loglog(freqs, psds_uV2[ei, c], color='0.82', lw=0.6)
    ax.loglog(freqs, psds_uV2[ei, wc], color='#2c7fb8', lw=1.3, label=f'{ch[wc]} (worst R²)')
    mae, r2, off, exp, _sm = fits[wc]
    if not np.isnan(off):
        fm = freqs >= 2
        ax.loglog(freqs[fm], 10 ** (off - exp * np.log10(freqs[fm])),
                  color='#c0392b', lw=1.3, ls='--', label='aperiodic fit')
    ax.set_title(f'PSD — worst-R² channel {ch[wc]}: MAE={mae:.3f}, R²={r2:.3f}, exp={exp:.2f}',
                 fontsize=9)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('PSD (µV²/Hz)')
    ax.legend(fontsize=7)

    # (b) mean band power + 50 Hz ratio
    ax = axes[0, 1]
    bp = band_powers(freqs, psds_uV2[ei].mean(axis=0))
    ax.bar([b[0] for b in BANDS], [bp[b[0]] for b in BANDS], color='#807dba')
    ax.set_ylabel('Band power (µV²)')
    ax.set_yscale('log')
    ax.set_title('Mean band power (across channels)', fontsize=9)
    if sf / 2 > 52:
        f2, p2 = _welch(data[ei, wc], fs=sf, nperseg=min(int(4 * sf), data.shape[-1]))
        def _bp(lo, hi):
            m = (f2 >= lo) & (f2 < hi)
            return float(np.trapz(p2[m], f2[m])) if m.sum() > 1 else np.nan
        nb = _bp(40, 47)
        ratio = _bp(48, 52) / nb if nb and nb > 0 else np.nan
        ax.text(0.98, 0.96, f'50 Hz / 40–47 Hz = {ratio:.2f}  ({ch[wc]})',
                transform=ax.transAxes, ha='right', va='top', fontsize=8)
    else:
        ax.text(0.98, 0.96, '50 Hz ratio n/a (Nyquist ≤ 52 Hz)',
                transform=ax.transAxes, ha='right', va='top', fontsize=8)

    # (c) epoch spectrogram
    ax = axes[1, 0]
    f3, t3, Sxx = sp_spectrogram(data[ei, spectro_ch_idx], fs=sf,
                                 nperseg=int(1.5 * sf), noverlap=int(0.75 * sf))
    fm = f3 <= 40
    ax.pcolormesh(t3, f3[fm], 10 * np.log10(Sxx[fm] + 1e-12), shading='auto', cmap='viridis')
    ax.set_ylabel('Frequency (Hz)')
    ax.set_xlabel('Time (s)')
    ax.set_title(f'Epoch spectrogram — {ch[spectro_ch_idx]}', fontsize=9)

    # (d) per-channel metric table
    ax = axes[1, 1]
    ax.set_axis_off()
    amp_thr = thresholds['amplitude_ptp_uV'].get(stg, 250.0)
    cur = data[ei]
    table = [['ch', 'p-p (thr)', 'grad (thr)', 'MAE (thr)', 'R² (thr)']]
    for c in range(n_ch):
        ptp = float(np.ptp(cur[c]))
        grad = float(np.max(np.abs(np.diff(cur[c]))))
        mae_c, r2_c = fits[c][0], fits[c][1]
        table.append([ch[c],
                      f'{ptp:.0f} ({amp_thr:.0f})',
                      f'{grad:.0f} ({thresholds["gradient_uV_per_sample"]:.0f})',
                      f'{mae_c:.3f} ({thresholds["1f_mae_max"]:.2f})',
                      f'{r2_c:.3f} ({thresholds["1f_r2_min"]:.2f})'])
    tb = ax.table(cellText=table, loc='center', cellLoc='center')
    tb.auto_set_font_size(False)
    tb.set_fontsize(8)
    tb.scale(1, 1.5)
    ax.set_title(f'Per-channel metrics (stage {stg})', fontsize=9)
    fig.tight_layout()
    return fig


def plot_review_strip(P, final_reject, overridden, custom_stages=()):
    """Section-4 review strip: hypnogram step-line on top; below, a per-epoch bar coloured by the
    final keep/reject decision, with overridden epochs outlined. Returns a Figure."""
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
    ax.bar(np.arange(n_ep)[~fr], np.ones((~fr).sum()), width=1.0, color='#2ecc71', align='edge')
    ax.bar(np.arange(n_ep)[fr], np.ones(fr.sum()), width=1.0, color='#c0392b', align='edge')
    for ei in sorted(overridden):
        ax.plot([ei + 0.5], [1.15], marker='v', color='k', ms=4)
    ax.set_ylim(0, 1.35)
    ax.set_yticks([])
    ax.set_xlim(0, n_ep)
    ax.set_xlabel('Epoch index', fontsize=9)
    ax.set_title(f'Final decision — green = keep ({int((~fr).sum())}), red = reject ({int(fr.sum())}), '
                 f'▼ = manually overridden ({len(overridden)})', fontsize=9)
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    return fig


def save_report_html(out_path, title, figs, table_html):
    """Write an mne.Report HTML with the Section-2 figures + rejection table."""
    report = mne.Report(title=title, verbose=False)
    report.add_html(html=table_html, title='Rejection summary', section='Overview')
    for name, fig in figs:
        report.add_figure(fig=fig, title=name, section='Per-stage', image_format='PNG')
        plt.close(fig)
    report.save(str(out_path), overwrite=True, open_browser=False, verbose=False)
