"""Preprocessing and flagging of tool 7 (7_preprocessing_voila), shared with tool 0 and with tools 8 / 8bis / 9
through qc_rejected_epochs_lib (the PSD smoothing).

Moved verbatim from tool 7 (its shared-function cell, and the computation lines of its run loop): tool 7 and
tool 0's step 6 call the same code, so they give the same flags by construction. Edit here, never in a
notebook. The Curry twin of tool 7 imports this file as well (everything here is format-agnostic: reading the
recording, the scored events and the context channels stays in the notebooks, where the Curry generator
rewrites it).

METHOD_ORDER is the single source of truth of the flagging methods (order, codes, labels).
"""
import json
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
from specparam import SpectralModel
from statsmodels.nonparametric.smoothers_lowess import lowess


# ---------------------------------------------------------------------------------------------------------
# Methods, stages, flagging, tables and figures (tool 7's shared-function cell)
# ---------------------------------------------------------------------------------------------------------
# ---- Rejection-method registry (single source of truth) ----
# 'event' is the optional 6th method (epochs containing selected scored events). The helpers
# below iterate METHOD_ORDER so 'event' is handled when present and silently ignored when absent.
METHOD_ORDER  = ['amplitude', 'gradient', 'flat', '1f_r2', '1f_error', 'event']
METHOD_CODE   = {m: i + 1 for i, m in enumerate(METHOD_ORDER)}   # 1..6  (0 = none)
MULTIPLE_CODE = len(METHOD_ORDER) + 1                            # 7 = multiple
# Human-readable method names for report table headers (keyed by METHOD_ORDER labels).
METHOD_LABEL  = {'amplitude': 'Amplitude', 'gradient': 'Gradient', 'flat': 'Flat',
                 '1f_r2': '1/f R²', '1f_error': '1/f error', 'event': 'Event'}


# --- Custom (non-AASM) sleep stages -------------------------------------------
# Declared in config_param/custom_stages.json (written by 3_remap_hypno) and recognised here
# so per-stage rejection thresholds, summaries and the heatmap include them.
BASE_STAGE_COLORS = {'W': '#969696', 'N1': '#9e9ac8', 'N2': '#807dba', 'N3': '#6a51a3', 'R': '#c994c7'}
# Distinct, non-red palette for custom stages (red is reserved for REM on the hypnogram).
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


def compute_event_epoch_mask(events_df, event_remap, selected_types, n_epochs,
                             epoch_sec=30.0):
    """Flag each epoch (epoch_sec long) that CONTAINS the onset of any selected canonical event type.

    events_df      : DataFrame Name/Start/Duration (seconds) from load_events, or None
    event_remap    : {raw_label: canonical_label_or_None} from event_remap.json
    selected_types : set/list of canonical labels chosen by the user
    Returns a (n_epochs,) bool array (epoch-level; replicated across channels later).
    Containment rule: an epoch is flagged when an event's ONSET (Start) falls inside it.
    The annotated Duration is intentionally IGNORED -- clinicians often score only the
    event onset, not a precise duration, so the duration field is not trusted to widen
    the rejection. Each matching event therefore flags exactly one epoch.
    """
    mask = np.zeros(n_epochs, dtype=bool)
    selected = set(selected_types)
    if events_df is None or not selected:
        return mask
    for _, ev in events_df.iterrows():
        canonical = event_remap.get(str(ev['Name']))
        if canonical is None or canonical not in selected:
            continue
        start = float(ev['Start'])
        if not np.isfinite(start) or start < 0:
            continue
        e = int(np.floor(start / epoch_sec))
        if 0 <= e < n_epochs:
            mask[e] = True
    return mask


def smooth_psd_median(psds_uV2, freqs, span_hz=3.0):
    """Running-median smooth of each (epoch, channel) linear PSD along frequency.
    Reproduces oscip.smooth_spectrum_median (MATLAB movmedian over a span given in Hz): a centred
    median filter that removes narrow spikes (residual line-noise harmonics, single-bin artefacts)
    BEFORE the LOWESS mean smoothing. The Hz span is converted to a point count
    (span_pts = round(span_hz / freq_res)); pandas rolling(center=True, min_periods=1) matches
    MATLAB\'s truncated-window edge behaviour. Non-fatal: a failing (epoch, channel) is left as-is."""
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
    Reproduces oscip.smooth_spectrum from the Snipes MATLAB pipeline: MATLAB\'s non-robust
    \'lowess\' (local linear regression, tricube weights) over a span given in Hz, applied to the
    LINEAR power (matching oscip; specparam logs it internally afterwards). The Hz span is converted
    to a fraction of points because statsmodels expresses the window as a fraction (frac), whereas
    MATLAB uses a point count. it=0 = non-robust (matches \'lowess\', not \'rlowess\').
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


def compute_rejection_masks(epochs, sfreq, hypno_epochs, ptp_thresholds,
                             flat_thresh, grad_thresh,
                             psd_freqs, psds_data_uV2, thresh_error, thresh_r2,
                             event_epoch_mask=None, fit_fmin=2.0, fit_fmax=45.0, progress=None,
                             do_amplitude=True, do_flat=True, do_gradient=True):
    """
    Compute per-(epoch, channel) boolean rejection masks for each method.

    Per-channel variant (memory): reads the signal ONE CHANNEL AT A TIME from `epochs`
    (epochs.get_data(picks=[ci])) so a high-density montage is never materialised as a second
    full (n_epochs, n_channels, n_times) float64 array. The formulas are identical to the former
    full-array version, so the per-(epoch, channel) masks are byte-identical.

    progress(done, total, msg) is called periodically for a UI progress bar (per channel for the
    time-domain pass, then per epoch for the 1/f fits).

    Returns dict of (n_epochs, n_channels) bool arrays, one key per ENABLED method:
        'amplitude'/'flat'/'gradient' when do_amplitude/do_flat/do_gradient (default True),
        '1f_error'/'1f_r2' when a PSD is provided (psds_data_uV2 not None), plus 'event' if provided.
    """
    n_epochs = len(epochs)
    n_channels = len(epochs.ch_names)

    ptp_thresh_arr = np.array([ptp_thresholds.get(s, 250.0) for s in hypno_epochs])  # (n_epochs,)

    mask_amplitude = np.zeros((n_epochs, n_channels), dtype=bool)
    mask_flat      = np.zeros((n_epochs, n_channels), dtype=bool)
    mask_gradient  = np.zeros((n_epochs, n_channels), dtype=bool)

    do_1f = psds_data_uV2 is not None and psd_freqs is not None
    do_time = do_amplitude or do_flat or do_gradient   # any time-domain method enabled?
    total = (n_channels if do_time else 0) + (n_epochs if do_1f else 0)
    done = 0

    # Time-domain methods, one channel at a time — only the enabled ones are computed (a disabled
    # method keeps its all-False mask and is left out of mask_dict below).
    if do_time:
        for ci in range(n_channels):
            ch = epochs.get_data(picks=[ci])[:, 0, :] * 1e6   # (n_epochs, n_times) in uV
            if do_amplitude or do_flat:
                ptp = np.ptp(ch, axis=-1)
            if do_amplitude:
                mask_amplitude[:, ci] = ptp > ptp_thresh_arr
            if do_flat:
                mask_flat[:, ci]      = ptp < flat_thresh
            if do_gradient:
                mask_gradient[:, ci]  = np.max(np.abs(np.diff(ch, axis=-1)), axis=-1) > grad_thresh
            del ch
            done += 1
            if progress is not None:
                progress(done, total, 'scanning channels')

    # 1/f fit quality via specparam (per epoch x channel) — same model as the full-array version.
    mask_1f_error = np.zeros((n_epochs, n_channels), dtype=bool)
    mask_1f_r2    = np.zeros((n_epochs, n_channels), dtype=bool)
    if do_1f:
        freq_mask = (psd_freqs >= fit_fmin) & (psd_freqs <= fit_fmax)
        freqs_reduced = psd_freqs[freq_mask]
        sp_model = SpectralModel(
            peak_width_limits=[0.5, 20],
            aperiodic_mode='fixed',
            min_peak_height=0.3,
            max_n_peaks=8,   # cap the peak search: bounds specparam's runaway peak-fitting on
                             # noisy short-epoch PSDs (8 >= the realistic neural peak count on 2-45 Hz,
                             # so it is NOT the max_n_peaks=0 case that would over-reject N2/REM).
        )
        for ei in range(n_epochs):
            for ci in range(n_channels):
                psd_reduced = psds_data_uV2[ei, ci, freq_mask]
                try:
                    sp_model.fit(freqs_reduced, psd_reduced)
                    error = sp_model.get_metrics('error', 'mae')
                    r2    = sp_model.get_metrics('gof', 'squared')
                    if error > thresh_error:
                        mask_1f_error[ei, ci] = True
                    if r2 < thresh_r2:
                        mask_1f_r2[ei, ci] = True
                except Exception:
                    mask_1f_error[ei, ci] = True
                    mask_1f_r2[ei, ci]    = True
            done += 1
            if progress is not None and (ei % 25 == 0 or ei == n_epochs - 1):
                progress(done, total, '1/f fit')

    # Only enabled methods enter mask_dict, so every downstream table/summary that iterates
    # [m for m in METHOD_ORDER if m in mask_dict] drops a disabled method automatically (like 'event').
    mask_dict = {}
    if do_amplitude:
        mask_dict['amplitude'] = mask_amplitude
    if do_flat:
        mask_dict['flat'] = mask_flat
    if do_gradient:
        mask_dict['gradient'] = mask_gradient
    if do_1f:
        mask_dict['1f_error'] = mask_1f_error
        mask_dict['1f_r2']    = mask_1f_r2
    if event_epoch_mask is not None:
        mask_dict['event'] = np.tile(np.asarray(event_epoch_mask, dtype=bool)[:, np.newaxis],
                                     (1, n_channels))
    return mask_dict


def build_combined_method_matrix(mask_dict):
    """
    Build an (n_epochs, n_channels) integer matrix encoding which method first
    flagged each cell.  0=none, then METHOD_CODE (1..6), MULTIPLE_CODE (7) for multi-flag cells.
    """
    methods   = [m for m in METHOD_ORDER if m in mask_dict]
    shape     = next(iter(mask_dict.values())).shape
    combined  = np.zeros(shape, dtype=int)
    for name in methods:
        combined[mask_dict[name]] = METHOD_CODE[name]
    n_flags = sum(mask_dict[name].astype(int) for name in methods)
    combined[n_flags > 1] = MULTIPLE_CODE  # overwrite with 'multiple'
    return combined


def plot_rejection_heatmap(combined_matrix, ch_names, hypno_epochs, file_id, custom_stages=()):
    """
    Plot a channels x epochs heatmap coloured by rejection method.
    Top strip: hypnogram as a YASA-style step-line (W at top, N3 at bottom).
    Gray line for all stages; REM highlighted in red.
    Returns a matplotlib Figure.
    """
    # Index = method code: 0 none, 1 amplitude, 2 gradient, 3 flat, 4 1/f R2, 5 1/f error,
    # 6 event, 7 multiple (matches METHOD_CODE / MULTIPLE_CODE). 'event' colour unused when off.
    # CVD-validated palette (in sync with qc_rejected_epochs_lib.HEATMAP_COLORS): 6 method hues picked
    # for colour-blind separation; 'multiple' = cyan #22d3d3, distinct from all six on the dark bg.
    COLORS = ['#1c0a3b', '#e34948', '#4a3aa7', '#2a78d6', '#008300', '#eda100', '#e87ba4', '#22d3d3']
    LABELS = ['none', 'amplitude', 'gradient', 'flat', '1/f R2', '1/f error', 'event', 'multiple']
    cmap   = mcolors.ListedColormap(COLORS)
    bounds = np.arange(-0.5, 8.5, 1)
    norm   = mcolors.BoundaryNorm(bounds, cmap.N)

    # Hypnogram heights: W at top (4), N3 at bottom (0); custom stages stacked below N3.
    STAGE_H, _stage_colors, _ytick_pos, _ytick_labels = custom_stage_style(custom_stages)
    _floor = min(_ytick_pos) - 1
    hypno_y = np.array([STAGE_H.get(s, _floor) for s in hypno_epochs])

    n_channels   = len(ch_names)
    n_epochs     = combined_matrix.shape[0]
    # Pair-based: fraction of (epoch x channel) cells flagged (not '>=1 bad channel per epoch',
    # which over-counts). Matches the per-stage table denominator (n_epochs * n_channels).
    pct_rejected = 100.0 * (combined_matrix > 0).mean()

    fig_w = max(8, min(n_epochs / 6, 20))
    # Fixed-ish height: rows get thinner as the montage grows (capped so a 32-channel Curry
    # montage stays compact instead of a ~16 in tall figure). Unchanged for <=14 channels.
    fig_h = max(3, min(0.45 * n_channels + 2.0, 8.0))
    fig, axes = plt.subplots(
        2, 1, figsize=(fig_w, fig_h),
        gridspec_kw={'height_ratios': [1.2, n_channels], 'hspace': 0.04}
    )

    # --- Hypnogram strip ---
    ax_hyp = axes[0]
    x_step = np.arange(n_epochs + 1)
    y_step = np.append(hypno_y, hypno_y[-1])

    # Full step line in gray
    ax_hyp.step(x_step, y_step, where='post', color='#555555', linewidth=1.2)

    # REM in red; each custom stage highlighted in its own colour (at the line's width: the toolkit's
    # discreet convention).
    for _s in (['R'] + list(custom_stages)):
        _col = '#c0392b' if _s == 'R' else _stage_colors.get(_s)
        for ei in range(n_epochs):
            if hypno_epochs[ei] == _s:
                ax_hyp.plot([ei, ei + 1], [STAGE_H[_s], STAGE_H[_s]],
                            color=_col, linewidth=1.2, solid_capstyle='butt')

    ax_hyp.set_yticks(_ytick_pos)
    ax_hyp.set_yticklabels(_ytick_labels, fontsize=7)
    ax_hyp.set_ylim(min(_ytick_pos) - 0.8, max(_ytick_pos) + 0.8)
    ax_hyp.set_xlim(0, n_epochs)
    ax_hyp.set_xticks([])
    for sp in ['top', 'right', 'bottom']:
        ax_hyp.spines[sp].set_visible(False)

    # --- Main heatmap ---
    ax = axes[1]
    im = ax.pcolormesh(
        np.arange(n_epochs + 1), np.arange(n_channels + 1),
        combined_matrix.T,  # (n_channels, n_epochs)
        cmap=cmap, norm=norm, rasterized=True
    )
    ax.set_yticks(np.arange(n_channels) + 0.5)
    ax.set_yticklabels(ch_names, fontsize=max(5, min(8, 200 / n_channels)))
    ax.set_xlabel('Epoch index', fontsize=9)
    ax.set_xlim(0, n_epochs)
    ax.set_ylim(0, n_channels)
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)

    cbar = fig.colorbar(im, ax=ax, orientation='vertical',
                        fraction=0.02, pad=0.01, ticks=np.arange(len(LABELS)))
    cbar.ax.set_yticklabels(LABELS, fontsize=7)

    fig.suptitle(f'{file_id}: {pct_rejected:.0f}% of flagged epochs',
                 fontsize=11, y=1.01)
    fig.tight_layout()
    return fig


def build_rejection_summary(mask_dict, hypno_epochs, file_id, ch_names, custom_stages=()):
    """Per-stage, per-method rejection counts DataFrame. Counts are over (epoch, channel) PAIRS:
    n_total = n_epochs_in_stage * n_channels and n_rejected = number of flagged pairs, so
    pct_rejected is the fraction of that stage's (epoch x channel) cells that were flagged. This
    replaces the old '>=1 bad channel => whole epoch rejected' rule, which reached ~100% as soon
    as every epoch had one bad channel (even on different channels) and so over-counted."""
    methods    = [m for m in METHOD_ORDER if m in mask_dict]
    n_channels = len(ch_names)
    rows       = []
    for stage in ['W', 'N1', 'N2', 'N3', 'R'] + list(custom_stages):
        stage_mask = (hypno_epochs == stage)
        n_stage_ep = int(stage_mask.sum())
        n_total    = n_stage_ep * n_channels                          # (epoch x channel) pairs
        any_pair   = np.zeros((n_stage_ep, n_channels), dtype=bool)   # OR of methods, per pair
        for name in methods:
            any_pair |= mask_dict[name][stage_mask, :]
        n_rej_any = int(any_pair.sum())
        rows.append({
            'file_id': file_id, 'stage': stage, 'method': 'any',
            'n_total': n_total, 'n_rejected': n_rej_any,
            'pct_rejected': round(100 * n_rej_any / n_total, 1) if n_total > 0 else float('nan'),
        })
        for name in methods:
            n_rej = int(mask_dict[name][stage_mask, :].sum())         # flagged pairs for this method
            rows.append({
                'file_id': file_id, 'stage': stage, 'method': name,
                'n_total': n_total, 'n_rejected': n_rej,
                'pct_rejected': round(100 * n_rej / n_total, 1) if n_total > 0 else float('nan'),
            })
    return pd.DataFrame(rows)


def build_epoch_rejection_log(mask_dict, hypno_epochs, file_id):
    """Per-epoch rejection log: one row per epoch with epoch-level (any-channel) method flags."""
    methods  = [m for m in METHOD_ORDER if m in mask_dict]
    n_epochs = len(hypno_epochs)
    rows = []
    for ei in range(n_epochs):
        row = {
            'file_id':     file_id,
            'epoch_idx':   ei,
            'stage':       hypno_epochs[ei],
            'reject_flag': bool(any(mask_dict[m][ei].any() for m in methods)),
        }
        for m in methods:
            col = 'flag_' + m.replace('/', '').replace(' ', '')
            row[col] = bool(mask_dict[m][ei].any())
        rows.append(row)
    return pd.DataFrame(rows)


def build_epoch_channel_rejection_log(mask_dict, hypno_epochs, ch_names, file_id):
    """Per-(epoch, channel) rejection log: one row per (epoch, channel) pair, one 0/1 column per
    method (event is replicated across all channels of an epoch — mask_dict['event'] is already
    tiled), plus reject_any (OR over the methods present). Persists the fine granularity the
    per-epoch _epoch_rejection.tsv collapses away, so a downstream decision tool can rebuild any
    per-channel / per-epoch rejection rule without reloading the signal. Rows are epoch-major
    (all channels of epoch 0, then epoch 1, ...). Built from the (n_epochs, n_channels) masks."""
    methods    = [m for m in METHOD_ORDER if m in mask_dict]
    n_epochs   = len(hypno_epochs)
    n_channels = len(ch_names)
    cols = {
        'file_id':   file_id,
        'epoch_idx': np.repeat(np.arange(n_epochs), n_channels),
        'stage':     np.repeat(np.asarray(hypno_epochs), n_channels),
        'channel':   np.tile(np.asarray(ch_names), n_epochs),
    }
    any_pair = np.zeros((n_epochs, n_channels), dtype=bool)
    for m in methods:
        col = 'flag_' + m.replace('/', '').replace(' ', '')
        cols[col] = mask_dict[m].reshape(-1).astype(int)
        any_pair |= mask_dict[m]
    cols['reject_any'] = any_pair.reshape(-1).astype(int)
    return pd.DataFrame(cols)


def build_global_summary_table(all_rejection_summaries, custom_stages=()):
    """
    Build a dataset-level pivot table from per-participant rejection summaries.
    One row per sleep stage; columns: n_total, n_rejected, pct_rejected (pooled),
    then n_rej_{method} and pct_rej_{method} for each rejection method. Counts are
    over (epoch x channel) PAIRS (n_total = sum of n_epochs_in_stage * n_channels
    across participants), inherited from build_rejection_summary — no code change needed,
    the pooled pct is automatically pair-based.
    """
    df_all  = pd.concat(all_rejection_summaries, ignore_index=True)
    # Methods present in the loaded summaries (so 'event' columns appear only when used and
    # mixed runs — some with, some without event rejection — never raise).
    methods_list = [m for m in METHOD_ORDER if m in set(df_all['method'])]
    rows = []
    for stage in ['W', 'N1', 'N2', 'N3', 'R'] + list(custom_stages):
        stage_data = df_all[df_all['stage'] == stage]
        any_data   = stage_data[stage_data['method'] == 'any']
        n_total    = int(any_data['n_total'].sum())
        n_rej      = int(any_data['n_rejected'].sum())
        n_parts    = int(len(any_data))
        row = {
            'stage':          stage,
            'n_participants': n_parts,
            'n_total':        n_total,
            'n_rejected':     n_rej,
            'pct_rejected':   round(100 * n_rej / n_total, 1) if n_total > 0 else float('nan'),
        }
        for m in methods_list:
            m_data  = stage_data[stage_data['method'] == m]
            n_m_rej = int(m_data['n_rejected'].sum())
            row[f'n_rej_{m}']   = n_m_rej
            row[f'pct_rej_{m}'] = round(100 * n_m_rej / n_total, 1) if n_total > 0 else float('nan')
            # Canary for mixed-method databases: how many participants actually ran method m. When it
            # is < n_participants, pct_rej_{m} is diluted (its denominator n_total spans ALL files).
            row[f'n_part_{m}']  = int(m_data['file_id'].nunique())
        rows.append(row)
    return pd.DataFrame(rows)


def build_stage_method_html(rej_summary):
    """Per-stage flagging table for the HTML report: one row/stage, one column/method.

    The 'epochs×channels' column is the any-method flagged rate for that stage; each method
    column is that method's rate. Every cell = % (in front) of that stage's (epoch × channel)
    PAIRS flagged, with the raw pair count in parentheses. 'event' appears only when event
    flagging ran (methods present in rej_summary). % is relative to the stage total n_total =
    n_epochs_in_stage * n_channels, matching global_rejection_by_stage.tsv.
    """
    if rej_summary is None or rej_summary.empty:
        return '<p><em>Per-stage summary unavailable.</em></p>'
    methods = [m for m in METHOD_ORDER if m in set(rej_summary['method'])]
    def _cell(pct, n):
        if pd.isna(pct):
            return '<td style="padding:3px 10px;text-align:right;">—</td>'
        return (f'<td style="padding:3px 10px;text-align:right;"><b>{pct:.1f}%</b>'
                f'<small style="color:#888;"> ({int(n)})</small></td>')
    head = ('<tr><th style="padding:3px 10px;text-align:left;">Stage</th>'
            '<th style="padding:3px 10px;text-align:right;">epochs&times;channels</th>'
            + ''.join(f'<th style="padding:3px 10px;text-align:right;">{METHOD_LABEL.get(m, m)}</th>'
                      for m in methods) + '</tr>')
    body = ''
    for stage in list(rej_summary['stage'].drop_duplicates()):
        rows_s  = rej_summary[rej_summary['stage'] == stage]
        any_row = rows_s[rows_s['method'] == 'any']
        if any_row.empty:
            continue
        body += (f'<tr><td style="padding:3px 10px;">{stage}</td>'
                 + _cell(any_row['pct_rejected'].iloc[0], any_row['n_rejected'].iloc[0]))
        for m in methods:
            m_row = rows_s[rows_s['method'] == m]
            body += (_cell(m_row['pct_rejected'].iloc[0], m_row['n_rejected'].iloc[0])
                     if not m_row.empty else '<td style="padding:3px 10px;text-align:right;">—</td>')
        body += '</tr>'
    return f'<table style="border-collapse:collapse;font-size:.9em;">{head}{body}</table>'


def build_channel_stage_html(mask_dict, hypno_epochs, ch_names, custom_stages=()):
    """Per-channel × per-stage rejection table for the HTML report: channels as rows, stages
    as columns. Each cell = % of that stage's epochs where the channel was flagged by ANY
    method, with the raw count in parentheses; the last 'All' column pools all stages. Built
    from the per-(epoch, channel) masks (not persisted to TSV) — report-only."""
    methods = [m for m in METHOD_ORDER if m in mask_dict]
    if not methods:
        return '<p><em>Per-channel summary unavailable.</em></p>'
    hyp = np.asarray(hypno_epochs)
    any_ec = np.zeros(mask_dict[methods[0]].shape, dtype=bool)   # (n_epochs, n_channels) any-method
    for m in methods:
        any_ec |= mask_dict[m]
    stages = ['W', 'N1', 'N2', 'N3', 'R'] + list(custom_stages)
    def _cell(n_rej, n_tot):
        if not n_tot:
            return '<td style="padding:3px 10px;text-align:right;">—</td>'
        return (f'<td style="padding:3px 10px;text-align:right;"><b>{100 * n_rej / n_tot:.1f}%</b>'
                f'<small style="color:#888;"> ({int(n_rej)})</small></td>')
    head = ('<tr><th style="padding:3px 10px;text-align:left;">Channel</th>'
            + ''.join(f'<th style="padding:3px 10px;text-align:right;">{s}</th>' for s in stages)
            + '<th style="padding:3px 10px;text-align:right;">All</th></tr>')
    body = ''
    for ci, ch in enumerate(ch_names):
        cells = ''
        for s in stages:
            stage_mask = (hyp == s)
            cells += _cell(int((any_ec[:, ci] & stage_mask).sum()), int(stage_mask.sum()))
        cells += _cell(int(any_ec[:, ci].sum()), len(hyp))
        body += f'<tr><td style="padding:3px 10px;">{ch}</td>{cells}</tr>'
    return f'<table style="border-collapse:collapse;font-size:.9em;">{head}{body}</table>'


def build_event_type_html(event_type_masks, hypno_epochs, custom_stages=()):
    """Per-event-type rejection table for the HTML report (event rejection only).

    One row per selected canonical event type + a bold '(any selected)' union row; columns:
    n epochs flagged, % of all epochs, then % (n) per sleep stage. Helps spot event types that
    reject too many epochs. Onset-only containment, same rule as compute_event_epoch_mask.
    """
    if not event_type_masks:
        return None
    stages   = ['W', 'N1', 'N2', 'N3', 'R'] + list(custom_stages)
    hypno    = np.asarray(hypno_epochs)
    n_epochs = len(hypno)
    stage_totals = {s: int((hypno == s).sum()) for s in stages}
    def _pct(n, tot):
        return f'{100 * n / tot:.1f}%' if tot else '—'
    def _row(label, mask, bold=False):
        mask = np.asarray(mask, dtype=bool)
        n    = int(mask.sum())
        weight = 'font-weight:600;' if bold else ''
        cells = ''
        for s in stages:
            if stage_totals[s]:
                ns = int((mask & (hypno == s)).sum())
                cells += (f'<td style="padding:3px 10px;text-align:right;">{_pct(ns, stage_totals[s])}'
                          f'<small style="color:#888;"> ({ns})</small></td>')
            else:
                cells += '<td style="padding:3px 10px;text-align:right;">—</td>'
        return (f'<tr style="{weight}"><td style="padding:3px 10px;">{label}</td>'
                f'<td style="padding:3px 10px;text-align:right;">{n}</td>'
                f'<td style="padding:3px 10px;text-align:right;">{_pct(n, n_epochs)}</td>{cells}</tr>')
    any_mask = np.zeros(n_epochs, dtype=bool)
    for t in sorted(event_type_masks):
        any_mask |= np.asarray(event_type_masks[t], dtype=bool)
    head = ''.join(f'<th style="padding:3px 10px;text-align:right;">{s}</th>' for s in stages)
    body = ''.join(_row(t, event_type_masks[t]) for t in sorted(event_type_masks))
    body += _row('(any selected)', any_mask, bold=True)
    return ('<p style="font-size:.85em;color:#666;">Onset-only: each event flags the epoch '
            'containing its Start; annotated durations are ignored.</p>'
            '<table style="border-collapse:collapse;font-size:.9em;">'
            '<tr><th style="padding:3px 10px;text-align:left;">Event type</th>'
            '<th style="padding:3px 10px;text-align:right;">n epochs</th>'
            '<th style="padding:3px 10px;text-align:right;">%</th>'
            f'{head}</tr>{body}</table>')


def build_event_epoch_flags(event_type_masks, hypno_epochs, file_id):
    """Per-epoch event-type flag log: one row per epoch, one 0/1 column evt_<canonical> per selected
    event type (from event_type_masks, each a (n_epochs,) bool mask). The per-(epoch, channel)
    _epoch_channel_rejection.tsv only keeps the OR'd flag_event, so which event TYPE flagged which epoch
    is lost there; this sidecar preserves it so tool 8bis can re-select a subset of event types without
    reloading the raw events. Written only when event flagging ran on a file that has events."""
    n_epochs = len(hypno_epochs)
    cols = {
        'file_id':   file_id,
        'epoch_idx': np.arange(n_epochs),
        'stage':     np.asarray(hypno_epochs),
    }
    for t in sorted(event_type_masks):
        cols['evt_' + str(t)] = np.asarray(event_type_masks[t], dtype=bool).astype(int)
    return pd.DataFrame(cols)


def build_event_counts(events_df, event_remap, event_type_masks=None):
    """Per-canonical-type raw event count for one file: one row per canonical type present in events_df,
    n_events = number of scored events of that category (events_df rows whose event_remap[Name] maps to
    it), plus n_epochs_flagged from event_type_masks when available (reference only). Counts ALL events
    in the file, independent of which types were selected for flagging, so a downstream threshold can key
    on any category. Drives tool 8bis's per-event participant-exclusion threshold (raw count, not epochs)."""
    counts = {}
    if events_df is not None:
        for _, ev in events_df.iterrows():
            canonical = event_remap.get(str(ev['Name']))
            if canonical is None:
                continue
            counts[canonical] = counts.get(canonical, 0) + 1
    rows = []
    for t in sorted(counts):
        if event_type_masks and t in event_type_masks:
            n_flag = int(np.asarray(event_type_masks[t], dtype=bool).sum())
        else:
            n_flag = np.nan   # type present in the file but not selected for flagging -> no mask
        rows.append({'event_type': t, 'n_events': counts[t], 'n_epochs_flagged': n_flag})
    return pd.DataFrame(rows, columns=['event_type', 'n_events', 'n_epochs_flagged'])


def build_event_onsets(events_df, event_remap):
    """Per-event onset table for one file: one row per scored event that maps to a canonical type, with
    its onset (seconds from recording start) and annotated duration. Persisted as {file_id}_event_onsets.tsv
    so tool 8 can draw event-onset markers on the per-epoch montage WITHOUT reloading the raw recording
    (it otherwise has only per-epoch boolean flags). Counts ALL mapped events (not just the selected
    types), like build_event_counts. Returns a DataFrame [type, onset_s, duration_s] (empty if none)."""
    rows = []
    if events_df is not None:
        for _, ev in events_df.iterrows():
            canonical = event_remap.get(str(ev['Name']))
            if canonical is None:
                continue
            try:
                onset = float(ev['Start'])
            except Exception:
                continue
            try:
                dur = float(ev['Duration'])
            except Exception:
                dur = float('nan')
            rows.append({'type': canonical, 'onset_s': onset, 'duration_s': dur})
    return pd.DataFrame(rows, columns=['type', 'onset_s', 'duration_s'])


# ---------------------------------------------------------------------------------------------------------
# One recording's preprocessing and flagging (lifted from run_preprocessing of 7_preprocessing_voila: the
# lines are tool 7's own; its progress-bar updates became the optional `status` / `step_done` callbacks, its
# run-log prints the `warn` callback, and a fatal step raises StepFailed instead of `continue`)
# ---------------------------------------------------------------------------------------------------------
class StepFailed(Exception):
    """A fatal step of tool 7 for one recording: `reason` is its preprocessing_failed.tsv entry, `message`
    its run-log line. The recording is skipped, the run goes on."""

    def __init__(self, reason, message):
        super().__init__(message)
        self.reason = reason
        self.message = message


def select_montage_channels(raw, selected_channels, deselected_channels):
    """Tool 7's channel selection on a raw whose channels already carry their harmonised names: the selected
    channels are kept, the deselected ones kept too but marked bad (out of the reference and of the flagging),
    any other channel dropped. Returns the selected channels present (empty: nothing to process, `raw`
    untouched). These are lines of tool 7's load block, which the Curry generator rewrites: tool 7 keeps its own
    copy until it reads through the shared libraries (plan phase 8); tool 0 calls this one."""
    present = [ch for ch in selected_channels if ch in raw.ch_names]
    if not present:
        return present
    bad_present = [ch for ch in deselected_channels if ch in raw.ch_names]
    to_drop = [ch for ch in raw.ch_names if ch not in present and ch not in bad_present]
    if to_drop:
        raw.drop_channels(to_drop)
    raw.info['bads'] = bad_present
    return present


def preprocess_raw(raw, file_id, sub_config, do_resample, target_freq, do_notch, notch_freq_val, do_filter,
                   l_freq_val, h_freq_val, warn=print, status=None, step_done=None):
    """Tool 7's signal preparation, in place on the loaded `raw` (the montage channels, the deselected ones
    marked bad): resampling, re-reference (tool 2's `ref_channels`; the bad channels get the same reference
    subtracted), notch, band-pass. `status(text)` before a step, `step_done(name)` after it ('resample',
    'notch', 'filter'). Resampling, notch and band-pass are fatal (StepFailed), the re-reference is not."""
    # [C] Resampling — fatal
    if do_resample and target_freq is not None:
        if status:
            status('resampling…')
        try:
            raw.resample(target_freq, npad='auto', verbose=False)
        except Exception as e:
            raise StepFailed(f'resampling: {e}', f'[{file_id}] Error during resampling: {e}')
        if step_done:
            step_done('resample')

    # [D] Re-référencement — non-fatal (warn + continuer sans re-ref)
    ref_channels = sub_config.get('ref_channels', [])
    try:
        bads = list(raw.info['bads'])
        ref_signal = None   # the reference subtracted from the good channels (to apply to the bads)
        if ref_channels == 'average':
            # MNE averages the GOOD channels only: a dead channel never enters the reference
            if bads:
                ref_signal = raw.get_data(picks=[c for c in raw.ch_names if c not in bads]).mean(axis=0)
            raw.set_eeg_reference(ref_channels='average', verbose=False)
        elif isinstance(ref_channels, list) and ref_channels:
            # a deselected reference channel is not used (as before, when it was dropped)
            ref_present = [ch for ch in ref_channels if ch in raw.ch_names and ch not in bads]
            if ref_present:
                if bads:
                    ref_signal = raw.get_data(picks=ref_present).mean(axis=0)
                raw.set_eeg_reference(ref_channels=ref_present, verbose=False)
                raw.drop_channels(ref_present)
        # MNE leaves the bad channels in their original reference: give them the same one, so
        # the features computed on them downstream are comparable with the good channels'.
        if bads and ref_signal is not None:
            raw.apply_function(lambda x: x - ref_signal, picks=bads, channel_wise=True)
    except Exception as e:
        warn(f'[{file_id}] ⚠ Error during re-referencing: {e}, skipped.')

    # [D2] Notch filter — fatal
    if do_notch:
        if status:
            status('notch filtering…')
        try:
            raw.notch_filter(freqs=notch_freq_val, verbose=False)
        except Exception as e:
            raise StepFailed(f'notch filtering: {e}', f'[{file_id}] Error during notch filtering: {e}')
        if step_done:
            step_done('notch')

    # [E] Filtrage — fatal
    if do_filter:
        if status:
            status('band-pass filtering…')
        try:
            raw.filter(
                l_freq=l_freq_val, h_freq=h_freq_val,
                l_trans_bandwidth='auto', h_trans_bandwidth='auto',
                filter_length='auto', method='fir', phase='zero-double',
                fir_window='hamming', fir_design='firwin', verbose=False
            )
        except Exception as e:
            raise StepFailed(f'filtering: {e}', f'[{file_id}] Error during filtering: {e}')
        if step_done:
            step_done('filter')
    return raw


def prepare_hypnogram(hypno_path, raw, sf, epoch_factor, custom_stages, file_id, warn=print):
    """Tool 7's hypnogram: read (fatal on failure), then hypnogram_epochs.
    Returns (hypno_sub, stage_mapping, hypno_vec)."""
    try:
        expert_hypno = np.loadtxt(str(hypno_path), dtype=str).astype('<U10')
    except Exception as e:
        raise StepFailed(f'hypno loading: {e}', f'[{file_id}] Error loading hypnogram: {e}')
    return hypnogram_epochs(expert_hypno, raw, sf, epoch_factor, custom_stages, file_id, warn)


def hypnogram_epochs(expert_hypno, raw, sf, epoch_factor, custom_stages, file_id, warn=print):
    """Tool 7's rules on the scored labels (one per 30 s epoch): checked against the complete 30 s epochs of
    the signal (one extra score is dropped with a note, any other mismatch is fatal), each 30 s score repeated
    over its sub-epochs, and every label given its MNE event code (W 0 ... R 4, then the custom stages, then any
    other label). Tool 0 passes the hypnogram in use (rescorings included). Returns (hypno_sub, stage_mapping,
    hypno_vec)."""
    expected_epochs = int(np.floor(raw.n_times / sf / 30))
    n_hyp = len(expert_hypno)
    n_over = n_hyp - expected_epochs
    if 0 < n_over <= 1:
        # Recording ends mid-epoch: the last epoch is partial (< 30 s) but was still
        # scored. Drop the overhanging score and epoch only the complete 30 s epochs.
        last_sec = raw.n_times / sf - expected_epochs * 30
        expert_hypno = expert_hypno[:expected_epochs]
        warn(f'[{file_id}] Last epoch {last_sec:.1f}s (< 30s) skipped; hypnogram '
             f'trimmed {n_hyp} -> {expected_epochs} to match complete 30s epochs.')
    elif n_over != 0:
        raise StepFailed('hypno/EEG length mismatch', f'[{file_id}] Hypnogram/EEG length mismatch ({n_hyp} '
                                                      f'vs {expected_epochs} epochs): skipped.')

    # Re-cut the 30 s scored hypnogram into sub-epochs: each sub-epoch inherits the sleep
    # stage of its parent 30 s epoch (epoch_factor == 1 for the classic 30 s case → no change).
    hypno_sub = np.repeat(expert_hypno, epoch_factor)

    stage_mapping = {'W': 0, 'N1': 1, 'N2': 2, 'N3': 3, 'R': 4}
    # Every other label present (declared custom stages first, then e.g. MT) gets its own code from 5:
    # MNE keeps only the epochs whose code is in event_id when it READS a .fif back, so a stage left
    # without a code would make the saved file unreadable by tools 8 and 9.
    _extra = [s for s in custom_stages if s in set(hypno_sub)]
    _extra += sorted(set(hypno_sub) - set(stage_mapping) - set(_extra))
    for _code, _stage in enumerate(_extra, start=len(stage_mapping)):
        stage_mapping[str(_stage)] = _code
    hypno_vec = np.full(len(hypno_sub), -1, dtype=float)
    for stage, value in stage_mapping.items():
        hypno_vec[hypno_sub == stage] = value
    return hypno_sub, stage_mapping, hypno_vec


def make_epochs(raw, epoch_sec, hypno_sub, hypno_vec, stage_mapping, file_id):
    """Tool 7's epochs: fixed length from t = 0, as many as the (expanded) hypnogram, each epoch's event code
    = its stage. Fatal on failure or when no epoch is left (StepFailed). Returns (epochs, n_epochs); the
    caller frees `raw` right after (memory)."""
    # [F] Epoching — fatal
    try:
        epochs = mne.make_fixed_length_epochs(raw, duration=epoch_sec, preload=True, verbose=False)
        epochs = epochs[:len(hypno_sub)]
        epochs.events[:, 2] = hypno_vec
        epochs.event_id     = stage_mapping
        n_epochs = len(epochs)
    except Exception as e:
        raise StepFailed(f'epoching: {e}', f'[{file_id}] Error creating epochs: {e}')

    if n_epochs == 0:
        raise StepFailed('no epochs', f'[{file_id}] No epochs: skipped.')
    return epochs, n_epochs


def set_bads_aside(epochs):
    """Tool 7's [F2]: the deselected (bad) channels leave the epochs here and come back only at save time,
    so the PSD, the flagging, the tables and the report work on the flagged channels only.
    Returns (bad_names, bad_epochs or None, all_channel_order or None)."""
    # [F2] The deselected (bad) channels leave the epochs here and come back only at save time:
    # everything below (PSD, flagging, tables, report) works on the channels actually flagged,
    # exactly as when they were dropped. Only the bad channels' data is copied aside.
    bad_names  = list(epochs.info['bads'])
    bad_epochs = None
    all_channel_order = None
    if bad_names:
        all_channel_order = list(epochs.ch_names)
        bad_epochs = mne.EpochsArray(
            epochs.get_data(picks=bad_names),
            mne.pick_info(epochs.info, mne.pick_channels(epochs.ch_names, bad_names, ordered=True)),
            events=epochs.events, tmin=epochs.tmin, event_id=epochs.event_id, verbose=False)
        epochs.drop_channels(bad_names)
    return bad_names, bad_epochs, all_channel_order


def epoch_psd_for_1f(epochs, sf, epoch_sec, do_1f, fit_fmax_val, ss_enabled_val, ss_median_val,
                     ss_median_span_val, ss_span_val, file_id, warn=print):
    """Tool 7's Welch PSD of every (epoch, channel), µV²/Hz, smoothed before the 1/f fit when asked (median
    then LOWESS, oscip order). Only when the 1/f method is on; non-fatal (on failure the 1/f method is off
    for this recording). Returns (psds_data_uV2, psd_freqs), both None when not computed."""
    # PSD — only when 1/f flagging is enabled (skipping it is the main speed-up on short
    # epochs). Non-fatal: on failure the 1/f method is disabled for this file.
    psds_data_uV2 = None
    psd_freqs     = None
    if do_1f:
        try:
            n_per_seg = int(min(4, epoch_sec) * sf)   # 4 s Welch window, capped to the epoch length
            n_overlap = int(n_per_seg / 2)
            fmax_psd  = min(fit_fmax_val, sf / 2 - 0.5)  # PSD spans up to the 1/f fit max
            psds_obj  = epochs.compute_psd(
                method='welch', fmin=0.5, fmax=fmax_psd,
                n_fft=n_per_seg, n_overlap=n_overlap, n_per_seg=n_per_seg,
                window='hann', verbose=False
            )
            psds_data_uV2 = psds_obj.get_data() * 1e12  # µV²/Hz
            psd_freqs     = psds_obj.freqs
            # PSD smoothing: smooth the PSD before the 1/f fit (oscip parity) — tames the noisy
            # few-segment Welch PSD on short epochs so specparam still fits it. Two stages,
            # in oscip order: a running median (kills narrow spikes) then the LOWESS mean.
            if ss_enabled_val and psds_data_uV2 is not None:
                if ss_median_val:
                    psds_data_uV2 = smooth_psd_median(psds_data_uV2, psd_freqs, span_hz=ss_median_span_val)
                psds_data_uV2 = smooth_psd_lowess(psds_data_uV2, psd_freqs, span_hz=ss_span_val)
        except Exception as e:
            warn(f'[{file_id}] PSD warning: {e}, 1/f flagging disabled.')
            psds_data_uV2 = None
            psd_freqs     = None
    return psds_data_uV2, psd_freqs


def event_flags(events_df, file_remap, event_types, n_epochs, epoch_sec):
    """Tool 7's event method for one recording: the epochs containing the onset of a selected event type
    (all selected types, and each type alone), the raw count per canonical type (tool 8bis's thresholds) and
    the onset table (tool 8's montage). `file_remap`: event_remap.json with the manual annotations mapped to
    themselves. Returns (event_epoch_mask, event_type_masks, event_counts_df, event_onsets_df)."""
    event_epoch_mask = compute_event_epoch_mask(
        events_df, file_remap, event_types, n_epochs, epoch_sec
    )
    event_type_masks = {
        t: compute_event_epoch_mask(events_df, file_remap, [t], n_epochs, epoch_sec)
        for t in sorted(event_types)
    }
    # Raw per-type event count for the whole file (all canonical types, not just the
    # selected ones) — persisted so tool 8bis can exclude a participant that has more
    # events of a category than a user threshold.
    event_counts_df = build_event_counts(events_df, file_remap, event_type_masks)
    # Per-event onset table (all mapped events) -> {file_id}_event_onsets.tsv so tool 8
    # can mark event onsets on the montage without reloading the raw recording.
    event_onsets_df = build_event_onsets(events_df, file_remap)
    return event_epoch_mask, event_type_masks, event_counts_df, event_onsets_df


def epoch_metadata(mask_dict, hypno_epochs, event_type_masks, n_epochs):
    """Tool 7's per-epoch metadata of the .fif: stage, reject_flag (any method on any channel: the naive
    default, not the decision), reject_method (the method, or 'multiple'), one flag_<method> per method and,
    when event flagging ran, one evt_<type> per event type. Returns (reject_epoch, metadata_df)."""
    reject_epoch = np.zeros(n_epochs, dtype=bool)
    for name in mask_dict:
        reject_epoch |= mask_dict[name].any(axis=1)

    method_priority = [m for m in ['amplitude', 'gradient', 'flat', '1f_r2', '1f_error', 'event']
                       if m in mask_dict]
    reject_method   = np.full(n_epochs, 'none', dtype=object)
    for name in reversed(method_priority):
        reject_method[mask_dict[name].any(axis=1)] = name
    n_methods = sum(mask_dict[name].any(axis=1).astype(int) for name in method_priority)
    reject_method[n_methods > 1] = 'multiple'

    meta_rows = []
    for ei in range(n_epochs):
        row = {
            'epoch_idx':     ei,
            'stage':         hypno_epochs[ei],
            'reject_flag':   bool(reject_epoch[ei]),
            'reject_method': reject_method[ei],
        }
        for name in method_priority:
            col = 'flag_' + name.replace('/', '').replace(' ', '')
            row[col] = bool(mask_dict[name][ei].any())
        # Per-event-type epoch flags (only when event flagging ran): lets tool 8's per-epoch
        # inspector show WHICH event type flagged the epoch (flag_event alone is the OR).
        if event_type_masks:
            for t in sorted(event_type_masks):
                row['evt_' + str(t)] = bool(np.asarray(event_type_masks[t])[ei])
        meta_rows.append(row)
    metadata_df = pd.DataFrame(meta_rows)
    return reject_epoch, metadata_df
