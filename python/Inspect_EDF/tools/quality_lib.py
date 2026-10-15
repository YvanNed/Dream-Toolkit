"""Signal-quality computation of tool 6 (6_quality_overview_voila), shared with tool 0.

Moved verbatim from tool 6 (its shared-function cell, and the computation lines of its run loop): tool 6 and
tool 0's step 5 call the same code, so they give the same numbers by construction. Edit here, never in a
notebook. The Curry twin of tool 6 does not import this file: tools_curry/_make_tool6_curry.py copies it into
the twin and strips the EDF physical-range metric (bounds_pct) from that copy, so re-run the generator after
any edit here.

Format-agnostic: everything works on an already-loaded mne Raw (µV = get_data() * 1e6). Reading the recording
(EDF here, Curry in the twin) stays in the notebooks.
"""
import json
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, savgol_filter
from scipy.stats import kurtosis as scipy_kurtosis, skew as scipy_skew


# ---------------------------------------------------------------------------------------------------------
# Metrics, flags and figures (tool 6's shared-function cell)
# ---------------------------------------------------------------------------------------------------------
def odd(x):
    """Return nearest odd integer to x (required by savgol_filter)."""
    n = int(round(x))
    if n % 2 != 0:
        return n
    return n - 1 if abs(x - (n - 1)) <= abs((n + 1) - x) else n + 1


# Robust percentiles bounding the histogram used for the Savitzky-Golay curve + peak detection
# (see compute_signal_metrics); keeps the amplitude range from being stretched by rare extremes.
HIST_CLIP_PCT = (0.5, 99.5)


def compute_signal_metrics(sig_uV, phys_min_uV, phys_max_uV):
    """Compute all quality metrics for one EEG channel (unfiltered signal, in uV)."""
    cur_std = float(np.std(sig_uV))
    u = np.unique(sig_uV)

    # Flat signal: fraction of consecutive sample pairs with no change.
    # Tolerance = 2x ADC resolution to handle rounding, minimum 0.06 uV.
    if len(u) > 1:
        res = float(np.min(np.abs(np.diff(u))))
        flat_atol = max(2 * res, 0.06)
    else:
        flat_atol = 0.06
    flat_pct = float(np.isclose(np.diff(sig_uV), 0, atol=flat_atol, rtol=0.0).mean()) * 100

    # EDF physical bounds: fraction of samples at/near the header declared limits.
    # Detects saturation at the EDF dynamic range boundary (hard clipping).
    if not (np.isnan(phys_min_uV) or np.isnan(phys_max_uV)):
        bounds_pct = float(
            ((sig_uV <= phys_min_uV + 0.5) | (sig_uV >= phys_max_uV - 0.5)).mean()
        ) * 100
    else:
        bounds_pct = 0.0

    # Histogram + Savitzky-Golay smoothing + peak detection.
    # The histogram feeding the SG curve + peaks is bounded to robust percentiles (HIST_CLIP_PCT).
    # Without it, rare extreme samples (DC drift/artifacts on DC-coupled data, occasional spikes on
    # EDF) stretch the amplitude range, inflating the Scott bin count and the 0.02*n_bins smoothing
    # window -> the curve over-smooths into a dome/flat line. Clipping concentrates the bins on the
    # bulk of the signal so the curve tracks the distribution and resolves real modes.
    # SG window = 2% of bin count (at least 5); peak prominence = 25% of the histogram maximum.
    # hist_extreme_pct stays on the FULL-range histogram (it is a clipping proxy on the outer bins).
    if cur_std < 1e-10 or len(u) < 3:
        histo_full, edges_full = np.histogram(sig_uV, bins=10)
    else:
        histo_full, edges_full = np.histogram(sig_uV, bins='scott')

    if cur_std < 1e-10 or len(u) < 3:
        histo, edges = histo_full, edges_full          # degenerate/flat channel: no clip
        sg_window = 5
    else:
        lo, hi = np.percentile(sig_uV, HIST_CLIP_PCT)
        if hi > lo:
            histo, edges = np.histogram(sig_uV, bins='scott', range=(float(lo), float(hi)))
        else:
            histo, edges = histo_full, edges_full      # degenerate robust range: fall back to full
        sg_window = max(5, odd(0.02 * len(histo)))
    bin_centers = (edges[:-1] + edges[1:]) / 2
    y_smooth = savgol_filter(histo.astype(float), window_length=sg_window, polyorder=3)
    prom_thresh = 0.25 * float(np.max(y_smooth)) if np.max(y_smooth) > 0 else 0.0
    peaks, _ = find_peaks(y_smooth, prominence=prom_thresh)

    # Extreme histogram bins: fraction of samples in the outermost bins.
    # Proxy for in-range clipping (signal saturating within the declared EDF physical range).
    total = float(np.sum(histo_full))
    hist_extreme_pct = float((histo_full[0] + histo_full[-1]) / total * 100) if total > 0 else 0.0

    # Percentile-based amplitude metrics: robust characterization of the signal tail.
    # p99 / p99.9 of |amplitude| capture sustained artifacts (movement bursts, electrode pops)
    # without being dominated by a single sample like the max would be.
    abs_uV = np.abs(sig_uV)

    # Distribution statistics.
    # Fischer kurtosis: 0 = normal; >0 = heavy tails/spikes; <0 = flat/platykurtic.
    return {
        'mean_uV': float(np.mean(sig_uV)),
        'std_uV': cur_std,
        'kurtosis': float(scipy_kurtosis(sig_uV, fisher=True)),
        'skewness': float(scipy_skew(sig_uV)),
        'p99_abs_uV': float(np.percentile(abs_uV, 99.0)),
        'p999_abs_uV': float(np.percentile(abs_uV, 99.9)),
        'flat_pct': flat_pct,
        'bounds_pct': bounds_pct,
        'hist_extreme_pct': hist_extreme_pct,
        'n_peaks': int(len(peaks)),
        'histo': histo,
        'edges': edges,
        'bin_centers': bin_centers,
        'y_smooth': y_smooth,
        'peaks': peaks,
    }


def flag_channel(metrics, thresholds):
    """Return list of flag reasons for a channel, or empty list if within all thresholds."""
    reasons = []
    if metrics['flat_pct'] > thresholds['flat_pct']:
        reasons.append(f"flat_pct={metrics['flat_pct']:.1f}% > {thresholds['flat_pct']}%")
    if metrics['bounds_pct'] > thresholds['bounds_pct']:
        reasons.append(f"bounds_pct={metrics['bounds_pct']:.2f}% > {thresholds['bounds_pct']}%")
    if metrics['n_peaks'] >= thresholds['n_peaks']:
        shape = 'bimodal' if metrics['n_peaks'] == 2 else 'multimodal'
        reasons.append(f"n_peaks={metrics['n_peaks']} ({shape})")
    # Low std detects dead or near-dead channels (normal EEG std >> 5 µV).
    # Kurtosis was removed: EEG PSG has physiologically high kurtosis (spindles, K-complexes)
    # causing too many false positives (baseline Fp1 kurtosis ~466).
    if metrics['std_uV'] < thresholds['std_low']:
        reasons.append(f"std={metrics['std_uV']:.2f} µV < {thresholds['std_low']} µV (dead/flat channel)")
    # Extreme histogram bins proxy for in-range clipping (signal saturating within
    # the declared EDF physical range, undetectable via bounds_pct alone).
    if metrics['hist_extreme_pct'] > thresholds['hist_extreme_pct']:
        reasons.append(f"hist_extreme={metrics['hist_extreme_pct']:.2f}% > {thresholds['hist_extreme_pct']}%")
    return reasons


def plot_histogram_figure(metrics, ch_name, y_max=None, x_lim=None):
    """Return a matplotlib figure with the signal amplitude distribution."""
    edges = metrics['edges']
    histo = metrics['histo']
    bin_centers = metrics['bin_centers']
    y_smooth = metrics['y_smooth']
    peaks = metrics['peaks']
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.bar(edges[:-1], histo, width=edges[1] - edges[0], align='edge',
           color='#d4e8f5', edgecolor='#7ab4d0')
    ax.plot(bin_centers, y_smooth, color='#6a51a3', linewidth=2,
            label='Savitzky-Golay smooth')
    ax.plot(bin_centers[peaks], y_smooth[peaks], 'r+', markersize=10, markeredgewidth=2,
            label=f'{len(peaks)} peak(s) detected')
    ax.set_xlabel('Amplitude (μV)', fontsize=11)
    ax.set_ylabel('Counts', fontsize=11)
    if y_max is not None:
        ax.set_ylim(0, y_max * 1.05)
    if x_lim is not None:
        ax.set_xlim(-x_lim, x_lim)
    for sp in ['right', 'top']:
        ax.spines[sp].set_visible(False)
    ax.legend(fontsize=9)
    plt.tight_layout()
    return fig


def plot_psd_figure(sig_uV, sf):
    """Global PSD (Welch, 4 s windows, 50 % overlap) on the unfiltered signal."""
    # Frequency bands for visual annotation (non-overlapping, covers 0.5–30 Hz).
    BANDS = [
        (0.5,  4,  'δ',  '#cce0f5'),
        (4,    8,  'θ',  '#ccf0d8'),
        (8,   12,  'α',  '#fff2cc'),
        (12,  16,  'σ',  '#f5d0f0'),
        (16,  30,  'β',  "#ffbaba"),
    ]
    # n_overlap derived from n_per_seg (not independently from sf) to guarantee
    # n_overlap < n_per_seg regardless of floating-point rounding in sf.
    # n_fft passed explicitly to prevent MNE from defaulting it to 256 and silently
    # capping n_per_seg, which would cause n_overlap > n_per_seg.
    n_per_seg = int(4 * sf)
    n_overlap = n_per_seg // 2
    psds, freqs = mne.time_frequency.psd_array_welch(
        sig_uV[np.newaxis, :], sfreq=sf,
        fmin=0.5, fmax=min(100.0, sf / 2 - 0.5),
        n_per_seg=n_per_seg, n_overlap=n_overlap, n_fft=n_per_seg,
        verbose=False
    )
    psd_db = 10 * np.log10(np.maximum(psds[0], 1e-10))
    fig, ax = plt.subplots(figsize=(7, 3))
    for flo, fhi, label, color in BANDS:
        fhi_clipped = min(fhi, float(freqs[-1]))
        ax.axvspan(flo, fhi_clipped, alpha=0.20, color=color, zorder=0, linewidth=0)
        ax.text((flo + fhi_clipped) / 2, 0.97, label, ha='center', va='top',
                fontsize=8, color='#555', transform=ax.get_xaxis_transform())
    ax.plot(freqs, psd_db, linewidth=1.0, color='#2a6099', zorder=2)
    ax.axvline(50, color='#c0392b', linewidth=0.8, linestyle='--', alpha=0.6, zorder=3)
    ax.text(50, 0.97, '50 Hz', ha='center', va='top', fontsize=7.5,
            color='#c0392b', transform=ax.get_xaxis_transform())
    ax.set_xlabel('Frequency (Hz)', fontsize=11)
    ax.set_ylabel('Power (dB/Hz)', fontsize=11)
    ax.set_xlim(0.5, float(freqs[-1]))
    band_ticks = [t for t in [0.5, 4, 8, 12, 16, 30, 50, 70, 90] if t <= float(freqs[-1])]
    ax.set_xticks(band_ticks)
    ax.set_xticklabels([str(t) for t in band_ticks], fontsize=9)
    for sp in ['right', 'top']:
        ax.spines[sp].set_visible(False)
    plt.tight_layout()
    return fig


def plot_timeseries_figure(sig_uV, sf, y_lim):
    """Raw time series (downsampled to ~10 Hz for display), fixed Y-axis."""
    step = max(1, int(sf / 10))
    t_h = np.arange(0, len(sig_uV), step) / sf / 3600
    fig, ax = plt.subplots(figsize=(10, 2))
    ax.plot(t_h, sig_uV[::step], linewidth=0.3, color='#2a6099', rasterized=True)
    ax.set_xlim(0, len(sig_uV) / sf / 3600)
    ax.set_ylim(-y_lim, y_lim)
    ax.set_xlabel('Time (h)', fontsize=11)
    ax.set_ylabel('µV', fontsize=11)
    ax.axhline(0, color='#bbb', linewidth=0.5, linestyle='--')
    for sp in ['right', 'top']:
        ax.spines[sp].set_visible(False)
    plt.tight_layout()
    return fig


# --- Overview overlays: all channels superimposed, one colour per channel -----------------
# Feed the report's first section: a raw.plot_psd-style multi-channel view (time series,
# PSD, amplitude distribution) plus an across-channels metric summary. Format-agnostic
# (no EDF physical-bounds dependency), so the Curry twin reuses them unchanged.
def _channel_colors(n):
    """n visually distinct colours (turbo colormap)."""
    cmap = plt.get_cmap('turbo')
    return [cmap(i / max(1, n - 1)) for i in range(n)] if n > 1 else [cmap(0.5)]


def _overlay_legend(ax):
    """Compact per-channel legend placed to the right of the axes."""
    handles, labels = ax.get_legend_handles_labels()
    ncol = 1 if len(labels) <= 16 else 2
    ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.01, 1.0),
              ncol=ncol, fontsize=6.5, frameon=False, handlelength=1.2, columnspacing=0.8)


def plot_overlay_timeseries(raw, sf, y_lim=None):
    """Butterfly view: every channel's time series overlaid on one uV axis.
    Reads one channel at a time from `raw` (a cheap slice of the preloaded data) so the
    whole multi-channel recording is never held as a second float64 array."""
    ch_names = raw.ch_names
    n_times = raw.n_times
    step = max(1, int(sf / 10))   # ~10 Hz display resolution
    t_h = np.arange(0, n_times, step) / sf / 3600
    colors = _channel_colors(len(ch_names))
    _alpha = max(0.3, min(0.8, 10.0 / max(1, len(ch_names))))  # more channels -> more transparent
    fig, ax = plt.subplots(figsize=(10, 3.4))
    for i, ch in enumerate(ch_names):
        sig = raw.get_data(picks=[i])[0] * 1e6
        ax.plot(t_h, sig[::step], lw=0.3, color=colors[i], alpha=_alpha, label=ch, rasterized=True)
        del sig
    ax.set_xlim(0, n_times / sf / 3600)
    if y_lim is not None:
        ax.set_ylim(-y_lim, y_lim)
    ax.axhline(0, color='#bbb', lw=0.5, ls='--')
    ax.set_xlabel('Time (h)', fontsize=11)
    ax.set_ylabel('µV', fontsize=11)
    for sp in ['right', 'top']:
        ax.spines[sp].set_visible(False)
    _overlay_legend(ax)
    fig.subplots_adjust(right=0.82)
    return fig


def plot_overlay_distribution(raw, x_lim=None):
    """Amplitude distribution (density) of every channel overlaid as outline curves
    (one channel at a time)."""
    ch_names = raw.ch_names
    colors = _channel_colors(len(ch_names))
    _alpha = max(0.3, min(0.8, 10.0 / max(1, len(ch_names))))  # more channels -> more transparent
    fig, ax = plt.subplots(figsize=(9, 3.6))
    for i, ch in enumerate(ch_names):
        sig = raw.get_data(picks=[i])[0] * 1e6
        bins = 'scott' if float(np.std(sig)) > 1e-10 else 10
        histo, edges = np.histogram(sig, bins=bins)
        width = edges[1] - edges[0]
        total = histo.sum()
        dens = histo / (total * width) if (total > 0 and width > 0) else histo.astype(float)
        centers = (edges[:-1] + edges[1:]) / 2
        ax.plot(centers, dens, lw=1.0, color=colors[i], alpha=_alpha, label=ch)
        del sig
    if x_lim is not None:
        ax.set_xlim(-x_lim, x_lim)
    ax.set_xlabel('Amplitude (µV)', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    for sp in ['right', 'top']:
        ax.spines[sp].set_visible(False)
    _overlay_legend(ax)
    fig.subplots_adjust(right=0.82)
    return fig


def plot_overlay_psd(raw, sf):
    """Welch PSD (4 s windows) of every channel overlaid on one axis (one channel at a time)."""
    ch_names = raw.ch_names
    n_per_seg = int(4 * sf)
    n_overlap = n_per_seg // 2
    colors = _channel_colors(len(ch_names))
    _alpha = max(0.3, min(0.8, 10.0 / max(1, len(ch_names))))  # more channels -> more transparent
    fig, ax = plt.subplots(figsize=(9, 3.6))
    freqs = None
    for i, ch in enumerate(ch_names):
        sig = raw.get_data(picks=[i])[0] * 1e6
        psds, freqs = mne.time_frequency.psd_array_welch(
            sig[np.newaxis, :], sfreq=sf, fmin=0.5, fmax=min(100.0, sf / 2 - 0.5),
            n_per_seg=n_per_seg, n_overlap=n_overlap, n_fft=n_per_seg, verbose=False)
        psd_db = 10 * np.log10(np.maximum(psds[0], 1e-10))
        ax.plot(freqs, psd_db, lw=0.9, color=colors[i], alpha=_alpha, label=ch)
        del sig, psds
    if freqs is not None:
        ax.set_xlim(0.5, float(freqs[-1]))
    ax.set_xlabel('Frequency (Hz)', fontsize=11)
    ax.set_ylabel('Power (dB/Hz)', fontsize=11)
    for sp in ['right', 'top']:
        ax.spines[sp].set_visible(False)
    _overlay_legend(ax)
    fig.subplots_adjust(right=0.82)
    return fig


 

def build_avg_metrics_html(ch_metrics, ch_names):
    """Table of per-channel metrics summarised across channels (mean/median/min/max).
    Each metric is guarded on key presence, so bounds_pct (EDF-only) is skipped for
    float-based recordings without any code change."""
    METRICS = [
        ('std_uV', 'Std dev (\u00b5V)'), ('flat_pct', 'Flat signal (%)'),
        ('bounds_pct', 'At EDF bounds (%)'), ('hist_extreme_pct', 'Extreme histogram (%)'),
        ('p99_abs_uV', 'p99 |amplitude| (\u00b5V)'), ('p999_abs_uV', 'p99.9 |amplitude| (\u00b5V)'),
        ('kurtosis', 'Kurtosis (Fischer)'), ('skewness', 'Skewness'), ('n_peaks', 'n_peaks'),
    ]
    body = ''
    for key, label in METRICS:
        vals = [ch_metrics[ch][key] for ch in ch_names if key in ch_metrics[ch]]
        if not vals:
            continue
        arr = np.asarray(vals, dtype=float)
        body += (f'<tr><td style="padding:4px 10px;">{label}</td>'
                 f'<td style="padding:4px 10px;text-align:right;">{arr.mean():.3g}</td>'
                 f'<td style="padding:4px 10px;text-align:right;">{np.median(arr):.3g}</td>'
                 f'<td style="padding:4px 10px;text-align:right;">{arr.min():.3g}</td>'
                 f'<td style="padding:4px 10px;text-align:right;">{arr.max():.3g}</td></tr>')
    return (
        f'<p style="margin:6px 0;">Summary across the {len(ch_names)} plotted channel(s).</p>'
        '<table style="border-collapse:collapse;font-size:0.9em;">'
        '<tr style="background:#eee;">'
        '<th style="padding:4px 10px;text-align:left;">Metric</th>'
        '<th style="padding:4px 10px;text-align:right;">Mean</th>'
        '<th style="padding:4px 10px;text-align:right;">Median</th>'
        '<th style="padding:4px 10px;text-align:right;">Min</th>'
        '<th style="padding:4px 10px;text-align:right;">Max</th></tr>'
        f'{body}</table>'
    )


# Per-metric interpretation shown as a third column in the metrics table.
METRIC_INTERPRETATIONS = {
    'Kurtosis (Fischer)': 'PSG EEG: typically ~100&ndash;300; compare across channels, a channel significantly lower than the others&nbsp;=&nbsp;suspicious',
    'Skewness': '~0&nbsp;=&nbsp;symmetric; |skew|&nbsp;&gt;&nbsp;1&nbsp;=&nbsp;notable asymmetry',
    'n_peaks (distribution)': '1&nbsp;=&nbsp;unimodal; &ge;2&nbsp;=&nbsp;irregular / suspect distribution (bimodal or multimodal: inspect the channel)',
    'p99 |amplitude|': '99% of samples are below this value&nbsp;; compare across channels&nbsp;',
    'p99.9 |amplitude|': '99.9% of samples are below this value&nbsp;; compare across channels&nbsp;',
}


# --- Custom (non-AASM) sleep stages -------------------------------------------
# Projects may intentionally keep labels outside the AASM set (W/N1/N2/N3/R). They are
# declared in config_param/custom_stages.json (written by 3_remap_hypno) and recognised
# here so the hypnospectrogram and per-stage tables include them.
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
    AASM keeps YASA heights (W top ... N3 bottom); custom stages stack below N3 (-1, -2, ...)
    in declaration order with auto-assigned colours. Returns (stage_y, stage_colors, ytick_pos, ytick_labels)."""
    stage_y = {'W': 4, 'R': 3, 'N1': 2, 'N2': 1, 'N3': 0}
    stage_colors = dict(BASE_STAGE_COLORS)
    for i, cs in enumerate(custom_stages):
        stage_y[cs] = -1 - i
        stage_colors[cs] = CUSTOM_STAGE_PALETTE[i % len(CUSTOM_STAGE_PALETTE)]
    ordered = sorted(stage_y.items(), key=lambda kv: kv[1])
    ytick_pos = [y for _, y in ordered]
    ytick_labels = ['REM' if s == 'R' else s for s, _ in ordered]
    return stage_y, stage_colors, ytick_pos, ytick_labels


def compute_hypnospectrogram_sxx(sig_uV, sf, fmin=0.5, fmax=40):
    """Multitaper spectrogram of one channel, in dB, restricted to [fmin, fmax].
    Split out of plot_hypnospectrogram so the same Sxx can be computed once and used twice:
    for the per-channel figure and for the across-channel median in the Overview section.
    Returns (f, t_hours, Sxx_dB)."""
    from yasa.plotting import spectrogram_lspopt

    # Multitaper spectrogram (mirrors yasa.plot_spectrogram internals).
    nperseg = int(30 * sf)
    f, t, Sxx = spectrogram_lspopt(sig_uV, sf, nperseg=nperseg, noverlap=0)
    with np.errstate(divide='ignore'):
        Sxx = 10 * np.log10(Sxx)
    good = np.logical_and(f >= fmin, f <= fmax)
    Sxx, f = Sxx[good, :], f[good]
    t = t / 3600  # seconds -> hours
    return f, t, Sxx


# --- Display conventions of every hypnogram / hypnospectrogram of the toolkit -------------------------------
REM_COLOR = '#c0392b'          # REM in red on every hypnogram, at the line's own width (the usual PSG convention)
SPECTROGRAM_CMAP = 'viridis'   # the palette of every hypnospectrogram (tool 6, tool 0's navigator and step 5)
DEAD_EPOCH_COLOR = '#d9d9d9'   # spectrogram columns without signal (flat / disconnected epochs)


def spectrogram_colour_scale(Sxx, trimperc=2.5):
    """Colour scale of a dB spectrogram (frequency x time columns), robust to flat / disconnected stretches:
    the dead columns are left out of the percentiles and drawn grey, the span is capped at 45 dB. Shared by
    every hypnospectrogram (tool 6's figures, tool 0's navigator). Returns (Sxx with the dead columns set to
    NaN, the Normalize, the colormap)."""
    from matplotlib.colors import Normalize

    # Flat / disconnected epochs (whole 30 s columns) have near-zero power -> extreme negative
    # dB; left in, they drag the percentile colour scale down and wash the plot to one colour,
    # so exclude them from vmin/vmax and render them grey ('no signal').
    Sxx = np.where(np.isfinite(Sxx), Sxx, np.nan)
    with np.errstate(invalid='ignore'):
        col_max = np.nanmax(np.where(np.isfinite(Sxx), Sxx, -np.inf), axis=0)
    col_max[~np.isfinite(col_max)] = np.nan
    ref = np.nanmedian(col_max)
    valid = np.isfinite(col_max) & (col_max > ref - 60)   # 60 dB below the median epoch peak
    pix = Sxx[:, valid]
    pix = pix[np.isfinite(pix)]
    if pix.size < 10:        # almost everything dead -> graceful fallback
        pix = Sxx[np.isfinite(Sxx)]
    vmin, vmax = np.percentile(pix, [0 + trimperc, 100 - trimperc])
    vmin = max(vmin, vmax - 45)   # cap colour span (dB): real EEG structure fits within ~35 dB, so a few degraded low-power epochs can't stretch the scale and wash the plot to one colour
    norm = Normalize(vmin=vmin, vmax=vmax)
    Sxx[:, ~valid] = np.nan   # dead epochs rendered as the colormap 'bad' colour
    spec_cmap = plt.get_cmap(SPECTROGRAM_CMAP).copy()
    spec_cmap.set_bad(DEAD_EPOCH_COLOR)
    return Sxx, norm, spec_cmap


def render_hypnospectrogram(f, t, Sxx, hypno_labels, custom_stages, trimperc=2.5):
    """Draw the hypnogram strip + spectrogram from an already-computed dB spectrogram.
    Split out of plot_hypnospectrogram so the Overview median spectrogram reuses the exact
    same rendering (colour scaling, dead-epoch greying) as the per-channel figures.
    Returns a Figure."""
    hypno_labels = np.asarray(hypno_labels)
    stage_y, stage_colors, ytick_pos, ytick_labels = custom_stage_style(custom_stages)
    Sxx, norm, spec_cmap = spectrogram_colour_scale(Sxx, trimperc)

    # Two rows (hypnogram strip, spectrogram) and a thin right column for the colour scale, used only beside
    # the spectrogram: both plots keep the same width (aligned in time) and the scale is the spectrogram's height.
    fig = plt.figure(figsize=(10, 4.2))
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 2], width_ratios=[1, 0.022], hspace=0.1, wspace=0.03)
    ax0 = fig.add_subplot(grid[0, 0])
    ax1 = fig.add_subplot(grid[1, 0])
    cax = fig.add_subplot(grid[1, 1])

    # Hypnogram strip: one y per 30s epoch; black line, REM red, each custom stage in its colour, all at the
    # line's own width (discreet).
    n_ep = len(hypno_labels)
    ep_hours = np.arange(n_ep + 1) * 30 / 3600
    y = np.array([stage_y.get(s, np.nan) for s in hypno_labels], dtype=float)
    y_step = np.append(y, y[-1] if n_ep else np.nan)
    ax0.step(ep_hours, y_step, where='post', color='k', lw=1.0)
    for s in (['R'] + list(custom_stages)):
        col = REM_COLOR if s == 'R' else stage_colors.get(s)
        for ei in np.where(hypno_labels == s)[0]:
            ax0.plot([ep_hours[ei], ep_hours[ei + 1]], [stage_y[s], stage_y[s]],
                     color=col, lw=1.0, solid_capstyle='butt')
    ax0.set_yticks(ytick_pos)
    ax0.set_yticklabels(ytick_labels, fontsize=8)
    ax0.set_ylim(min(ytick_pos) - 0.8, max(ytick_pos) + 0.8)
    ax0.set_xlim(0, ep_hours[-1] if n_ep else 1)
    ax0.xaxis.set_visible(False)
    for sp in ['top', 'right', 'bottom']:
        ax0.spines[sp].set_visible(False)

    # Spectrogram.
    im = ax1.pcolormesh(t, f, Sxx, norm=norm, cmap=spec_cmap, antialiased=True, shading='auto')
    ax1.set_xlim(0, t.max() if len(t) else 1)
    ax1.set_ylabel('Frequency [Hz]', fontsize=11)
    ax1.set_xlabel('Time [hrs]', fontsize=11)
    ax1.tick_params(labelsize=8)
    # Colorbar in absolute dB. The per-channel colour scale is relative (auto-scaled to each
    # channel), so a globally attenuated / low-gain channel looks structurally normal; showing the
    # absolute dB axis exposes it -- its whole scale sits tens of dB below a healthy channel's.
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label('Power (dB, µV²/Hz)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    return fig


def plot_hypnospectrogram(sig_uV, sf, hypno_labels, custom_stages, fmin=0.5, fmax=40, trimperc=2.5):
    """Multitaper spectrogram (YASA core) with a custom hypnogram strip that supports AASM +
    project-specific custom stages. The spectrogram is stage-agnostic (reuses yasa's
    spectrogram_lspopt); only the top hypnogram band is drawn here, so any custom label is
    plottable (YASA's own Hypnogram overlay hard-rejects non-AASM labels). Returns a Figure.
    Thin wrapper: compute then render (see compute_hypnospectrogram_sxx / render_hypnospectrogram)."""
    f, t, Sxx = compute_hypnospectrogram_sxx(sig_uV, sf, fmin=fmin, fmax=fmax)
    return render_hypnospectrogram(f, t, Sxx, hypno_labels, custom_stages, trimperc=trimperc)


def plot_median_hypnospectrogram(f, t, sxx_stack, hypno_labels, custom_stages, trimperc=2.5):
    """Across-channel median spectrogram for the Overview section of a dense montage.

    All channels are kept, including the flagged ones — the Overview plots show everything, and
    the median (not the mean) is what makes that safe: it tolerates up to half the channels being
    aberrant, so one dead or saturated electrode cannot dominate the single aggregate image.
    The median is taken on the dB spectrogram (log power), which is the scale-robust quantity;
    a linear-power mean would be driven by the noisiest channels. No per-channel normalisation:
    a dead channel has ~zero variance, so z-scoring would explode instead of being absorbed.

    The dead-epoch greying of render_hypnospectrogram then applies to the median, which gives
    majority-of-channels semantics for free: a column is grey only when most channels are dead.
    """
    Sxx_med = np.median(np.stack(sxx_stack, axis=0), axis=0).astype(np.float64)
    return render_hypnospectrogram(f, t, Sxx_med, hypno_labels, custom_stages, trimperc=trimperc)


# --- Per-participant "All electrodes" figures (Overview section) --------------------------------
# The dataset_overview.html "All electrodes" section pools every participant's channels together.
# The two figures below are the per-participant analogue, restricted to a single recording's own
# channels, so an outlier channel is spotted at a glance in that participant's report. The key
# metrics and their labels are shared with generate_dataset_overview through these module-level
# constants so the dataset-level and per-participant views stay in sync.
OVERVIEW_KEY_METRICS = [
    "std_uV", "flat_pct", "bounds_pct", "hist_extreme_pct",
    "p99_abs_uV", "p999_abs_uV",
]
OVERVIEW_METRIC_LABELS = {
    "std_uV": "Std dev (µV)",
    "flat_pct": "Flat signal (%)",
    "bounds_pct": "At EDF bounds (%)",
    "hist_extreme_pct": "Extreme histogram (%)",
    "p99_abs_uV": "p99 |amplitude| (µV)",
    "p999_abs_uV": "p99.9 |amplitude| (µV)",
    "mean_uV": "Mean (µV)",
    "kurtosis": "Kurtosis (Fischer)",
    "skewness": "Skewness",
}


def _overview_flag_legend(fig):
    """Flagged (red) / not-flagged (grey) dot legend, upper-right of the figure."""
    handles = [
        matplotlib.lines.Line2D([0], [0], marker="o", color="w",
                                markerfacecolor="#c0392b", markersize=8, label="Flagged"),
        matplotlib.lines.Line2D([0], [0], marker="o", color="w",
                                markerfacecolor="#555555", markersize=8, label="Not flagged"),
    ]
    fig.legend(handles=handles, loc="upper right", fontsize=9, framealpha=0.9)


def plot_participant_pooled_boxplots(df_sub, df_stage_sub=None, custom_stages=None):
    """Pooled distribution of each metric across THIS participant's channels (one box per metric).

    Each dot is one channel of the montage (red = flagged, grey = within thresholds). When the
    participant has a valid hypnogram (df_stage_sub non-empty), a per-stage inset in the upper-right
    of each subplot shows the same metric split by sleep stage (a box across channels per stage).
    Per-participant analogue of generate_dataset_overview's pooled boxplots.
    """
    custom_stages = custom_stages or []
    _sy, SS_COLORS, _ytp, _ytl = custom_stage_style(custom_stages)
    STAGE_ORDER = ["W", "N1", "N2", "N3", "R"] + list(custom_stages)
    file_id = str(df_sub["file_id"].iloc[0]) if "file_id" in df_sub.columns and len(df_sub) else ""

    present_stages = []
    if df_stage_sub is not None and len(df_stage_sub) > 0:
        present_stages = [s for s in STAGE_ORDER if s in df_stage_sub["stage"].values]

    ncols = min(3, len(OVERVIEW_KEY_METRICS))
    nrows = int(np.ceil(len(OVERVIEW_KEY_METRICS) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows), squeeze=False)
    axes = axes.flatten()
    for ax, m in zip(axes, OVERVIEW_KEY_METRICS):
        if m not in df_sub.columns:
            ax.set_visible(False)
            continue
        valid = df_sub[m].notna()
        vals = df_sub.loc[valid, m].values
        flags = (df_sub.loc[valid, "exclude"].fillna(False).astype(bool).values
                 if "exclude" in df_sub.columns else np.zeros(len(vals), dtype=bool))
        if len(vals) == 0:
            ax.set_visible(False)
            continue
        bp = ax.boxplot(
            [vals], patch_artist=True, showfliers=False,
            medianprops=dict(color="#c0392b", linewidth=2),
            whiskerprops=dict(linewidth=1.2), capprops=dict(linewidth=1.2),
        )
        bp["boxes"][0].set_facecolor("#d4e8f5")
        bp["boxes"][0].set_alpha(0.8)
        # One dot per channel, jittered around x=1, red when the channel is flagged.
        # Fixed seed: the same jitter on every run, so a report can be compared with an earlier one.
        jitter = np.random.default_rng(0).uniform(-0.07, 0.07, len(vals))
        for v, j, f in zip(vals, jitter, flags):
            ax.scatter(1 + j, v, color="#c0392b" if f else "#555555", s=28, alpha=0.75,
                       zorder=3, linewidths=0)
        ax.set_title(OVERVIEW_METRIC_LABELS.get(m, m), fontsize=11, fontweight="bold")
        ax.set_xticks([])
        for sp in ["right", "top"]:
            ax.spines[sp].set_visible(False)

        # Per-stage inset (upper-right), same layout as the dataset overview: widen the x-limits so
        # the main box (x=1) is pushed left and the inset gets the right ~half of the axes.
        if present_stages and m in df_stage_sub.columns:
            ax.set_xlim(0.65, 2.6)
            ax_in = ax.inset_axes([0.42, 0.20, 0.56, 0.74])
            groups_in = [df_stage_sub[df_stage_sub["stage"] == s][m].dropna().values
                         for s in present_stages]
            bp_in = ax_in.boxplot(
                groups_in, labels=present_stages, patch_artist=True, showfliers=False,
                medianprops=dict(color="#333", linewidth=1.5),
                whiskerprops=dict(linewidth=1.0), capprops=dict(linewidth=1.0),
            )
            for patch, stage in zip(bp_in["boxes"], present_stages):
                patch.set_facecolor(SS_COLORS[stage])
                patch.set_alpha(0.85)
            ax_in.tick_params(axis="x", labelsize=7)
            ax_in.tick_params(axis="y", labelsize=7)
            for sp in ["right", "top"]:
                ax_in.spines[sp].set_visible(False)
            ax_in.patch.set_facecolor("#fafafa")
            ax_in.patch.set_alpha(0.9)

    for ax in axes[len(OVERVIEW_KEY_METRICS):]:
        ax.set_visible(False)
    _overview_flag_legend(fig)
    fig.suptitle(f"{file_id}: All electrodes pooled (n={len(df_sub)} channels)",
                 fontsize=12, y=1.02)
    plt.tight_layout()
    return fig


def plot_participant_electrode_scatter(df_sub, custom_stages=None):
    """One subplot per metric; x-axis = electrode (montage order), one point per electrode.

    Point is red when the channel is flagged, grey otherwise, so an outlier channel is obvious at a
    glance for this participant. Per-participant analogue of generate_dataset_overview's grouped
    boxplots (a boxplot degenerates here since each electrode carries a single value).
    """
    file_id = str(df_sub["file_id"].iloc[0]) if "file_id" in df_sub.columns and len(df_sub) else ""
    channels = list(df_sub["channel"])
    x = np.arange(len(channels))
    ncols = min(3, len(OVERVIEW_KEY_METRICS))
    nrows = int(np.ceil(len(OVERVIEW_KEY_METRICS) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows), squeeze=False)
    axes = axes.flatten()
    for ax, m in zip(axes, OVERVIEW_KEY_METRICS):
        if m not in df_sub.columns:
            ax.set_visible(False)
            continue
        vals = df_sub[m].values
        flags = (df_sub["exclude"].fillna(False).astype(bool).values
                 if "exclude" in df_sub.columns else np.zeros(len(vals), dtype=bool))
        colors = ["#c0392b" if f else "#555555" for f in flags]
        ax.scatter(x, vals, color=colors, s=40, alpha=0.85, zorder=3, linewidths=0)
        ax.set_title(OVERVIEW_METRIC_LABELS.get(m, m), fontsize=11, fontweight="bold")
        ax.set_xticks(x)
        rot = 30 if len(channels) > 4 else 0
        ax.set_xticklabels(channels, rotation=rot, fontsize=9, ha="right" if rot else "center")
        ax.margins(x=0.05)
        for sp in ["right", "top"]:
            ax.spines[sp].set_visible(False)
    for ax in axes[len(OVERVIEW_KEY_METRICS):]:
        ax.set_visible(False)
    _overview_flag_legend(fig)
    fig.suptitle(f"{file_id}: metric by electrode", fontsize=12, y=1.02)
    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------------------------------------
# One recording's computation (lifted from run_analysis of 6_quality_overview_voila: the lines are tool 6's
# own, its progress-bar updates became the optional `status` / `on_channel` callbacks and its notebook
# prints the `warn` callback)
# ---------------------------------------------------------------------------------------------------------
def hypno_vector(expert_hypno, custom_stages):
    """Stage labels -> numeric stage vector (W 0, N1 1, N2 2, N3 3, R 4, custom stages 5, 6, ...).
    Returns (hypno_vec, hypno_warning, hypno_warning_plain); hypno_vec is None when more than 10 % of the
    epochs carry an unrecognised label (the per-stage metrics and the spectrogram are then skipped)."""
    hypno_vec = None
    hypno_warning = None
    hypno_warning_plain = None
    try:
        stage_map = {'W': 0, 'N1': 1, 'N2': 2, 'N3': 3, 'R': 4}
        for _i, _cs in enumerate(custom_stages):
            stage_map[_cs] = 5 + _i
        hypno_vec = np.full(len(expert_hypno), np.nan)
        for stage, val in stage_map.items():
            hypno_vec[expert_hypno == stage] = val
        # MT (movement time) is tolerated: YASA treats it as artifact (NaN), no warning needed.
        # Any other unrecognised label means the hypnogram was not converted to AASM convention.
        _tolerated = {'MT'}
        _unknown_mask = np.array([s not in stage_map and s not in _tolerated for s in expert_hypno])
        _unknown_labels = sorted(set(expert_hypno[_unknown_mask]))
        _n_unknown = int(_unknown_mask.sum())
        _pct_unknown = _n_unknown / len(expert_hypno) * 100
        if _pct_unknown > 10:
            hypno_vec = None
            hypno_warning = (
                f'Hypnogram labels not AASM: {_n_unknown}/{len(expert_hypno)} epochs '
                f'({_pct_unknown:.1f}%) have unrecognised labels '
                f'({", ".join(_unknown_labels)}). Spectrogram skipped. '
                f'Run <i>3_remap_hypno_voila</i> to convert labels to W/N1/N2/N3/R.'
            )
            hypno_warning_plain = (
                f'Hypnogram labels not AASM: {_n_unknown}/{len(expert_hypno)} epochs '
                f'({_pct_unknown:.1f}%) have unrecognised labels '
                f'({", ".join(_unknown_labels)}). Spectrogram skipped.'
            )
        elif _unknown_labels:
            hypno_warning = (
                f'{_n_unknown} epoch(s) ({_pct_unknown:.1f}%) have unrecognised label(s): '
                f'{", ".join(_unknown_labels)}. These epochs will be treated as artifact '
                f'in the spectrogram.'
            )
            hypno_warning_plain = (
                f'{_n_unknown} epoch(s) ({_pct_unknown:.1f}%) have unrecognised label(s): '
                f'{", ".join(_unknown_labels)}. Treated as artifact in the spectrogram.'
            )
    except Exception as e:
        hypno_warning = f'Error loading hypnogram: {e}'
        hypno_warning_plain = f'Error loading hypnogram: {e}'
    return hypno_vec, hypno_warning, hypno_warning_plain


def prepare_signal(raw, sf, file_id, do_resample, target_freq, do_highpass, hp_freq_val, warn=print, status=None):
    """Tool 6's two optional signal steps, in place on `raw`: resampling (never upsamples) THEN high-pass
    (analysis only). `status(text)` is called before each step, `warn(text)` when a step fails (the run goes
    on without it). Returns the sampling frequency after the resampling."""
    # --- Optional resampling (speeds up processing; off by default) ---
    # Resampling changes the signal, so shape metrics (flat_pct, histograms) and the
    # spectrogram are computed on the resampled data. Guarded to never upsample.
    # Done BEFORE the high-pass: mne.io.Raw.resample already anti-aliases (FFT method),
    # so resampling first introduces no aliasing, and any resampling edge transients
    # stay at the recording ends (wake). See the AFTER-resampling note below on order.
    if do_resample and target_freq is not None and target_freq < sf:
        if status:
            status(f'{file_id} · resampling to {target_freq} Hz…')
        try:
            raw.resample(target_freq, npad='auto', verbose=False)
            sf = raw.info['sfreq']
        except Exception as _rs_err:
            warn(f'  ⚠ {file_id}: resample to {target_freq} Hz failed: {_rs_err}')

    # --- Optional high-pass (analysis only; applied AFTER resampling) ---
    # Removes DC offset / slow drift. Off by default; useful to harmonise a heterogeneous
    # dataset (or for DC-coupled data, e.g. Curry). Order vs resampling does not change the
    # result: both operations are linear, so swapping them leaves the signal interior
    # identical to <0.2 uV (verified empirically) — the resampler's anti-alias step rings
    # the same on an abrupt amplitude change either way, since a low-corner high-pass does
    # not touch the sharp edge that causes the ring. High-passing the already-downsampled
    # signal just runs on far fewer samples, so it is markedly faster.
    if do_highpass and hp_freq_val is not None:
        if status:
            status(f'{file_id} · high-pass {hp_freq_val:g} Hz…')
        try:
            raw.filter(l_freq=hp_freq_val, h_freq=None, verbose=False)
        except Exception as _hp_err:
            warn(f'  ⚠ {file_id}: high-pass ({hp_freq_val:g} Hz) failed: {_hp_err}')
    return sf


def check_hypno_length(raw, sf, expert_hypno, hypno_vec, hypno_warning, hypno_warning_plain):
    """Tool 6's hypnogram length rule: a score past the last complete 30 s epoch is dropped, any other
    mismatch disables the per-stage metrics and the spectrogram (hypno_vec None).
    Returns (expert_hypno, hypno_vec, hypno_warning, hypno_warning_plain)."""
    if hypno_vec is not None:
        expected = int(np.floor(raw.n_times / sf / 30))
        n_hyp = len(hypno_vec)
        overhang = n_hyp - expected
        if 0 < overhang <= 1:
            # Recording ends mid-epoch: the last epoch is partial (< 30 s) but was still
            # scored. Drop the overhanging score (no full 30 s of signal to analyse).
            last_sec = raw.n_times / sf - expected * 30
            expert_hypno = expert_hypno[:expected]
            hypno_vec = hypno_vec[:expected]
            hypno_warning = (
                f'Recording ends mid-epoch: the last epoch is only {last_sec:.1f}s '
                f'(< 30s) so its score is skipped. Hypnogram trimmed from {n_hyp} to '
                f'{expected} epochs to match the complete 30s signal epochs.'
            )
            hypno_warning_plain = (
                f'Last epoch {last_sec:.1f}s (< 30s) skipped; hypnogram trimmed '
                f'{n_hyp} -> {expected} to match complete 30s signal epochs.'
            )
        elif overhang != 0:
            hypno_warning = (
                f'Length mismatch: EEG has {expected} 30s-epochs, '
                f'hypnogram has {n_hyp}. Spectrogram skipped.'
            )
            hypno_warning_plain = (
                f'Length mismatch: EEG has {expected} 30s-epochs, '
                f'hypnogram has {n_hyp}. Spectrogram skipped.'
            )
            hypno_vec = None
    return expert_hypno, hypno_vec, hypno_warning, hypno_warning_plain


def compute_channel_metrics(raw, sf, file_id, thresholds, hypno_vec, custom_stages, phys_bounds_by_name,
                            final_to_orig, on_channel=None):
    """Tool 6's metrics and flags of every channel on the whole night, plus the per-stage metrics when a
    valid hypnogram vector is given (the 30 s epochs of each stage, concatenated). `phys_bounds_by_name`
    {original channel name: (min, max) µV} and `final_to_orig` {harmonised name: original name} give the EDF
    physical range of each channel (bounds_pct). `on_channel(ci, ch)` is called before each channel.
    Returns (ch_metrics, ch_flags, file_rows, file_stage_rows): the rows of {file_id}_quality_metrics.tsv
    and of {file_id}_quality_by_stage.tsv."""
    file_rows = []        # this file's per-channel metric rows → {file_id}_quality_metrics.tsv
    file_stage_rows = []  # this file's per-stage rows → {file_id}_quality_by_stage.tsv (if hypnogram)

    # --- Compute metrics for all channels (unfiltered signal) ---
    # Memory: a full-night, high-rate recording with many channels (e.g. 32 EEG @ 1024 Hz)
    # is ~10 GiB as one float64 array. We never materialise the whole thing — each channel
    # is read on demand with raw.get_data(picks=[ci]) (a cheap slice of the preloaded data),
    # and the shared amplitude limits are derived from the per-channel p99.9 metric below.

    ch_metrics = {}
    ch_flags = {}
    epoch_len = int(30 * sf)
    for ci, ch in enumerate(raw.ch_names):
        if on_channel:
            on_channel(ci, ch)
        sig_uV = raw.get_data(picks=[ci])[0] * 1e6
        orig_base = final_to_orig.get(ch, ch)
        phys_min, phys_max = phys_bounds_by_name.get(orig_base, (np.nan, np.nan))
        m = compute_signal_metrics(sig_uV, phys_min, phys_max)
        ch_metrics[ch] = m
        ch_flags[ch] = flag_channel(m, thresholds)
        file_rows.append({
            'file_id': file_id,
            'channel': ch,
            'mean_uV': m['mean_uV'],
            'std_uV': m['std_uV'],
            'kurtosis': m['kurtosis'],
            'skewness': m['skewness'],
            'p99_abs_uV': m['p99_abs_uV'],
            'p999_abs_uV': m['p999_abs_uV'],
            'flat_pct': m['flat_pct'],
            'bounds_pct': m['bounds_pct'],
            'hist_extreme_pct': m['hist_extreme_pct'],
            'n_peaks': m['n_peaks'],
            'suspect_reason': '; '.join(ch_flags[ch]),
            'exclude': bool(ch_flags[ch]),
        })

        # Per-stage key metrics: split the signal into 30s epochs by hypnogram stage.
        # Skipped if no valid hypnogram (not found, length mismatch, or bad labels).
        if hypno_vec is not None:
            _stage_pairs = [(0, 'W'), (1, 'N1'), (2, 'N2'), (3, 'N3'), (4, 'R')]
            _stage_pairs += [(5 + _i, _cs) for _i, _cs in enumerate(custom_stages)]
            for stage_int, stage_str in _stage_pairs:
                stage_indices = np.where(hypno_vec == stage_int)[0]
                if len(stage_indices) == 0:
                    continue
                stage_samples = []
                for ep_idx in stage_indices:
                    start_samp = ep_idx * epoch_len
                    end_samp = start_samp + epoch_len
                    if end_samp <= len(sig_uV):
                        stage_samples.append(sig_uV[start_samp:end_samp])
                if not stage_samples:
                    continue
                stage_sig = np.concatenate(stage_samples)
                sm = compute_signal_metrics(stage_sig, phys_min, phys_max)
                file_stage_rows.append({
                    'file_id': file_id,
                    'channel': ch,
                    'stage': stage_str,
                    'mean_uV': sm['mean_uV'],
                    'std_uV': sm['std_uV'],
                    'flat_pct': sm['flat_pct'],
                    'bounds_pct': sm['bounds_pct'],
                    'hist_extreme_pct': sm['hist_extreme_pct'],
                    'p99_abs_uV': sm['p99_abs_uV'],
                    'p999_abs_uV': sm['p999_abs_uV'],
                })
    return ch_metrics, ch_flags, file_rows, file_stage_rows


def display_limits(ch_metrics, ch_flags):
    """Shared amplitude limits of the time-series and histogram figures.
    Returns (y_lim_ts, x_lim_hist, hist_y_max)."""
    # Shared time-series / histogram amplitude limits from the per-channel p99.9 metric
    # (avoids materialising the whole recording just to take a global percentile).
    _p999 = [m['p999_abs_uV'] for m in ch_metrics.values()]
    y_lim_ts = float(max(_p999)) if _p999 else 1.0
    x_lim_hist = y_lim_ts

    # Shared Y-axis for histograms: computed from non-suspect channels only.
    # A flat or dead channel concentrates all samples in 1-2 bins, giving a
    # y_max that would crush the distributions of healthy channels.
    # Flagged channels are plotted with their own auto-scale (y_max=None) to
    # avoid their distribution being clipped by the healthy-channel scale.
    healthy_chs = [ch for ch, f in ch_flags.items() if not f]
    if healthy_chs:
        hist_y_max = max(float(np.max(ch_metrics[ch]['histo'])) for ch in healthy_chs)
    elif ch_metrics:
        hist_y_max = max(float(np.max(m['histo'])) for m in ch_metrics.values())
    else:
        hist_y_max = None
    return y_lim_ts, x_lim_hist, hist_y_max


def compute_channel_spectrograms(raw, sf, file_id, hypno_vec, warn=print, on_channel=None):
    """Each channel's dB hypnospectrogram (display band-pass 0.1-40 Hz, 30 s multitaper windows, float32),
    computed once for the per-channel figures and the across-channel median. Nothing without a valid
    hypnogram vector. `on_channel(ci, ch, done)` is called before (done False) and after (done True) each
    channel. Returns (ch_sxx, sxx_f, sxx_t)."""
    # --- Spectrograms: computed once here, used twice ---
    # Each channel's dB spectrogram is needed both for its own report figure and for the
    # across-channel median in the Overview. The median needs every channel at once, so the
    # stack has to be held in memory anyway (~300 MB float32 on a 44 ch / 12 h montage);
    # given that, recomputing it in the report loop would be pure waste (~8 min per file).
    # float32 halves the stack and is far below the display's dB resolution.
    # Frequency resolution is set by the 30 s window (1/30 Hz), not by sfreq, so resampling
    # shrinks the compute time but not this stack.
    ch_sxx = {}
    sxx_f, sxx_t = None, None
    if hypno_vec is not None:
        for ci, ch in enumerate(raw.ch_names):
            if on_channel:
                on_channel(ci, ch, False)
            try:
                _sig = raw.get_data(picks=[ci])[0] * 1e6
                _sig_filt = mne.filter.filter_data(
                    _sig[np.newaxis, :], sfreq=sf, l_freq=0.1, h_freq=40,
                    method='fir', phase='zero-double', fir_window='hamming', verbose=False
                )[0]
                del _sig
                sxx_f, sxx_t, _sxx = compute_hypnospectrogram_sxx(
                    _sig_filt, sf, fmin=0.5, fmax=40)
                del _sig_filt
                ch_sxx[ch] = _sxx.astype(np.float32)
                del _sxx
            except Exception as _sx_err:
                warn(f'  ⚠ {file_id}/{ch}: spectrogram computation failed: {_sx_err}')
            if on_channel:
                on_channel(ci, ch, True)
    return ch_sxx, sxx_f, sxx_t
