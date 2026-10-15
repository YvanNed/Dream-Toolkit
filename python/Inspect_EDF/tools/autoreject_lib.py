"""Automatic epoch rejection of tool 8bis (8bis_reject_automatically_voila), shared with tool 0.

Moved verbatim from tool 8bis (its shared-function cell, and a few lines of its Scan and Run cells): tool 8bis
and tool 0's step 7 call the same code, so they take the same decision by construction. Edit here, never in a
notebook. The Curry twin of tool 8bis (a verbatim copy of the notebook) imports this file as well.

The rule is channel-first then epoch: drop a channel flagged in more than channel_pct of the reference-stage
epochs, then reject every epoch flagged in more than epoch_pct of the remaining good channels.
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

import qc_rejected_epochs_lib as L          # METHOD_ORDER, custom-stage styling (same folder)


# ---------------------------------------------------------------------------------------------------------
# Decision rule, tables and figures (tool 8bis's shared-function cell)
# ---------------------------------------------------------------------------------------------------------
AASM_STAGES = ['W', 'N1', 'N2', 'N3', 'R']

# ---- Editable defaults (channel-first, then epoch) ----
DEFAULT_CHANNEL_REJECT_PCT = 20.0   # drop a channel flagged in > this % of its reference-stage epochs

DEFAULT_EPOCH_REJECT_PCT   = 20.0   # drop an epoch flagged in > this % of the remaining good channels

def read_event_counts(path):
    """Return {canonical_type: n_events} from a tool-7 {file_id}_event_counts.tsv (raw per-type event
    count for the whole file), or {} if the sidecar is absent/unreadable. Feeds the per-event
    participant-exclusion threshold."""
    try:
        return event_counts_from_table(pd.read_csv(path, sep='\t'))
    except Exception:
        return {}


def event_counts_from_table(d):
    """{canonical_type: n_events} from a tool-7 event-counts table (read_event_counts reads it from its file)."""
    return {str(r['event_type']): int(r['n_events']) for _, r in d.iterrows()}

def build_pair_matrix(tsv_path, methods=None, event_master=True, event_types=None, event_flags_path=None):
    """Recompose a boolean (epoch x channel) rejection matrix from a tool-7
    {file_id}_epoch_channel_rejection.tsv, keeping ONLY the selected flagging methods (their per-method
    flag_<method> columns OR'd together), plus a per-epoch stage Series and a small meta dict.

    methods          : methods to include (subset of L.METHOD_ORDER); None = every method present.
    event_master     : include the 'event' method at all (only when 'event' is also in `methods`).
    event_types      : canonical event types allowed to reject (sub-selection); needs event_flags_path.
    event_flags_path : the tool-7 {file_id}_event_epoch_flags.tsv (one evt_<type> column per epoch).
                       Absent -> fall back to the TSV's OR'd flag_event (all event types, no sub-select).

    With methods=None + event_master=True the union of every flag_<method> equals the stored reject_any,
    so a default (all-methods) run stays byte-identical to the previous tool."""
    df = pd.read_csv(tsv_path, sep='\t')
    event_flags = None
    if event_master and (methods is None or 'event' in methods):
        if event_flags_path is not None and Path(event_flags_path).exists():
            event_flags = pd.read_csv(event_flags_path, sep='\t')
    return pair_matrix_from_tables(df, methods, event_master, event_types, event_flags)


def pair_matrix_from_tables(df, methods=None, event_master=True, event_types=None, event_flags=None):
    """build_pair_matrix on tables already in memory (tool 0 passes its step 6's): `df` = the
    {file_id}_epoch_channel_rejection.tsv table, `event_flags` = the {file_id}_event_epoch_flags.tsv table,
    or None when there is none (the OR'd flag_event is then used: all event types, no sub-selection)."""
    channels = list(pd.unique(df['channel']))               # epoch-major rows -> tool-7 channel order
    epoch_index = sorted(df['epoch_idx'].unique())
    stage_by_epoch = (df.drop_duplicates('epoch_idx').set_index('epoch_idx')['stage']
                      .astype(str).reindex(epoch_index))
    if methods is None:
        methods = [m for m in L.METHOD_ORDER if ('flag_' + m) in df.columns]
    piv = pd.DataFrame(False, index=epoch_index, columns=channels)

    def _or_flag(col):
        p = (df.pivot(index='epoch_idx', columns='channel', values=col)
             .reindex(index=epoch_index, columns=channels).fillna(0).astype(bool))
        return piv | p

    used_methods = []
    for m in methods:
        if m == 'event':
            continue                                        # event handled separately (epoch-level)
        if ('flag_' + m) in df.columns:
            piv = _or_flag('flag_' + m)
            used_methods.append(m)

    warn = None
    used_event_types = []
    if event_master and 'event' in methods:
        if event_flags is not None:
            ef = event_flags.set_index('epoch_idx')
            evt_cols = [f'evt_{t}' for t in (event_types or []) if f'evt_{t}' in ef.columns]
            used_event_types = [c[len('evt_'):] for c in evt_cols]
            ev_mask = (ef[evt_cols].astype(bool).any(axis=1) if evt_cols
                       else pd.Series(False, index=ef.index))
            # broadcast the epoch-level event flag across all channels (matches flag_event's storage)
            ev_col = ev_mask.reindex(epoch_index).fillna(False).values[:, np.newaxis]
            ev_df = pd.DataFrame(np.repeat(ev_col, len(channels), axis=1),
                                 index=epoch_index, columns=channels)
            piv = piv | ev_df
            used_methods.append('event')
        elif 'flag_event' in df.columns:
            piv = _or_flag('flag_event')                    # old data: no per-type sidecar
            used_methods.append('event')
            warn = 'no _event_epoch_flags.tsv (used flag_event = all event types)'

    meta = {'used_methods': used_methods, 'used_event_types': used_event_types, 'warn': warn}
    return piv, stage_by_epoch, meta

def auto_reject_decision(piv, stage_by_epoch, ref_stages, ch_pct, ep_pct):
    """Channel-first then epoch decision (ch_pct / ep_pct are fractions in [0, 1]):
    (1) drop a channel flagged in > ch_pct of the REFERENCE-stage epochs (judged there so that wake
        movement / EMG artefacts do not drop a channel that is clean in sleep);
    (2) reject ANY epoch, whatever its stage, flagged in > ep_pct of the remaining GOOD channels: the
        decision covers the whole night, the stages to analyse are chosen later (tool 9, statistics)."""
    channels = list(piv.columns)
    is_ref = stage_by_epoch.isin(ref_stages).values
    ref_epochs = list(piv.index[is_ref])
    n_ref = len(ref_epochs)
    badness = piv.loc[ref_epochs].mean(axis=0) if n_ref else pd.Series(0.0, index=channels)
    rejected_channels = [c for c in channels if float(badness[c]) > ch_pct]
    good = [c for c in channels if c not in rejected_channels]
    if good:
        frac = piv[good].mean(axis=1)
    else:
        frac = pd.Series(1.0, index=piv.index)                   # no good channel -> reject every epoch
    reject_epoch = frac > ep_pct                                  # every epoch of the file
    kept_epochs = [int(e) for e, r in reject_epoch.items() if not bool(r)]
    return {
        'channels': channels, 'badness': badness, 'rejected_channels': rejected_channels,
        'good_channels': good, 'ref_epochs': [int(e) for e in ref_epochs], 'n_ref': n_ref,
        'n_rejected_ref': int(reject_epoch.loc[ref_epochs].sum()) if n_ref else 0,
        'reject_epoch': reject_epoch, 'kept_epochs': kept_epochs, 'n_epochs': len(piv.index),
        'n_rejected_epochs': int(reject_epoch.sum()), 'stage_by_epoch': stage_by_epoch,
    }

def stages_in_order(dec, soi):
    """The reference stages first, then every other stage present (AASM order, then custom sorted)."""
    present = set(dec['stage_by_epoch'].dropna().astype(str))
    others = ([s for s in AASM_STAGES if s in present and s not in soi]
              + sorted(s for s in present if s not in AASM_STAGES and s not in soi))
    return list(soi) + others

def decision_summary_row(fid, dec, soi, ch_pct, ep_pct, methods_used='', event_types_used='',
                         excluded=False, exclude_reason=''):
    """One-row DataFrame = the durable per-participant record AND the global-summary row source.
    methods_used / event_types_used trace which flagging methods (and event types) drove the decision;
    excluded / exclude_reason mark a participant dropped by a per-event count threshold."""
    badness = dec['badness']
    badness_str = ';'.join(f'{c}:{float(badness[c]):.2f}' for c in dec['channels'])
    n_all, n_rej = dec['n_epochs'], dec['n_rejected_epochs']
    n_ref, n_rej_ref = dec['n_ref'], dec['n_rejected_ref']
    row = {
        'file_id': fid,
        'n_channels': len(dec['channels']),
        'n_channels_rejected': len(dec['rejected_channels']),
        'rejected_channels': ';'.join(dec['rejected_channels']),
        'good_channels': ';'.join(dec['good_channels']),
        'channel_badness_pct': badness_str,
        'channel_reject_pct': round(ch_pct * 100, 3),
        'epoch_reject_pct': round(ep_pct * 100, 3),
        'stages_of_interest': '+'.join(soi),
        'methods_used': methods_used,
        'event_types_used': event_types_used,
        'excluded': bool(excluded),
        'exclude_reason': exclude_reason,
        'n_epochs': n_all,                                # every epoch: the decision covers the night
        'n_epochs_rejected': n_rej,
        'n_epochs_kept': len(dec['kept_epochs']),
        'pct_epochs_rejected': round(100 * n_rej / n_all, 2) if n_all else 0.0,
        # the same counts over the reference stages (stages_of_interest), the ones that judged the channels
        'n_epochs_ref': n_ref,
        'n_epochs_rejected_ref': n_rej_ref,
        'pct_epochs_rejected_ref': round(100 * n_rej_ref / n_ref, 2) if n_ref else 0.0,
    }
    # per-stage counts, every stage present
    sbe, rej = dec['stage_by_epoch'].astype(str), dec['reject_epoch']
    for st in stages_in_order(dec, soi):
        eps = [e for e in rej.index if sbe.get(e) == st]
        row[f'n_{st}'] = len(eps)
        row[f'n_rejected_{st}'] = int(sum(bool(rej[e]) for e in eps))
    return pd.DataFrame([row])

def read_tool7_bad_channels(fif_folder, fid):
    """Channels tool 7 saved in the .fif but marked bad (deselected by the user, e.g. after tool 6):
    read from its sidecar {id}_preprocessing_params.json (key bad_channels; absent = none)."""
    try:
        with open(Path(fif_folder) / f'{fid}_preprocessing_params.json', encoding='utf-8') as f:
            return list(json.load(f).get('bad_channels') or [])
    except Exception:
        return []

def build_epoch_decision(fid, dec, piv):
    """One row per epoch of the file: in_scope (always True: the decision covers every stage), rejected
    and why: how many of the remaining good channels flagged it."""
    good = dec['good_channels']
    rej = dec['reject_epoch']
    n_flag = piv[good].sum(axis=1) if good else pd.Series(0, index=piv.index)
    rows = []
    for e in piv.index:
        e = int(e)
        rejected = bool(rej.get(e, False))
        if not rejected:
            reason = ''
        elif not good:
            reason = 'no good channel left'
        else:
            reason = f'{int(n_flag[e])}/{len(good)} good channels flagged'
        rows.append({'file_id': fid, 'epoch_idx': e, 'stage': dec['stage_by_epoch'].get(e),
                     'in_scope': True, 'rejected': rejected, 'reject_reason': reason})
    return pd.DataFrame(rows, columns=['file_id', 'epoch_idx', 'stage', 'in_scope', 'rejected', 'reject_reason'])

def build_channel_decision(fid, dec, ch_pct, tool7_bads):
    """One row per channel: the channels tool 7 flagged (dropped when badness > ch_pct) plus the ones
    tool 7 kept as bad (always dropped). badness_pct = % of reference-stage epochs flagged (NaN for tool-7 bads,
    which were never flagged)."""
    rows = []
    for c in dec['channels']:
        b = 100.0 * float(dec['badness'][c])
        dropped = c in dec['rejected_channels']
        rows.append({'file_id': fid, 'channel': c, 'dropped': dropped,
                     'drop_reason': f'badness {b:.1f} % > {100 * ch_pct:g} %' if dropped else '',
                     'badness_pct': round(b, 3)})
    for c in tool7_bads:
        if c not in dec['channels']:
            rows.append({'file_id': fid, 'channel': c, 'dropped': True, 'drop_reason': 'tool7_deselected',
                         'badness_pct': np.nan})
    return pd.DataFrame(rows, columns=['file_id', 'channel', 'dropped', 'drop_reason', 'badness_pct'])

def wrap_scroll(html_table):
    """Wrap a wide HTML table so it scrolls horizontally inside its container instead of overflowing
    the page to the right (applied to the notebook display and the global HTML report)."""
    return ('<div style="overflow-x:auto; max-width:100%">'
            + html_table.replace('<table', '<table style="border-collapse:collapse"', 1)
            + '</div>')

def per_stage_table_html(dec, soi):
    """HTML: rejected / kept epoch counts per stage, every stage present (reference stages marked)."""
    sbe, rej = dec['stage_by_epoch'].astype(str), dec['reject_epoch']
    rows = []
    for st in stages_in_order(dec, soi):
        st_epochs = [e for e in rej.index if sbe.get(e) == st]
        n = len(st_epochs); r = sum(int(rej[e]) for e in st_epochs)
        pct = f'{100 * r / n:.1f}%' if n else '-'
        label = f'{st} (reference)' if st in soi else st
        rows.append(f'<tr><td>{label}</td><td>{n}</td><td>{r}</td><td>{pct}</td></tr>')
    return ('<table border="1" cellpadding="4" style="border-collapse:collapse">'
            '<tr><th>Stage</th><th>Epochs</th><th>Rejected</th><th>% rejected</th></tr>'
            + ''.join(rows) + '</table>')

def decision_overview_html(fid, dec, soi, ch_pct, ep_pct, methods_used='', event_types_used='',
                           exclude_reason=''):
    """HTML summary block for the per-participant report. Traces the flagging methods (and event types)
    that drove the decision, and shows a red banner when the participant was excluded by an event
    count threshold."""
    badness = dec['badness']
    ch_lines = ''.join(
        f'<li>{c}: {100 * float(badness[c]):.1f}% flagged'
        + (' <b style="color:#c0392b">-> dropped</b>' if c in dec['rejected_channels'] else '')
        + '</li>' for c in dec['channels'])
    pct = 100 * dec['n_rejected_epochs'] / dec['n_epochs'] if dec['n_epochs'] else 0.0
    pct_ref = 100 * dec['n_rejected_ref'] / dec['n_ref'] if dec['n_ref'] else 0.0
    banner = (f'<p style="padding:6px 10px;background:#fdecea;border:1px solid #c0392b;color:#c0392b;">'
              f'<b>EXCLUDED</b>: {exclude_reason} (exclusion registry, 8bis_auto).</p>') if exclude_reason else ''
    methods_line = (f'<p>Methods used: <b>{methods_used or "(none)"}</b>'
                    + (f', event types: <b>{event_types_used}</b>' if event_types_used else '')
                    + '</p>')
    return (
        f'<h3>{fid} - automatic rejection</h3>'
        + banner + methods_line +
        f'<p>Reference stages (judge the channels): <b>{"+".join(soi)}</b>, channel threshold &gt; {100 * ch_pct:.0f}% '
        f'&mdash; epoch threshold &gt; {100 * ep_pct:.0f}% of good channels.</p>'
        f'<p>Channels ({len(dec["good_channels"])}/{len(dec["channels"])} kept):</p><ul>{ch_lines}</ul>'
        f'<p>Epochs (every stage): <b>{dec["n_epochs"]}</b>, rejected: <b>{dec["n_rejected_epochs"]}</b> '
        f'&mdash; kept: <b>{len(dec["kept_epochs"])}</b> ({pct:.1f}% rejected). In the reference stages: '
        f'{dec["n_rejected_ref"]} / {dec["n_ref"]} rejected ({pct_ref:.1f}%)</p>'
        + per_stage_table_html(dec, soi))

def plot_decision_heatmap(piv, dec, fid, custom_stages):
    """Channels x epochs flagged-pair heatmap with a hypnogram strip and a rejected-epoch strip.
    Rejected-channel rows get a red bold label; the reject strip covers every epoch."""
    channels = dec['channels']
    epochs = list(piv.index)
    ep_pos = {e: i for i, e in enumerate(epochs)}
    n_ep, n_ch = len(epochs), len(channels)
    M = piv[channels].values.T.astype(float)                     # (n_ch, n_ep) 0/1

    stage_y, _sc, ytick_pos, ytick_labels = L.custom_stage_style(list(custom_stages))
    sbe = dec['stage_by_epoch']
    hyp = np.array([stage_y.get(str(sbe[e]), np.nan) for e in epochs], dtype=float)

    rej = dec['reject_epoch']
    strip = np.zeros(n_ep)
    for e in epochs:
        strip[ep_pos[e]] = 1.0 if bool(rej[e]) else 0.0          # 1 = rejected (red), 0 = kept (green)

    # The channels x epochs matrix height scales with the channel count (like tool 7's heatmap); the
    # hypnogram and epoch-decision strips keep a FIXED height so they stay readable on dense montages.
    fig_w   = max(8, min(n_ep / 6, 20))
    h_hyp   = 1.0                                                 # inches, fixed
    h_strip = 0.4                                                 # inches, fixed
    h_mat   = min(max(0.22 * n_ch, 0.8), 12.0)                   # inches, ~0.22 in / channel row
    fig, (ax_h, ax_s, ax_m) = plt.subplots(
        3, 1, figsize=(fig_w, h_hyp + h_strip + h_mat + 0.7), sharex=True, layout='constrained',
        gridspec_kw={'height_ratios': [h_hyp, h_strip, h_mat]})

    ax_h.step(range(n_ep), hyp, where='mid', color='#333333', lw=0.8)
    for i, e in enumerate(epochs):          # REM in red at the line's width (the toolkit's discreet convention)
        if str(sbe[e]) == 'R':
            ax_h.plot([i - 0.5, i + 0.5], [stage_y['R']] * 2, color='#c0392b', lw=0.8, solid_capstyle='butt')
    ax_h.set_yticks(ytick_pos); ax_h.set_yticklabels(ytick_labels, fontsize=7)
    ax_h.set_ylabel('stage', fontsize=8)
    ax_h.set_title(f'{fid} - flagged (epoch x channel); red = dropped channel / rejected epoch', fontsize=10)

    ax_s.imshow(strip[np.newaxis, :], aspect='auto', vmin=0, vmax=2,
                cmap=ListedColormap(['#2e7d32', '#c0392b', '#d9d9d9']))
    ax_s.set_yticks([]); ax_s.set_ylabel('epoch', fontsize=8, rotation=0, ha='right', va='center')

    # 0 = not flagged, 1 = flagged (good channel, blue), 2 = flagged in a DROPPED channel (grey).
    # Only the flagged cells of a dropped channel turn grey, so which epochs were flagged stays visible.
    disp = M.copy()
    for i, c in enumerate(channels):
        if c in dec['rejected_channels']:
            disp[i, M[i, :] > 0] = 2.0
    ax_m.imshow(disp, aspect='auto', vmin=0, vmax=2, interpolation='nearest',
                cmap=ListedColormap(['#f7f7f7', '#34495e', '#d9d9d9']))
    ax_m.set_yticks(range(n_ch)); ax_m.set_yticklabels(channels, fontsize=8)
    for tick, c in zip(ax_m.get_yticklabels(), channels):
        if c in dec['rejected_channels']:
            tick.set_color('#c0392b'); tick.set_fontweight('bold')
    ax_m.set_xlabel('epoch index', fontsize=8)
    return fig


# ---------------------------------------------------------------------------------------------------------
# Lines of 8bis's Scan and Run cells, as functions (tool 0's step 7 calls them too)
# ---------------------------------------------------------------------------------------------------------
def ordered_stages(stages):
    """The stages present, AASM order first, then the custom ones sorted (8bis's stage checkboxes)."""
    return [s for s in AASM_STAGES if s in stages] + sorted(s for s in stages if s not in AASM_STAGES)


def default_reference_stages(ordered):
    """8bis's default reference stages: N2 / N3 / R when present, else every stage."""
    prefer = [s for s in ordered if s in ('N2', 'N3', 'R')]
    return prefer if prefer else list(ordered)


def event_exclusion_reason(counts, event_thresholds):
    """Per-event participant exclusion: the reason when the file has more events of a ticked type than its
    threshold (raw event count, `counts` as read_event_counts returns it), else ''."""
    excl_reason = ''
    for t, thr in event_thresholds.items():
        n = counts.get(t)
        if n is not None and n > thr:
            excl_reason = f'{n} {t} events > {thr}'
            break
    return excl_reason


def registry_reason(excl_reason, dec):
    """Why 8bis excludes this participant in the registry (8bis_auto), '' when it does not: its event-count
    exclusion, else every channel dropped, else every reference-stage epoch rejected."""
    auto_reason = excl_reason
    if not auto_reason and not dec['good_channels']:
        auto_reason = 'all channels dropped'
    elif not auto_reason and dec['n_ref'] and dec['n_rejected_ref'] == dec['n_ref']:
        auto_reason = 'all reference-stage epochs rejected'
    return auto_reason
