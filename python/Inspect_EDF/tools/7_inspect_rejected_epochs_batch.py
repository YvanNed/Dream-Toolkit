#!/usr/bin/env python
"""
Tool 7 (batch twin) — Section-2 QC report of rejected epochs, per participant AND database-level.

Runs the same per-stage analysis as 7_inspect_rejected_epochs_voila.ipynb over a whole derivatives
tree (no interaction), then aggregates a database-level view. Reuses qc_rejected_epochs_lib so the
plots stay identical to the notebook. The interactive per-epoch navigator / override / clean-epo save
are notebook-only; this batch produces reports only and never modifies tool-6 outputs.

Usage
-----
    python 7_inspect_rejected_epochs_batch.py <derivatives_root> [--no-1f] [--limit N] [--out DIR]

    <derivatives_root>  folder holding the tool-6 *_all-epo.fif files (searched recursively)
    --no-1f             skip specparam 1/f fitting (much faster; MAE/R² columns left blank)
    --limit N           process only the first N participants (quick test)
    --out DIR           output folder (default: <derivatives_root>/qc2b_reports)

Outputs (in the output folder)
------------------------------
    {file_id}_qc2b_report.html            per-participant Section-2 report
    qc2b_database_report.html             database-level aggregate report
    qc2b_database_rejection_summary.tsv   one row per file_id x stage (counts + %)
    qc2b_failed.tsv                       participants that could not be processed (if any)
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qc_rejected_epochs_lib as L


def _default_thresholds():
    """A thresholds dict (tool-6 defaults) shaped like load_params()'s output, for the DB reference lines."""
    d = L.DEFAULT_THRESHOLDS
    return {
        'amplitude_ptp_uV': dict(d['amplitude_ptp_uV']),
        'flat_ptp_uV': d['flat_ptp_uV'],
        'gradient_uV_per_sample': d['gradient_uV_per_sample'],
        '1f_mae_max': d['1f_mae_max'],
        '1f_r2_min': d['1f_r2_min'],
    }


def _plot_per_participant_bar(summ):
    """Overall % rejected per participant (sorted). Returns a Figure."""
    agg = summ.groupby('file_id').agg(nr=('n_rejected', 'sum'), nt=('n_total', 'sum'))
    agg['pct'] = 100.0 * agg['nr'] / agg['nt']
    agg = agg.sort_values('pct', ascending=False)
    fig, ax = plt.subplots(figsize=(max(6, 0.35 * len(agg) + 2), 4.2))
    ax.bar(range(len(agg)), agg['pct'].values, color='#807dba')
    ax.set_xticks(range(len(agg)))
    ax.set_xticklabels(agg.index, rotation=90, fontsize=7)
    ax.set_ylabel('% epochs rejected')
    ax.set_title('Overall rejection rate per participant', fontsize=10)
    ax.grid(True, axis='y', alpha=0.2)
    fig.tight_layout()
    return fig


def _plot_per_stage_box(summ, stage_order):
    """Distribution of per-participant % rejected, per stage. Returns a Figure."""
    fig, ax = plt.subplots(figsize=(7, 4.2))
    data, labels = [], []
    for st in stage_order:
        vals = summ.loc[summ['stage'] == st, 'pct_rejected'].values
        if vals.size:
            data.append(vals)
            labels.append(st)
    if data:
        ax.boxplot(data, showfliers=False)
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels)
        for i, vals in enumerate(data, 1):
            ax.scatter(np.full(vals.size, i) + np.random.uniform(-0.12, 0.12, vals.size),
                       vals, s=8, color='#c0392b', alpha=0.5)
    ax.set_ylabel('% epochs rejected (per participant)')
    ax.set_title('Rejection rate by stage across participants', fontsize=10)
    ax.grid(True, axis='y', alpha=0.2)
    fig.tight_layout()
    return fig


def _db_rejection_table_html(summ, stage_order):
    """Mean ± sd of per-participant % rejected, per stage (+ overall). HTML string."""
    header = ['Stage', 'N participants', 'Mean % rejected', 'SD', 'Median']
    rows = []
    for st in stage_order:
        v = summ.loc[summ['stage'] == st, 'pct_rejected'].values
        if v.size:
            rows.append([st, str(v.size), f'{v.mean():.1f}%', f'{v.std():.1f}', f'{np.median(v):.1f}%'])
    agg = summ.groupby('file_id').agg(nr=('n_rejected', 'sum'), nt=('n_total', 'sum'))
    overall = 100.0 * agg['nr'] / agg['nt']
    rows.append(['<b>All (per participant)</b>', f'<b>{overall.size}</b>',
                 f'<b>{overall.mean():.1f}%</b>', f'<b>{overall.std():.1f}</b>',
                 f'<b>{np.median(overall):.1f}%</b>'])
    th = ''.join(f'<th style="padding:3px 8px;border:1px solid #ccc;">{h}</th>' for h in header)
    body = ''
    for r in rows:
        body += '<tr>' + ''.join(
            f'<td style="padding:3px 8px;border:1px solid #ccc;text-align:center;">{c}</td>' for c in r) + '</tr>'
    return f'<h4>Rejection rate by stage (across participants)</h4>' \
           f'<table style="border-collapse:collapse;font-size:.85em;"><tr>{th}</tr>{body}</table>'


def build_database_report(out_dir, alldf, summ, custom_stages):
    """Database-level aggregate report from the pooled per-epoch table + per-file/stage summary."""
    import mne
    stage_order = L.ordered_present_stages(alldf['stage'].values, custom_stages)
    thr = _default_thresholds()
    metrics = {k: alldf[k].values for k in ('ptp', 'gradient', 'mae', 'r2')}

    report = mne.Report(title='QC of rejected epochs — database-level', verbose=False)
    report.add_html(html=_db_rejection_table_html(summ, stage_order), title='Rejection summary',
                    section='Overview')

    fig = _plot_per_participant_bar(summ)
    report.add_figure(fig=fig, title='Rejection rate per participant', section='Overview',
                      image_format='PNG')
    plt.close(fig)

    fig = _plot_per_stage_box(summ, stage_order)
    report.add_figure(fig=fig, title='Rejection rate by stage', section='Overview', image_format='PNG')
    plt.close(fig)

    fig = L.plot_metric_distributions(metrics, alldf['stage'].values, alldf['reject_flag'].values,
                                      thr, stage_order, title_prefix='DB pooled — ')
    report.add_figure(fig=fig, title='Pooled metric distributions (all epochs, all participants)',
                      section='Pooled', image_format='PNG')
    plt.close(fig)

    fig = L.plot_metric_scatter(metrics, alldf['reject_flag'].values,
                                alldf['reject_method'].values, thr, title_prefix='DB pooled — ')
    report.add_figure(fig=fig, title='Pooled p-p vs gradient scatter', section='Pooled',
                      image_format='PNG')
    plt.close(fig)

    out = out_dir / 'qc2b_database_report.html'
    report.save(str(out), overwrite=True, open_browser=False, verbose=False)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('deriv_root', help='derivatives root holding the *_all-epo.fif files')
    ap.add_argument('--no-1f', action='store_true', help='skip specparam 1/f fitting (faster)')
    ap.add_argument('--limit', type=int, default=None, help='process only the first N participants')
    ap.add_argument('--out', default=None, help='output folder (default: <deriv_root>/qc2b_reports)')
    args = ap.parse_args()

    parts = L.find_participants(args.deriv_root)
    if args.limit:
        parts = parts[:args.limit]
    if not parts:
        print(f'No *_all-epo.fif found under {args.deriv_root}')
        return
    out_dir = Path(args.out) if args.out else Path(args.deriv_root) / 'qc2b_reports'
    out_dir.mkdir(parents=True, exist_ok=True)
    custom = L.load_custom_stages(args.deriv_root)

    pooled, summary_rows, failed = [], [], []
    for i, p in enumerate(parts, 1):
        fid = p['file_id']
        try:
            print(f'[{i}/{len(parts)}] {fid} …', flush=True)
            P = L.load_participant(p['fif'])
            thr, info = L.load_params(p['folder'], fid, custom_stages_fallback=custom)
            freqs, psds = L.compute_psds(P['epochs'])
            if args.no_1f:
                ptp = np.ptp(P['data_uV'], axis=-1).max(axis=1)
                grad = np.max(np.abs(np.diff(P['data_uV'], axis=-1)), axis=-1).max(axis=1)
                metrics = {'ptp': ptp, 'gradient': grad,
                           'mae': np.full(len(ptp), np.nan), 'r2': np.full(len(ptp), np.nan)}
            else:
                metrics = L.compute_epoch_metrics(P['data_uV'], freqs, psds)
            figs, html = L.build_participant_report_figs(P, metrics, freqs, psds, thr,
                                                         info['custom_stages'])
            L.save_report_html(out_dir / f'{fid}_qc2b_report.html',
                               f'{fid} — QC of rejected epochs', figs, html)
            pooled.append(pd.DataFrame({
                'file_id': fid, 'stage': P['stages'], 'reject_flag': P['reject_flag'],
                'reject_method': P['meta']['reject_method'].astype(str).values,
                'ptp': metrics['ptp'], 'gradient': metrics['gradient'],
                'mae': metrics['mae'], 'r2': metrics['r2'],
            }))
            for st in L.ordered_present_stages(P['stages'], info['custom_stages']):
                sel = (P['stages'] == st)
                n = int(sel.sum())
                if n == 0:
                    continue
                row = {'file_id': fid, 'stage': st, 'n_total': n,
                       'n_rejected': int(P['reject_flag'][sel].sum())}
                for m in P['methods_present']:
                    row['n_' + m] = int(P['meta'].loc[sel, 'flag_' + m].sum())
                summary_rows.append(row)
        except Exception as e:
            print(f'  FAILED: {e}')
            failed.append({'file_id': fid, 'reason': str(e)})

    if not pooled:
        print('No participants processed successfully.')
        if failed:
            pd.DataFrame(failed).to_csv(out_dir / 'qc2b_failed.tsv', sep='\t', index=False)
        return

    alldf = pd.concat(pooled, ignore_index=True)
    summ = pd.DataFrame(summary_rows)
    summ['pct_rejected'] = 100.0 * summ['n_rejected'] / summ['n_total']
    summ.to_csv(out_dir / 'qc2b_database_rejection_summary.tsv', sep='\t', index=False)
    db_report = build_database_report(out_dir, alldf, summ, custom)
    if failed:
        pd.DataFrame(failed).to_csv(out_dir / 'qc2b_failed.tsv', sep='\t', index=False)

    print(f'\nDone. {len(parts) - len(failed)} processed, {len(failed)} failed.')
    print(f'Per-participant reports + {db_report.name} in {out_dir}')


if __name__ == '__main__':
    main()
