"""Tool 9 after step 4: every epoch and channel computed from tool 7's *_all-epo.fif + the decision tables.

1. The central check of the refactor: on the kept epochs of the kept channels, every value equals what the
   former route (tool 9 on the clean-epo written by 8bis) produced: the golden snapshot.
2. Every epoch / channel is in the tables, with the decision columns; excluded participants flagged.
3. 'Re-aggregate only' after a decision change equals a full recomputation with that decision.
"""
import json
import shutil
from pathlib import Path

import mne
import numpy as np
import pandas as pd
import pytest

import minidb
import pipeline
import snapshot

EPOCH_TABLES = ['bandpower_epoch', 'aperiodic_epoch', 'periodic_peaks']
STAGE_TABLES = ['spectral_stage', 'aperiodic_stage', 'psd_stage']
SORT_KEYS = ['file_id', 'stage', 'band', 'third', 'channel', 'freq_hz', 'epoch_idx', 'peak_idx']


def tsv(path):
    return pd.read_csv(path, sep='\t', dtype={'file_id': str})


def sort(df):
    keys = [k for k in SORT_KEYS if k in df.columns]
    return df.sort_values(keys, kind='stable').reset_index(drop=True)


def same(new, old, ignore=('source',), atol=1e-9):
    cols = [c for c in old.columns if c not in ignore]
    return snapshot.compare_tables(new[cols].reset_index(drop=True), old[cols].reset_index(drop=True),
                                   atol=atol)


def kept_rows(df):
    """The rows the former clean-epo route had: kept epochs (in scope, not rejected) of kept channels."""
    keep = ~df['channel_dropped'].fillna(False).astype(bool)
    if 'in_scope' in df.columns:
        keep &= df['in_scope'].astype(bool) & ~df['rejected'].astype(bool)
    return df[keep]


def test_tool9_reproduces_the_clean_epochs_route(chain, golden):
    data, feat = chain['data'], 'derivatives/features_spectral'
    for fid in minidb.PARTICIPANTS[chain['dataset']]:
        for name in EPOCH_TABLES + STAGE_TABLES:
            old = snapshot._read_table(golden / feat / f'{fid}_{name}.tsv.gz')
            new = tsv(data / feat / f'{fid}_{name}.tsv')
            # 'third' changed on purpose: the night thirds now span the WHOLE night
            diffs = same(kept_rows(new), old, ignore=('source', 'third'))
            assert not diffs, f'{fid}_{name}: ' + '; '.join(diffs)
    for name in ['global_spectral_stage', 'global_aperiodic_stage', 'global_psd_stage']:
        old = snapshot._read_table(golden / 'reports_features_spectral' / f'{name}.tsv.gz')
        new = tsv(data / 'reports_features_spectral' / f'{name}.tsv')
        if 'n_epochs' in old.columns and name != 'global_psd_stage':
            old = old[old['n_epochs'] > 0]           # the former NaN padding of the dropped channels
        diffs = same(sort(kept_rows(new)), sort(old))
        assert not diffs, f'{name}: ' + '; '.join(diffs)


def test_tool9_computes_every_epoch_and_channel(chain):
    data = chain['data']
    p = minidb.PARTICIPANTS[chain['dataset']]
    feat, dec = data / 'derivatives' / 'features_spectral', data / 'derivatives' / 'rejection_auto'
    for fid in p:
        ep = mne.read_epochs(str(data / 'derivatives' / 'raw_epo' / f'{fid}_all-epo.fif'),
                             preload=False, verbose=False)
        ed = tsv(dec / f'{fid}_epoch_decision.tsv').set_index('epoch_idx')
        cd = tsv(dec / f'{fid}_channel_decision.tsv')
        dropped = set(cd.loc[cd['dropped'].astype(bool), 'channel']) | set(ep.info['bads'])
        ap = tsv(feat / f'{fid}_aperiodic_epoch.tsv')
        assert len(ap) == len(ep) * len(ep.ch_names)                      # every epoch x channel
        assert set(ap['channel']) == set(ep.ch_names)
        by_epoch = ap.drop_duplicates('epoch_idx').set_index('epoch_idx')
        assert (by_epoch['in_scope'] == ed['in_scope'].reindex(by_epoch.index)).all()
        assert (by_epoch['rejected'] == ed['rejected'].reindex(by_epoch.index)).all()
        assert set(ap.loc[ap['channel_dropped'], 'channel']) == dropped
        if dropped:                                                        # dropped: still computed
            assert ap.loc[ap['channel'].isin(dropped), 'exponent'].notna().any()
        z = np.load(feat / f'{fid}_psd_epoch.npz')
        assert z['psd'].shape[:2] == (len(ep), len(ep.ch_names))
        st = tsv(feat / f'{fid}_spectral_stage.tsv')
        assert {'channel_dropped', 'n_epochs_rejected', 'n_epochs_out_of_scope'} <= set(st.columns)
    g = tsv(data / 'reports_features_spectral' / 'global_spectral_stage.tsv')
    assert g.groupby('file_id')['participant_excluded'].first().to_dict() == {
        p[0]: False, p[1]: True, p[2]: False, p[3]: False}


@pytest.fixture
def data_copy(chain, tmp_path):
    dst = tmp_path / 'data'
    shutil.copytree(chain['data'], dst)
    return dst


def test_tool9_reaggregate_equals_a_full_run(chain, data_copy):
    """Reject 25 more epochs in the 8bis decision of one participant: 'Re-aggregate only' (no PSD, no
    fit) must give the same tables as recomputing everything with that decision."""
    fid = minidb.PARTICIPANTS[chain['dataset']][0]
    feat = data_copy / 'derivatives' / 'features_spectral'
    dec_path = data_copy / 'derivatives' / 'rejection_auto' / f'{fid}_epoch_decision.tsv'
    ed = tsv(dec_path)
    newly = ed.index[ed['in_scope'] & ~ed['rejected']][:25]
    ed.loc[newly, 'rejected'] = True
    ed.loc[newly, 'reject_reason'] = 'test: rejected by hand'
    ed.to_csv(dec_path, sep='\t', index=False)
    before = tsv(feat / f'{fid}_spectral_stage.tsv')

    raw_dir, dec_dir = data_copy / 'derivatives' / 'raw_epo', data_copy / 'derivatives' / 'rejection_auto'
    log = pipeline.run_tool9(data_copy, raw_dir, dec_dir, subset=[fid], reaggregate=True)
    assert 'Re-aggregate only' in log and '1 processed' in log and 'failed:' not in log, log[-2000:]
    reagg = {n: tsv(feat / f'{fid}_{n}.tsv') for n in EPOCH_TABLES + STAGE_TABLES + ['spectral_stage_third']}
    assert reagg['spectral_stage']['n_epochs'].sum() < before['n_epochs'].sum()   # the change took effect
    assert reagg['spectral_stage']['n_epochs_rejected'].sum() > before['n_epochs_rejected'].sum()

    log = pipeline.run_tool9(data_copy, raw_dir, dec_dir, subset=[fid], skip=False)  # full recomputation
    assert '1 processed' in log and 'failed:' not in log, log[-2000:]
    for name, a in reagg.items():
        b = tsv(feat / f'{fid}_{name}.tsv')
        # the per-epoch spectrum is stored in float32: psd_stage agrees to ~1e-7 relative
        diffs = same(a, b, ignore=(), atol=1e-5)
        assert not diffs, f'{name}: ' + '; '.join(diffs)


def test_tool9_night_thirds_count_their_own_epochs(chain):
    """n_epochs_rejected / n_epochs_out_of_scope / pct_rejected of a stage x third cell are THAT cell's,
    recounted here from the per-epoch table (which carries `third`, `in_scope`, `rejected`)."""
    feat = chain['data'] / 'derivatives' / 'features_spectral'
    for fid in minidb.PARTICIPANTS[chain['dataset']]:
        th = tsv(feat / f'{fid}_spectral_stage_third.tsv')
        ep = tsv(feat / f'{fid}_aperiodic_epoch.tsv').drop_duplicates('epoch_idx')
        for (stage, third), row in th.drop_duplicates(['stage', 'third']).set_index(['stage', 'third']).iterrows():
            cell = ep[(ep['stage'] == stage) & (ep['third'] == third)]
            n_rej = int((cell['in_scope'] & cell['rejected']).sum())
            n_kept = int((cell['in_scope'] & ~cell['rejected']).sum())
            assert row['n_epochs_rejected'] == n_rej, (fid, stage, third)
            assert row['n_epochs_out_of_scope'] == int((~cell['in_scope']).sum())
            assert row['pct_rejected'] == pytest.approx(100.0 * n_rej / (n_kept + n_rej), abs=0.01)
