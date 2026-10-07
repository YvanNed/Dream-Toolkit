"""Step 3 of the exclusion refactor: tools 7, 8bis and 8.

Tool 7  - deselected channels are kept in the .fif (marked bad), re-referenced like the others, never
          flagged; manual exclusion; exclusion columns + a per-stage summary over the non-excluded only.
          Re-referencing is checked EXACTLY through a property of linear processing: any common
          reference cancels out in a bipolar derivation (Fp1 - C3), and filtering is linear.
Tool 8bis - epoch / channel decision tables, equal to the former clean .fif (golden); automatic exclusions in the
          registry (and lifted again by a re-run).
Tool 8  - decision tables; the Exclude / Remove button.
"""
import json
import shutil
import sys
from pathlib import Path

import mne
import numpy as np
import pandas as pd
import pytest

import minidb
import pipeline
import snapshot
from nbdriver import NotebookSession
from pipeline import select_subset

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'tools'))
import participant_selection_lib as P     # noqa: E402


@pytest.fixture
def data_copy(chain, tmp_path):
    dst = tmp_path / 'data'
    shutil.copytree(chain['data'], dst)
    return dst


def ids(chain):
    return minidb.PARTICIPANTS[chain['dataset']]


def fif(data, fid, kind='raw_epo', suffix='all-epo'):
    return mne.read_epochs(str(Path(data) / 'derivatives' / kind / f'{fid}_{suffix}.fif'),
                           preload=True, verbose=False)


def tsv(path):
    return pd.read_csv(path, sep='\t', dtype={'file_id': str})


# =====================================================================================  tool 7
@pytest.mark.tool7
def test_tool7_keeps_deselected_channels_as_bad(chain):
    """The chain's tool-6 exclusions (clipped Fp1, dead C3) are deselected in tool 7: they are in the
    .fif, marked bad, listed in the sidecar, absent from the flagging tables."""
    p = ids(chain)
    for fid, bad in [(p[0], 'Fp1'), (p[1], 'C3')]:
        ep = fif(chain['data'], fid)
        assert ep.info['bads'] == [bad] and bad in ep.ch_names
        side = json.loads((chain['data'] / 'derivatives' / 'raw_epo'
                           / f'{fid}_preprocessing_params.json').read_text(encoding='utf-8'))
        assert side['bad_channels'] == [bad] and bad not in side['channels']
        flags = tsv(chain['data'] / 'reports_preprocessing' / f'{fid}_epoch_channel_rejection.tsv')
        assert bad not in set(flags['channel'])
    assert fif(chain['data'], p[2]).info['bads'] == []


@pytest.mark.tool7
def test_tool7_bad_channel_gets_the_m2_reference(chain, data_copy):
    """Re-run the clipping participant with Fp1 KEPT: every channel, Fp1 included, must equal the chain's
    run where Fp1 was deselected. Hence the bad Fp1 was re-referenced to M2 and filtered exactly like a
    good channel, and the good channels did not depend on it."""
    fid = ids(chain)[0]
    pipeline.run_tool7(data_copy, subset=[fid], skip=False, channels={fid: {'Fp1': True}})
    kept, deselected = fif(data_copy, fid), fif(chain['data'], fid)
    assert kept.info['bads'] == [] and deselected.info['bads'] == ['Fp1']
    for ch in ['Fp1', 'C3', 'O1']:
        np.testing.assert_allclose(kept.get_data(picks=[ch]), deselected.get_data(picks=[ch]),
                                   rtol=0, atol=1e-9, err_msg=ch)


@pytest.mark.tool7
def test_tool7_average_reference_ignores_the_bad_channel(chain, data_copy):
    """Average reference (high-density montages): MNE averages the GOOD channels only, and leaves the
    bad ones in their original reference, which tool 7 corrects. Checked on the clipping participant
    with Fp1 deselected (from tool 6) then kept."""
    fid = ids(chain)[0]
    cfg_path = data_copy / 'config_param' / 'remap_reref_persubject.json'
    cfg = json.loads(cfg_path.read_text(encoding='utf-8'))
    cfg[fid]['ref_channels'] = 'average'
    cfg_path.write_text(json.dumps(cfg, indent=2), encoding='utf-8')

    pipeline.run_tool7(data_copy, subset=[fid], skip=False)                 # Fp1 deselected
    bad_run = fif(data_copy, fid)
    pipeline.run_tool7(data_copy, subset=[fid], skip=False, channels={fid: {'Fp1': True}})
    good_run = fif(data_copy, fid)

    assert bad_run.info['bads'] == ['Fp1'] and 'M2' in bad_run.ch_names   # average ref keeps M2
    goods = [c for c in bad_run.ch_names if c != 'Fp1']
    # 1. the reference = the mean of the good channels only: they average to zero, Fp1 not included
    np.testing.assert_allclose(bad_run.get_data(picks=goods).mean(axis=1), 0, atol=1e-9)
    assert np.abs(bad_run.get_data().mean(axis=1)).max() > 1e-7
    np.testing.assert_allclose(good_run.get_data().mean(axis=1), 0, atol=1e-9)
    # 2. the bad Fp1 carries the same reference as the good channels: the bipolar Fp1 - C3 is the same
    #    whether Fp1 entered the average or not (a common reference cancels in a difference)
    bip = lambda ep: ep.get_data(picks=['Fp1']) - ep.get_data(picks=['C3'])
    np.testing.assert_allclose(bip(bad_run), bip(good_run), rtol=0, atol=1e-9)


@pytest.mark.tool7
def test_tool7_manual_exclusion_and_global_tables(chain, data_copy):
    p = ids(chain)

    def exclude_p2(nb):
        nb.set('dd_edit', p[2])
        nb.set('cb_edit_exclude', True)
        nb.click('btn_edit_save')
        assert 'Give a reason' in nb.value('lbl_edit.value')            # a reason is required
        nb.set('txt_edit_reason', 'tool 6: noisy night')
        nb.click('btn_edit_save')
        assert 'excluded' in nb.value('lbl_edit.value')
        assert p[2] in nb.value('html_edits.value')

    pipeline.run_tool7(data_copy, subset=[p[2]], skip=False, before_run=exclude_p2)
    reg, _ = P.load_registry(data_copy)
    st = P.effective_status(reg)
    assert P.status_of(st, p[2]) == (True, '7_manual: tool 6: noisy night')

    # every epoch row kept, the excluded participants flagged (p1 by 5bis, p2 here)
    g = tsv(data_copy / 'reports_preprocessing' / 'global_epoch_rejection.tsv')
    assert set(g['file_id']) == set(p)
    flagged = g.groupby('file_id')[P.COL_EXCLUDED].first()
    assert flagged.to_dict() == {p[0]: False, p[1]: True, p[2]: True, p[3]: False}

    # the per-stage summary is pooled over the 2 non-excluded participants only
    stage = tsv(data_copy / 'reports_preprocessing' / 'global_rejection_by_stage.tsv')
    assert (stage['n_participants'] == 2).all() and (stage['n_participants_excluded'] == 2).all()
    per_file = pd.concat([tsv(data_copy / 'reports_preprocessing' / f'{f}_rejection_summary.tsv')
                          for f in (p[0], p[3])])
    m0 = per_file['method'].iloc[0]
    expected = per_file[per_file['method'] == m0].groupby('stage')['n_total'].sum()
    assert dict(zip(stage['stage'], stage['n_total'])) == expected.to_dict()


@pytest.mark.tool7
def test_tool7_chain_stage_summary_leaves_out_the_5bis_exclusion(chain):
    p = ids(chain)
    stage = tsv(chain['data'] / 'reports_preprocessing' / 'global_rejection_by_stage.tsv')
    assert (stage['n_participants'] == 3).all() and (stage['n_participants_excluded'] == 1).all()
    per_file = pd.concat([tsv(chain['data'] / 'reports_preprocessing' / f'{f}_rejection_summary.tsv')
                          for f in (p[0], p[2], p[3])])
    m0 = per_file['method'].iloc[0]
    expected = per_file[per_file['method'] == m0].groupby('stage')['n_total'].sum()
    assert dict(zip(stage['stage'], stage['n_total'])) == expected.to_dict()


# =====================================================================================  tool 8bis
def golden_clean(golden, fid, kind):
    """The clean-epo the tools wrote BEFORE the refactor (golden snapshot): kept epoch indices + channels."""
    base = golden / 'derivatives' / kind / f'{fid}_clean-epo.fif'
    meta = snapshot._read_table(str(base) + '.metadata.tsv.gz')
    summary = json.loads(Path(str(base) + '.summary.json').read_text(encoding='utf-8'))
    return meta['epoch_idx'].astype(int).tolist(), summary['ch_names']


@pytest.mark.tool8bis
def test_tool8bis_decision_tables_match_the_clean_epochs(chain, golden):
    """INTENDED CHANGE (decisions cover every epoch): on the reference stages the decision tables keep
    exactly the epochs and channels the former clean-epo held (golden), and their `_ref` counts are the
    former in-scope counts; every other epoch is now decided too, by the same rule (rejected when more
    than 20 % of the good channels flag it, from tool 7's per-(epoch, channel) flags)."""
    p = ids(chain)
    dec_dir = chain['data'] / 'derivatives' / 'rejection_auto'
    rep_dir = chain['data'] / 'reports_rejection_auto'
    assert not (chain['data'] / 'derivatives' / 'clean_epo_auto').exists()     # no copy of the epochs
    for fid in p:
        ed = tsv(dec_dir / f'{fid}_epoch_decision.tsv')
        kept_old, chans_old = golden_clean(golden, fid, 'clean_epo_auto')
        rec = tsv(rep_dir / f'{fid}_autoreject_decision.tsv').iloc[0]
        ref = rec['stages_of_interest'].split('+')
        is_ref = ed['stage'].astype(str).isin(ref)
        assert len(ed) == len(fif(chain['data'], fid))                   # every epoch has a row
        assert ed['in_scope'].all()                                      # the decision covers every epoch
        assert ed[is_ref & ~ed['rejected']]['epoch_idx'].tolist() == kept_old
        assert (ed.loc[ed['rejected'], 'reject_reason'] != '').all()
        cd = tsv(dec_dir / f'{fid}_channel_decision.tsv')
        good = sorted(cd.loc[~cd['dropped'], 'channel'])
        assert good == sorted(chans_old)
        # the other stages: same epoch rule, on the good channels
        pairs = tsv(next((chain['data'] / 'reports_preprocessing').rglob(f'{fid}_epoch_channel_rejection.tsv')))
        frac = (pairs[pairs['channel'].isin(good)].groupby('epoch_idx')['reject_any']
                .mean().reindex(ed['epoch_idx']).fillna(0.0).values)
        expected = (frac > 0.20) if good else np.ones(len(ed), dtype=bool)
        assert (ed['rejected'].values == expected).all()
        assert (~is_ref).any()                                           # W/N1 epochs were decided
        # the reference-stage counts are the former in-scope ones; the every-epoch counts cover the night
        old = snapshot._read_table(str(golden / 'reports_rejection_auto' / f'{fid}_autoreject_decision.tsv.gz')).iloc[0]
        assert int(rec['n_epochs_ref']) == int(old['n_epochs'])
        assert int(rec['n_epochs_rejected_ref']) == int(old['n_epochs_rejected'])
        assert int(rec['n_epochs']) == len(ed) and int(rec['n_epochs_rejected']) == int(ed['rejected'].sum())
        assert sum(int(rec[f'n_{st}']) for st in ed['stage'].astype(str).unique()) == len(ed)
    for fid, bad in [(p[0], 'Fp1'), (p[1], 'C3')]:                      # tool-7 deselected channels
        cd = tsv(dec_dir / f'{fid}_channel_decision.tsv')
        assert cd.loc[cd['channel'] == bad, 'drop_reason'].item() == 'tool7_deselected'
    summary = tsv(chain['data'] / 'reports_rejection_auto' / 'global_autoreject_summary.tsv')
    assert summary.set_index('file_id')[P.COL_EXCLUDED].to_dict() == {
        p[0]: False, p[1]: True, p[2]: False, p[3]: False}
    reg, _ = P.load_registry(chain['data'])
    assert '8bis_auto' not in set(reg['source'])                         # nobody all-rejected here


@pytest.mark.tool8bis
def test_tool8bis_event_threshold_excludes_then_a_rerun_lifts_it(chain, data_copy):
    fid = ids(chain)[2]

    def run(threshold):
        with NotebookSession(pipeline.T8BIS) as nb:
            nb.pick('fc_raw', data_copy / 'derivatives' / 'raw_epo')
            nb.pick('fc_reports', data_copy / 'reports_preprocessing')
            nb.click('btn_scan')
            nb.run(f"S['event_type_rows']['arousal_spontaneous'][1].value = {threshold}")
            select_subset(nb, [fid], skip=False)
            nb.click('btn_run')
        reg, _ = P.load_registry(data_copy)
        return reg[reg['source'] == '8bis_auto']

    rows = run(1)
    assert rows['file_id'].tolist() == [fid]
    assert rows['reason'].iloc[0].endswith('arousal_spontaneous events > 1')
    assert run(0).empty                                                 # re-run without it: lifted


# =====================================================================================  tool 8
@pytest.mark.tool8
def test_tool8_decision_tables(chain, golden):
    fid = minidb.manual_review_participant(chain['dataset'])
    dec_dir = chain['data'] / 'derivatives' / 'rejection_manual'
    assert not (chain['data'] / 'derivatives' / 'clean_epo_manual').exists()
    ed = tsv(dec_dir / f'{fid}_epoch_decision.tsv')
    kept_old, chans_old = golden_clean(golden, fid, 'clean_epo_manual')
    assert ed['in_scope'].all()                       # the decision covers every epoch
    assert ed[ed['in_scope'] & ~ed['rejected']]['epoch_idx'].tolist() == kept_old
    cd = tsv(dec_dir / f'{fid}_channel_decision.tsv')
    assert sorted(cd.loc[~cd['dropped'], 'channel']) == sorted(chans_old)


@pytest.mark.tool8
def test_tool8_exclude_button(chain, data_copy):
    fid = ids(chain)[3]
    with NotebookSession(pipeline.T8) as nb:
        nb.pick('fc_data', data_copy)
        nb.pick('fc_raw', data_copy / 'derivatives' / 'raw_epo')
        nb.set('dd_part', fid)
        nb.click('btn_excl')
        assert 'Give a reason' in nb.value('lbl_excl.value')
        nb.set('txt_excl_reason', 'electrode off from 2 a.m.')
        nb.click('btn_excl')
        assert 'EXCLUDED' in nb.value('lbl_excl.value')
        assert nb.value('btn_excl.description') == 'Remove my exclusion'
        reg, _ = P.load_registry(data_copy)
        assert reg[reg['source'] == '8_manual']['file_id'].tolist() == [fid]
        nb.click('btn_excl')                                               # remove it again
        assert 'not excluded' in nb.value('lbl_excl.value')
    reg, _ = P.load_registry(data_copy)
    assert '8_manual' not in set(reg['source'])
