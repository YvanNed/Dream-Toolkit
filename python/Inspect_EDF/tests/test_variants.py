"""Full level: parameter variants on the synthetic database. No frozen reference: each variant runs on a
light copy of the database (one participant where possible) and is checked for what its option must do.

    tool 7   10 s epochs (then 8bis and 9 on them), resampling to 128 Hz, notch off, custom stage N4
    tool 6   resampling + high-pass: the dead channel stays flagged, no healthy channel gets flagged
    tool 9   multitaper PSD (same order as Welch), 'knee' aperiodic mode
    tool 8   a flagged epoch rescued / a kept epoch added in the navigator; an unticked stage stays in
             the decision (focus only, overrides kept);
             reloading clears the sections (after an unsaved-changes warning); a manual annotation
             re-read by tool 7
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
from nbdriver import NotebookSession
from pipeline import select_subset

pytestmark = pytest.mark.full
IDS = minidb.PARTICIPANTS['synthetic']
SIM_LINE, SIM_BURSTS = IDS[2], IDS[3]


def tsv(path):
    return pd.read_csv(path, sep='\t', dtype={'file_id': str})


def sidecar(data, fid):
    return json.loads((Path(data) / 'derivatives' / 'raw_epo' / f'{fid}_preprocessing_params.json')
                      .read_text(encoding='utf-8'))


@pytest.fixture
def src(chain):
    if chain['dataset'] != 'synthetic':
        pytest.skip('the variants run on the synthetic database only')
    return chain['data']


def inputs_copy(src, dst, dirs=()):
    """The chain's inputs (recordings, hypnograms, events, config_param, participants.tsv) + some output
    folders: the derivatives (~330 MB) are never copied whole."""
    dst.mkdir(parents=True)
    for f in src.iterdir():
        if f.is_file():
            shutil.copy2(f, dst / f.name)
    for d in ('config_param',) + tuple(dirs):
        shutil.copytree(src / d, dst / d)
    return dst


def run_tool7(data, fid, before_run=None):
    """Tool 7 on one participant (events on, Skip off), with `before_run(nb)` setting the variant."""
    return pipeline.run_tool7(data, subset=[fid], skip=False, before_run=before_run)


# =====================================================================================  tool 7
@pytest.mark.tool7
@pytest.mark.tool8bis
@pytest.mark.tool9
def test_tool7_10s_epochs_then_8bis_and_9(src, tmp_path):
    data = inputs_copy(src, tmp_path / 'data', ['reports_quality_overview'])
    run_tool7(data, SIM_LINE, lambda nb: nb.set('dd_epoch_len', 10))

    assert sidecar(data, SIM_LINE)['epoch_length_s'] == 10
    stages30 = np.loadtxt(data / f'{SIM_LINE}_Hypnogram_remapped.txt', dtype=str)
    rej = tsv(data / 'reports_preprocessing' / f'{SIM_LINE}_epoch_rejection.tsv')
    assert len(rej) == 3 * len(stages30)
    assert rej['stage'].tolist() == np.repeat(stages30, 3).tolist()     # each sub-epoch inherits its stage

    pipeline.run_tool8bis(data)
    dec = tsv(data / 'derivatives' / 'rejection_auto' / f'{SIM_LINE}_epoch_decision.tsv')
    assert len(dec) == 3 * len(stages30)

    pipeline.run_tool9(data, data / 'derivatives' / 'raw_epo', data / 'derivatives' / 'rejection_auto')
    bp = tsv(data / 'derivatives' / 'features_spectral' / f'{SIM_LINE}_bandpower_epoch.tsv')
    assert bp['epoch_idx'].nunique() == 3 * len(stages30)
    kept = bp[~bp['rejected'] & ~bp['channel_dropped']]
    assert np.isfinite(kept['power_db']).all()


@pytest.mark.tool7
def test_tool7_resampling(src, tmp_path):
    data = inputs_copy(src, tmp_path / 'data', ['reports_quality_overview'])

    def resample(nb):
        nb.set('cb_resample', True)
        nb.set('txt_target_freq', 128)
    run_tool7(data, SIM_LINE, resample)

    ep = mne.read_epochs(str(data / 'derivatives' / 'raw_epo' / f'{SIM_LINE}_all-epo.fif'), verbose=False)
    assert ep.info['sfreq'] == 128 and len(ep) == 480 and ep.get_data().shape[-1] == 30 * 128
    side = sidecar(data, SIM_LINE)
    assert side['sfreq_hz'] == 128 and side['resample'] == {'applied': True, 'target_freq_hz': 128.0}


def power_50hz_db(fif_path, channel):
    ep = mne.read_epochs(str(fif_path), verbose=False).pick([channel])
    psd = ep.compute_psd(method='welch', fmin=45, fmax=55, n_fft=1024, verbose=False)
    freqs, p = psd.freqs, psd.get_data().mean(axis=(0, 1))
    return 10 * np.log10(p[np.argmin(np.abs(freqs - 50))] * 1e12)


@pytest.mark.tool7
def test_tool7_notch_off_keeps_the_line_noise(src, tmp_path):
    """sim03 carries 50 Hz line noise on O1: without the notch it stays far above the chain's (notch on)."""
    data = inputs_copy(src, tmp_path / 'data', ['reports_quality_overview'])
    run_tool7(data, SIM_LINE, lambda nb: nb.set('cb_notch', False))

    assert sidecar(data, SIM_LINE)['notch'] == {'applied': False, 'freq_hz': None}
    off = power_50hz_db(data / 'derivatives' / 'raw_epo' / f'{SIM_LINE}_all-epo.fif', 'O1')
    on = power_50hz_db(src / 'derivatives' / 'raw_epo' / f'{SIM_LINE}_all-epo.fif', 'O1')
    assert off - on > 10, (off, on)


@pytest.mark.tool7
@pytest.mark.tool8bis
@pytest.mark.tool9
def test_tool7_custom_stage_n4_reaches_8bis_and_9(src, tmp_path):
    """The N3 epochs of the first half of sim03's night become 'N4', declared in custom_stages.json
    (what tool 3 writes): N4 gets its own amplitude threshold in 7 and is carried by 8bis and 9."""
    data = inputs_copy(src, tmp_path / 'data', ['reports_quality_overview'])
    hyp = data / f'{SIM_LINE}_Hypnogram_remapped.txt'
    stages = np.loadtxt(hyp, dtype=str)
    n4 = (stages == 'N3') & (np.arange(len(stages)) < len(stages) // 2)
    assert n4.sum() > 20
    stages[n4] = 'N4'
    hyp.write_text('\n'.join(stages) + '\n', encoding='utf-8')
    (data / 'config_param' / 'custom_stages.json').write_text(json.dumps({'custom_stages': ['N4']}),
                                                             encoding='utf-8')
    run_tool7(data, SIM_LINE)

    side = sidecar(data, SIM_LINE)
    assert side['custom_stages'] == ['N4'] and 'N4' in side['rejection_thresholds']['amplitude_ptp_uV']
    # The .fif must read back whole: a stage without an event code used to make MNE refuse the file.
    ep = mne.read_epochs(str(data / 'derivatives' / 'raw_epo' / f'{SIM_LINE}_all-epo.fif'), verbose=False)
    assert len(ep) == len(stages) and ep.event_id['N4'] == 5
    assert ep.metadata['stage'].tolist() == stages.tolist()
    rej = tsv(data / 'reports_preprocessing' / f'{SIM_LINE}_epoch_rejection.tsv')
    assert rej['stage'].tolist() == stages.tolist()
    summary = tsv(data / 'reports_preprocessing' / f'{SIM_LINE}_rejection_summary.tsv')
    assert 'N4' in set(summary['stage'])

    pipeline.run_tool8bis(data)
    dec = tsv(data / 'derivatives' / 'rejection_auto' / f'{SIM_LINE}_epoch_decision.tsv')
    assert (dec['stage'] == 'N4').sum() == n4.sum()

    pipeline.run_tool9(data, data / 'derivatives' / 'raw_epo', data / 'derivatives' / 'rejection_auto')
    bp = tsv(data / 'derivatives' / 'features_spectral' / f'{SIM_LINE}_bandpower_epoch.tsv')
    assert set(bp.loc[bp['stage'] == 'N4', 'epoch_idx']) == set(np.flatnonzero(n4))


# =====================================================================================  tool 6
@pytest.mark.tool6
def test_tool6_resampling_and_highpass(src, tmp_path):
    data = inputs_copy(src, tmp_path / 'data')
    with NotebookSession(pipeline.T6) as nb:
        nb.pick('fc_folder', data)
        nb.pick('fc_config', data / 'config_param' / 'remap_reref_persubject.json')
        nb.set('cb_resample', True)
        nb.set('txt_target_freq', 128)
        nb.set('hp_check', True)
        nb.set('hp_freq', 0.5)
        nb.click('btn_run')
    q = tsv(data / 'reports_quality_overview' / 'quality_summary.tsv')
    assert sorted(q['file_id'].unique()) == IDS and len(q) == 4 * 4
    flagged = sorted(map(tuple, q.loc[q['exclude'].astype(bool), ['file_id', 'channel']].values.tolist()))
    # The metrics are computed on the FILTERED signal on purpose (a high-pass can rescue a channel worth
    # keeping), so the filtering may smooth sim01's clipping plateaus below the limits: the dead C3 must
    # stay flagged, the clipped Fp1 may or may not, and nothing else may be flagged. The synthetic M2 (a
    # quiet mastoid, flat_pct 3.4 % unfiltered) can cross the 3.5 % flat limit once the slow drift is
    # high-passed away: it is left out of the comparison.
    flagged = [f for f in flagged if f[1] != 'M2']
    assert (IDS[1], 'C3') in flagged
    assert set(flagged) <= {(IDS[0], 'Fp1'), (IDS[1], 'C3')}


# =====================================================================================  tool 9
def copy_participant(src, data, fid, folders):
    """One participant's files of some derivatives/ folders (tool 9 writes its data beside the raw_epo folder
    it reads, so it must never be pointed at the chain's own)."""
    for folder in folders:
        (data / 'derivatives' / folder).mkdir(parents=True)
        for f in (src / 'derivatives' / folder).glob(f'{fid}_*'):
            shutil.copy2(f, data / 'derivatives' / folder / f.name)


def run_tool9_one(src, data, fid, set_option):
    copy_participant(src, data, fid, ['raw_epo', 'rejection_auto'])
    with NotebookSession(pipeline.T9) as nb:
        nb.pick('fc_data', data)
        nb.pick('fc_raw', data / 'derivatives' / 'raw_epo')
        nb.pick('fc_decision', data / 'derivatives' / 'rejection_auto')
        nb.click('btn_scan')
        set_option(nb)
        select_subset(nb, [fid], skip=False)
        nb.click('btn_run')
    return data / 'derivatives' / 'features_spectral'


@pytest.mark.tool9
def test_tool9_multitaper_matches_welch(src, tmp_path):
    data = inputs_copy(src, tmp_path / 'data')
    out = run_tool9_one(src, data, SIM_LINE, lambda nb: nb.set('dd_psd_method', 'multitaper'))
    keys = ['epoch_idx', 'channel', 'band']
    mt = tsv(out / f'{SIM_LINE}_bandpower_epoch.tsv')
    welch = tsv(src / 'derivatives' / 'features_spectral' / f'{SIM_LINE}_bandpower_epoch.tsv')
    kept = mt[~mt['rejected'] & ~mt['channel_dropped']]
    assert (kept['power_mean_uV2_Hz'] > 0).all() and np.isfinite(kept['power_db']).all()
    both = kept.merge(welch, on=keys, suffixes=('_mt', '_welch'))
    assert len(both) == len(kept)
    diff = (both['power_db_mt'] - both['power_db_welch']).groupby(both['band']).median()
    assert diff.abs().max() < 3.0, diff                  # same density units: a few dB at most


@pytest.mark.tool9
def test_tool9_knee_mode(src, tmp_path):
    data = inputs_copy(src, tmp_path / 'data')
    out = run_tool9_one(src, data, SIM_LINE, lambda nb: nb.set('dd_ap_mode', 'knee'))
    ap = tsv(out / f'{SIM_LINE}_aperiodic_epoch.tsv')
    kept = ap[~ap['rejected'] & ~ap['channel_dropped']]
    assert kept['knee'].notna().mean() > 0.9
    assert kept['fit_ok'].astype(bool).mean() > 0.8


# =====================================================================================  tool 8
def burst_copy(src, tmp_path):
    """Inputs + sim04's tool-7 outputs (the .fif tool 8 reviews)."""
    data = inputs_copy(src, tmp_path / 'data', ['reports_quality_overview'])
    raw = data / 'derivatives' / 'raw_epo'
    raw.mkdir(parents=True)
    for f in (src / 'derivatives' / 'raw_epo').glob(f'{SIM_BURSTS}_*'):
        shutil.copy2(f, raw / f.name)
    return data


def open_navigator(nb, data):
    nb.pick('fc_data', data)
    nb.pick('fc_raw', data / 'derivatives' / 'raw_epo')
    nb.set('dd_part', SIM_BURSTS)
    nb.click('btn_load')
    nb.click('btn_run_nav')
    return nb.value('[int(e) for e in S["nav_list"]]')


@pytest.mark.tool8
def test_tool8_rescue_in_the_navigator(src, tmp_path):
    """A flagged epoch set to Keep in the navigator is saved as rescued; the other flagged ones stay rejected."""
    data = burst_copy(src, tmp_path)
    with NotebookSession(pipeline.T8) as nb:
        nav = open_navigator(nb, data)
        in_scope = nb.value('[int(e) for e in np.flatnonzero(S["focus"])]')
        target = next(e for e in nav if e in in_scope)
        nb.set('sl_epoch', nav.index(target))
        nb.set('tgl_keep', 'keep')
        nb.click('btn_save')
    dec = tsv(data / 'derivatives' / 'rejection_manual' / f'{SIM_BURSTS}_epoch_decision.tsv').set_index('epoch_idx')
    assert dec.loc[target, 'reject_reason'] == 'manual_rescued' and not dec.loc[target, 'rejected']
    assert bool(dec.loc[target, 'manual_override'])
    others = [e for e in nav if e in in_scope and e != target]
    assert dec.loc[others, 'rejected'].all()


@pytest.mark.tool8
def test_tool8_add_in_the_navigator(src, tmp_path):
    """Show = kept lists the focus-stage epochs the selection does not reject; one set to Reject is saved as
    manual_added. The default 'flagged' list is the classic one."""
    data = burst_copy(src, tmp_path)
    with NotebookSession(pipeline.T8) as nb:
        flagged = open_navigator(nb, data)
        assert flagged == nb.value('[int(e) for e in np.flatnonzero(S["base_reject"])]')
        nb.set('dd_show', 'kept')
        kept = nb.value('[int(e) for e in S["nav_list"]]')
        assert kept and not set(kept) & set(flagged)
        target = kept[0]
        nb.set('sl_epoch', 0)
        nb.set('tgl_keep', 'reject')
        nb.click('btn_save')
    dec = tsv(data / 'derivatives' / 'rejection_manual' / f'{SIM_BURSTS}_epoch_decision.tsv').set_index('epoch_idx')
    assert dec.loc[target, 'reject_reason'] == 'manual_added' and dec.loc[target, 'rejected']
    assert bool(dec.loc[target, 'manual_override'])
    assert dec.loc[flagged, 'rejected'].all()
    log = tsv(data / 'reports_rejection_manual' / f'{SIM_BURSTS}_qc2b_review_log.tsv')
    assert log.set_index('epoch_idx').loc[target, 'action'] == 'added'


@pytest.mark.tool8
def test_tool8_unticked_stage_stays_in_the_decision(src, tmp_path):
    """A stage unticked (U = W, or N1 when no W epoch is flagged) = out of the review FOCUS only: the
    navigator no longer lists its epochs (no 'out of scope' mode any more), a manual override survives
    that Apply (the decision does not depend on the stages), the saved table still decides every U epoch
    (tool 7's flags: rejected), and the 'all' Apply-to-all needs a confirming second click. Changing a
    METHOD then does reset the overrides."""
    data = burst_copy(src, tmp_path)
    with NotebookSession(pipeline.T8) as nb:
        nav = open_navigator(nb, data)
        stages = nb.value('[str(s) for s in S["P"]["stages"]]')
        base = nb.value('[bool(b) for b in S["base_reject"]]')
        U = next(st for st in ('W', 'N1') if any(s == st and base[e] for e, s in enumerate(stages)))
        w_flagged = [e for e, s in enumerate(stages) if s == U and base[e]]
        target = next(e for e in nav if stages[e] != U)
        nb.set('sl_epoch', nav.index(target))
        nb.set('tgl_keep', 'keep')                       # an unsaved override on a non-W epoch
        nb.set(f"_stage_cb['{U}']", False)
        nb.click('btn_apply')                            # focus only: no unsaved-changes refusal
        assert nb.value('sorted(int(e) for e in S["overridden"])') == [target]
        assert nb.value('bool(S["final_reject"][%d])' % target) is False
        assert 'manual overrides kept' in nb.widget_text()
        assert 'out' not in nb.value('[v for _, v in dd_show.options]')
        nb.click('btn_run_nav')
        assert nb.value('[int(e) for e in S["nav_list"]]')
        assert nb.value(f'all(S["P"]["stages"][e] != "{U}" for e in S["nav_list"])')
        nb.set('dd_show', 'all')
        nb.set('tgl_keep', 'keep')                       # toggles the first epoch only; then bulk 'keep'
        nb.click('btn_apply_all')                        # first click only arms the confirmation
        assert nb.value('S["apply_all_armed"]') is False
        nb.click('btn_apply_all')
        assert nb.value('int(np.sum(np.asarray(S["final_reject"]) & np.asarray(S["focus"])))') == 0
        nb.click('btn_save')
        # a method change redecides: the overrides are reset (the save above cleared the unsaved flag)
        m = nb.value('S["methods_sel"][0]')
        nb.set(f"_method_cb['{m}']", False)
        nb.click('btn_apply')
        assert nb.value('len(S["overridden"])') == 0
    dec = tsv(data / 'derivatives' / 'rejection_manual' / f'{SIM_BURSTS}_epoch_decision.tsv').set_index('epoch_idx')
    assert dec['in_scope'].all()
    assert dec.loc[w_flagged, 'rejected'].all()          # never reviewed: tool 7's automatic decision
    assert not dec.loc[dec['stage'] != U, 'rejected'].any()
    rec = tsv(data / 'reports_rejection_manual' / f'{SIM_BURSTS}_manualreject_decision.tsv').iloc[0]
    assert rec['stages_used'] == '+'.join(s for s in ('W', 'N1', 'N2', 'N3', 'R') if s != U)
    assert int(rec['n_rejected_ref']) == 0
    assert int(rec['n_rejected']) == len(w_flagged) and int(rec['n_epochs']) == len(dec)


@pytest.mark.tool8
def test_tool8_reload_clears_the_sections(src, tmp_path):
    """Loading a participant with unsaved changes first only warns; the second click clears the
    navigator, so a toggle can never act on the previous participant's epoch list."""
    data = burst_copy(src, tmp_path)
    with NotebookSession(pipeline.T8) as nb:
        nav = open_navigator(nb, data)
        nb.set('tgl_keep', 'keep')
        nb.click('btn_load')                     # refused: unsaved change
        assert nb.value('[int(e) for e in S["nav_list"]]') == nav
        assert nb.value('len(S["overridden"])') == 1
        assert 'Unsaved manual changes' in nb.widget_text()
        nb.click('btn_load')                     # second click discards
        assert nb.value('S["nav_list"]') == [] and nb.value('len(S["overridden"])') == 0
        assert nb.value('S["dirty"]') is False


@pytest.mark.tool8
@pytest.mark.tool9
def test_tool9_names_a_stage_out_of_scope_in_a_legacy_decision(src, tmp_path):
    """Backward compatibility: decision tables written before tools 8/8bis covered every epoch can hold
    out-of-scope epochs. Simulated by marking sim04's N3 epochs in_scope False in its tool-8 table: tool 9
    keeps them in the per-epoch tables (in_scope False), writes no N3 row in its stage table, and says why in
    its log. Tool 9's scan only offers the stages some participant keeps, so a second participant keeping N3
    is needed (its 8bis decision tables, same schema, placed in the same decision folder)."""
    data = burst_copy(src, tmp_path)
    with NotebookSession(pipeline.T8) as nb:
        open_navigator(nb, data)
        nb.click('btn_save')
    manual = data / 'derivatives' / 'rejection_manual'
    legacy = manual / f'{SIM_BURSTS}_epoch_decision.tsv'
    ed = tsv(legacy)
    ed.loc[ed['stage'] == 'N3', ['in_scope', 'rejected']] = False
    ed.to_csv(legacy, sep='	', index=False)
    for folder, dst in (('raw_epo', data / 'derivatives' / 'raw_epo'), ('rejection_auto', manual)):
        for f in (src / 'derivatives' / folder).glob(f'{SIM_LINE}_*'):
            if folder == 'raw_epo' or f.name.endswith('_decision.tsv'):
                shutil.copy2(f, dst / f.name)
    with NotebookSession(pipeline.T9) as nb:
        nb.pick('fc_data', data)
        nb.pick('fc_raw', data / 'derivatives' / 'raw_epo')
        nb.pick('fc_decision', manual)
        nb.click('btn_scan')
        assert nb.value('"N3" in stage_checkboxes')
        select_subset(nb, [SIM_BURSTS, SIM_LINE], skip=False)
        nb.click('btn_run')
        text = nb.widget_text()
    assert 'stages requested but not analysed' in text and 'N3 (all ' in text and 'out of scope' in text
    out = data / 'derivatives' / 'features_spectral'
    assert 'N3' not in set(tsv(out / f'{SIM_BURSTS}_spectral_stage.tsv')['stage'])
    ep = tsv(out / f'{SIM_BURSTS}_bandpower_epoch.tsv')
    assert (ep['stage'] == 'N3').any() and not ep.loc[ep['stage'] == 'N3', 'in_scope'].any()


@pytest.mark.tool8
@pytest.mark.tool7
def test_tool8_annotation_is_read_back_by_tool7(src, tmp_path):
    """A 'hypopnea' annotated in tool 8 on an epoch with no scored event lands in sim04_manual_events.tsv
    beside the recording; tool 7 with 'Include tool-8 manual annotations' ticked lists it among the event
    onsets and flags its epoch."""
    data = burst_copy(src, tmp_path)
    scored = tsv(src / 'reports_preprocessing' / f'{SIM_BURSTS}_epoch_rejection.tsv').set_index('epoch_idx')
    with NotebookSession(pipeline.T8) as nb:
        nav = open_navigator(nb, data)
        target = next(e for e in nav if not scored.loc[e, 'flag_event'])
        nb.set('sl_epoch', nav.index(target))
        nb.set('dd_annot_label', 'hypopnea')
        nb.set('ft_annot_offset', 12.0)
        nb.click('btn_annot_add')
    manual = tsv(data / f'{SIM_BURSTS}_manual_events.tsv')
    onset = target * 30 + 12.0
    assert manual[['type', 'onset_s', 'epoch_idx']].values.tolist() == [['hypopnea', onset, target]]

    run_tool7(data, SIM_BURSTS, lambda nb: nb.set('cb_manual_events', True))
    onsets = tsv(data / 'derivatives' / 'raw_epo' / f'{SIM_BURSTS}_event_onsets.tsv')
    assert ((onsets['type'] == 'hypopnea') & np.isclose(onsets['onset_s'], onset)).sum() == 1
    rej = tsv(data / 'reports_preprocessing' / f'{SIM_BURSTS}_epoch_rejection.tsv').set_index('epoch_idx')
    assert rej.loc[target, 'flag_event']
    flags = tsv(data / 'reports_preprocessing' / f'{SIM_BURSTS}_event_epoch_flags.tsv').set_index('epoch_idx')
    assert flags.loc[target, 'evt_hypopnea']
