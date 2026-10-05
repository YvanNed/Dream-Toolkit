"""Full level: parameter variants on the synthetic database. No frozen reference: each variant runs on a
light copy of the database (one participant where possible) and is checked for what its option must do.

    tool 7   10 s epochs (then 8bis and 9 on them), resampling to 128 Hz, notch off, custom stage N4
    tool 6   resampling + high-pass: the dead channel stays flagged, no healthy channel gets flagged
    tool 9   multitaper PSD (same order as Welch), 'knee' aperiodic mode
    tool 8   a flagged epoch rescued in the navigator; a manual annotation re-read by tool 7
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
    """A flagged epoch set to Keep in the navigator is saved as rescued; the other flagged ones stay rejected.
    (The navigator lists flagged epochs only, by design, so 'manual_added' has no path in the interface.)"""
    data = burst_copy(src, tmp_path)
    with NotebookSession(pipeline.T8) as nb:
        nav = open_navigator(nb, data)
        in_scope = nb.value('[int(e) for e in np.flatnonzero(S["in_scope"])]')
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
