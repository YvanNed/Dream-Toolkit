"""Known-answer tests on the synthetic nights: the truth is recomputed from what tests/synthetic.py generated,
never read from a golden snapshot. A golden comparison only proves that nothing CHANGED (a bug present when the
golden was taken would pass forever); these tests prove the tools find what is really in the data.

    tool 5   sleep metrics recomputed from the generated hypnogram (exact)
    tool 6   flags exactly the two injected defective channels
    tool 7   flags the 10 injected movement bursts, and nothing else on that channel
    tool 9   the generated physiology: alpha in wake (occipital), spindles in N2, steeper 1/f in N3
"""
import sys
from pathlib import Path

import mne
import numpy as np
import pandas as pd
import pytest

import minidb
import synthetic

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'tools'))
import generate_test_data as gtd     # noqa: E402

EP_MIN = synthetic.EPOCH_S / 60.0


@pytest.fixture
def syn(chain):
    if chain['dataset'] != 'synthetic':
        pytest.skip('known answers exist for the synthetic nights only (the real ones are unscored truth)')
    return chain['data'], minidb.PARTICIPANTS['synthetic']


def tsv(path):
    return pd.read_csv(path, sep='\t', dtype={'file_id': str})


def hypnogram(data, fid):
    return np.loadtxt(data / f'{fid}_Hypnogram_remapped.txt', dtype=str)


# ----------------------------------------------------------------------------------------------- tool 5
@pytest.mark.tool5
def test_tool5_sleep_metrics_equal_the_generated_hypnogram(syn):
    """No lights source in the synthetic set: lights-off/on = the recording bounds. Sleep onset = first
    non-W epoch (the default rule); the metrics are counted between onset and the last sleep epoch."""
    data, ids = syn
    table = tsv(data / 'reports_features_macrostructure' / 'global_sleep_metrics.tsv').set_index('file_id')
    for fid in ids:
        st = hypnogram(data, fid)
        sleep = np.flatnonzero(st != 'W')
        onset, last = sleep[0], sleep[-1]
        seg = st[onset:last + 1]
        tst = (seg != 'W').sum() * EP_MIN
        expected = {'tib_min': len(st) * EP_MIN, 'sol_min': onset * EP_MIN,
                    'spt_min': (last - onset + 1) * EP_MIN, 'tst_min': tst,
                    'waso_min': (seg == 'W').sum() * EP_MIN, 'se_pct': 100 * tst / (len(st) * EP_MIN)}
        for stage, key in [('N1', 'n1'), ('N2', 'n2'), ('N3', 'n3'), ('R', 'rem')]:
            expected[f'{key}_min'] = (seg == stage).sum() * EP_MIN
            expected[f'{key}_pct'] = 100 * (seg == stage).sum() * EP_MIN / tst
            expected[f'lat_{key}_min'] = np.flatnonzero(seg == stage)[0] * EP_MIN   # from sleep onset
        for key, value in expected.items():
            assert table.loc[fid, key] == pytest.approx(value, abs=1e-3), (fid, key)


# ----------------------------------------------------------------------------------------------- tool 6
@pytest.mark.tool6
def test_tool6_flags_exactly_the_injected_defects(syn):
    data, ids = syn
    q = tsv(data / 'reports_quality_overview' / 'quality_summary.tsv')
    flagged = sorted(map(tuple, q.loc[q['exclude'].astype(bool), ['file_id', 'channel']].values.tolist()))
    assert flagged == [(ids[0], 'Fp1'), (ids[1], 'C3')]          # clipped Fp1, dead C3, nothing else


# ----------------------------------------------------------------------------------------------- tool 7
@pytest.mark.tool7
def test_tool7_flags_the_injected_movement_bursts_and_nothing_else(syn):
    """The bursts are re-created with the generator's own function and seed on a silent signal, which gives
    their exact samples: every epoch holding >= 1 s of burst must be flagged on Fp1 (amplitude or gradient),
    and no other Fp1 epoch (the synthetic background never reaches those thresholds)."""
    data, ids = syn
    fid = ids[3]
    n = synthetic.N_EPOCHS * synthetic.EPOCH_S * int(synthetic.SFREQ)
    silent = mne.io.RawArray(np.zeros((len(synthetic.CHANNELS), n)),
                             mne.create_info(synthetic.CHANNELS, synthetic.SFREQ, 'eeg'), verbose=False)
    bursts = gtd.inject_movement_bursts(silent, 'Fp1', n_bursts=10, amp_uv=300.0).get_data(picks=['Fp1'])[0]
    seconds = (np.abs(bursts) > 0).reshape(synthetic.N_EPOCHS, -1).sum(axis=1) / synthetic.SFREQ
    burst_epochs = set(np.flatnonzero(seconds >= 1.0).tolist())
    assert len(burst_epochs) == 10

    flags = tsv(data / 'reports_preprocessing' / f'{fid}_epoch_channel_rejection.tsv')
    fp1 = flags[flags['channel'] == 'Fp1'].set_index('epoch_idx')
    hit = fp1['flag_amplitude'].astype(bool) | fp1['flag_gradient'].astype(bool)
    assert set(fp1.index[hit].tolist()) == burst_epochs


# ----------------------------------------------------------------------------------------------- tool 9
def per_stage(data, fid, table, value):
    """Median of a per-epoch measure per (stage, channel[, band]) over the epochs that were not rejected and the
    channels not dropped. The per-epoch tables hold EVERY epoch, wake included, which the stage aggregates do
    not (tool 8bis keeps N2 / N3 / REM by default)."""
    t = tsv(data / 'derivatives' / 'features_spectral' / f'{fid}_{table}.tsv')
    t = t[~t['rejected'].astype(bool) & ~t['channel_dropped'].astype(bool)]
    keys = ['band', 'stage', 'channel'] if 'band' in t.columns else ['stage', 'channel']
    return t.groupby(keys)[value].median()


@pytest.mark.tool9
def test_tool9_finds_the_alpha_of_wake_on_the_occipital_channel(syn):
    data, ids = syn
    for fid in ids:
        p = per_stage(data, fid, 'bandpower_epoch', 'power_db')
        assert p[('alpha', 'W', 'O1')] - p[('alpha', 'N2', 'O1')] > 6.0, fid           # measured ~15 dB
        if ('alpha', 'W', 'Fp1') in p.index:                                            # Fp1 dropped in sim01
            assert p[('alpha', 'W', 'O1')] - p[('alpha', 'W', 'Fp1')] > 5.0, fid       # measured ~11 dB


@pytest.mark.tool9
def test_tool9_finds_the_spindle_peak_of_n2(syn):
    """Spindles (12-14 Hz) are generated strongest on C3: in N2 the mean PSD there stands above both its
    lower (9-11 Hz) and upper (16-20 Hz) neighbours. C3 is the dead channel in sim02, skipped there."""
    data, ids = syn
    checked = 0
    for fid in ids:
        ps = tsv(data / 'derivatives' / 'features_spectral' / f'{fid}_psd_stage.tsv')
        c3 = ps[(ps['stage'] == 'N2') & (ps['channel'] == 'C3') & ~ps['channel_dropped'].astype(bool)]
        if not len(c3):
            continue

        def band(lo, hi):
            return c3.loc[(c3['freq_hz'] >= lo) & (c3['freq_hz'] < hi), 'psd_db_from_log'].mean()

        assert band(12, 14) - band(9, 11) > 1.0, fid                                    # measured ~3.7 dB
        assert band(12, 14) - band(16, 20) > 3.0, fid                                   # measured ~10 dB
        checked += 1
    assert checked == 3


@pytest.mark.tool9
def test_tool9_finds_a_steeper_aperiodic_slope_in_n3(syn):
    """The generated N3 carries a steeper (1/f^3) slow component: its aperiodic exponent exceeds wake's."""
    data, ids = syn
    for fid in ids:
        e = per_stage(data, fid, 'aperiodic_epoch', 'exponent').groupby(level='stage').median()
        assert e['N3'] - e['W'] > 0.05, fid                                             # measured ~0.13
