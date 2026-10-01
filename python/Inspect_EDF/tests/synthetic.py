"""Synthetic sleep recordings for the tests: no real participant's data, fully reproducible.

Each night is generated from a fixed seed: a synthetic hypnogram (sleep cycles of ~90 min, N3 early in
the night, REM late), then an EEG whose content follows the stages, so every tool sees a plausible
signal (the 1/f fit converges, stage-wise PSDs differ, spindles show in N2):

    background   power-law noise (1/f^2), ~10 uV, on every channel
    W            alpha ~10 Hz (strongest on O1) + some fast activity
    N1, R        theta ~5.5 Hz
    N2           sleep spindles (12-14 Hz, ~1 s bursts) + some slow activity
    N3           steeper 1/f^3 slow activity (~25 uV sd, i.e. 100-200 uV peak-to-peak)

Every rhythm has a SMOOTH spectral shape: a brick-wall band would put a step in the PSD that no 1/f
model fits, and tool 7's 1/f-quality flag would then reject every epoch.

One defect is then injected per recording with the very functions that built tools/test_data
(tools/generate_test_data.py), so the synthetic set mirrors the real one. Scored events (arousals,
hypopneas, apneas, desaturations) are drawn at random inside sleep, with the Compumedics raw labels,
so tool 4's event_remap.json maps them.
"""
import sys
from pathlib import Path

import edfio
import mne
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / 'tools'))
import generate_test_data as gtd     # noqa: E402  (defect injectors shared with the real fixtures)

SFREQ = 256.0
EPOCH_S = 30
N_EPOCHS = 480                       # 4 h: long enough for 2-3 sleep cycles, short enough to run fast
CHANNELS = ['Fp1', 'C3', 'A2', 'O1']
PHYS_RANGE_UV = (-500.0, 500.0)

# file_id -> (seed, defect builder, physical-range overrides): the defects of the real 731/734/736/737
DEFECTS = {
    # clipped at 50 uV (75 in the real 731): the synthetic EEG is a bit smaller, and 50 uV puts the
    # same ~1 % of samples on the bounds, enough for tool 6 to flag the channel
    'sim01_clipping': (101, lambda r: gtd.inject_clipping(r, 'Fp1', threshold_uv=50.0),
                       {'Fp1': (-50.0, 50.0)}),
    'sim02_dead_channel': (102, lambda r: gtd.inject_dead_channel(r, 'C3', scale=0.01), None),
    'sim03_line_noise': (103, lambda r: gtd.inject_line_noise(r, 'O1', amp_uv_pp=50.0), None),
    'sim04_movement_bursts': (104, lambda r: gtd.inject_movement_bursts(r, 'Fp1', n_bursts=10,
                                                                        amp_uv=300.0), None),
}

# raw label -> share of the events (labels as exported by Compumedics, see event_remap.json)
EVENT_LABELS = [('Arousal (ARO SPONT)', 30), ('Arousal (ARO RES)', 6), ('Hypopnea', 12),
                ('Obstructive Apnea', 3), ('SpO2 desaturation', 4)]


def make_hypnogram(rng):
    """W at lights-off, then ~90 min cycles (N1 -> N2 -> N3 -> N2 -> R) until N_EPOCHS, short W at the
    end. N3 shrinks and REM grows across the night, as in a healthy adult; durations jitter per seed."""
    stages = ['W'] * int(rng.integers(12, 40))                       # sleep latency 6-20 min
    cycle = 0
    while len(stages) < N_EPOCHS - 10:
        n3 = max(0, 50 - 20 * cycle + int(rng.integers(-5, 6)))
        rem = 15 + 12 * cycle + int(rng.integers(-3, 4))
        stages += (['N1'] * int(rng.integers(4, 10)) + ['N2'] * int(rng.integers(30, 50)) + ['N3'] * n3
                   + ['N2'] * int(rng.integers(15, 30)) + ['R'] * rem)
        if rng.random() < 0.5:
            stages += ['W'] * int(rng.integers(1, 6))                 # brief awakening
        cycle += 1
    stages = stages[:N_EPOCHS - 10] + ['W'] * 10
    return np.array(stages)


def _shaped_noise(rng, n, shape):
    """Gaussian noise whose amplitude spectrum is `shape(f)` (smooth shapes only: a brick-wall band
    would put a step in the PSD that no 1/f model can fit, and tool 7 would flag every epoch)."""
    spec = np.fft.rfft(rng.standard_normal(n))
    f = np.fft.rfftfreq(n, 1 / SFREQ)
    spec *= shape(f)
    spec[0] = 0
    x = np.fft.irfft(spec, n)
    return x / x.std()


def _peak(center, width):
    """A rhythm = a smooth spectral bump (alpha, theta) on top of the 1/f background."""
    return lambda f: np.exp(-0.5 * ((f - center) / width) ** 2)


def _power_law(exponent):
    """Aperiodic activity: power ~ 1/f**exponent (amplitude ~ f**-exponent/2), floored at 0.3 Hz.
    Sleep EEG is close to a power law over 2-45 Hz, its exponent steepening with depth (N3 > N2 > W):
    a steeper component that grows with the stage is how deep sleep is modelled here, with no knee
    in the fitted band."""
    return lambda f: np.maximum(f, 0.3) ** (-exponent / 2.0)


def make_signal(stages, rng):
    """(n_channels, n_samples) in volts, following the hypnogram epoch by epoch."""
    n_ep = int(EPOCH_S * SFREQ)
    n = len(stages) * n_ep
    t = np.arange(n_ep) / SFREQ
    stage_per_sample = np.repeat(stages, n_ep)
    data = np.zeros((len(CHANNELS), n))
    # rhythm generators shared by the cortical channels (with channel-specific gains): the same
    # oscillation seen by several electrodes, as in a real montage
    alpha = _shaped_noise(rng, n, _peak(10.0, 0.8))
    theta = _shaped_noise(rng, n, _peak(5.5, 1.2))
    slow = _shaped_noise(rng, n, _power_law(3.0))
    fast = _shaped_noise(rng, n, _peak(25.0, 6.0))
    gains = {'Fp1': dict(alpha=0.4, slow=1.2, spindle=0.6),
             'C3': dict(alpha=0.7, slow=1.0, spindle=1.0),
             'O1': dict(alpha=1.5, slow=0.8, spindle=0.5)}
    spindles = np.zeros(n)
    for ei in np.flatnonzero(stages == 'N2'):
        for _ in range(int(rng.integers(1, 4))):                     # 1-3 spindles per N2 epoch
            dur = rng.uniform(0.7, 1.5)
            m = int(dur * SFREQ)
            start = ei * n_ep + int(rng.integers(0, n_ep - m))
            freq = rng.uniform(12.0, 14.0)
            spindles[start:start + m] += np.sin(2 * np.pi * freq * t[:m]) * np.hanning(m)
    # standard deviation in uV per stage (alpha, theta, slow waves, fast); the 1/f background (8 uV)
    # comes on top. N3 ~ 25 uV sd gives the usual 100-200 uV peak-to-peak slow waves.
    amp = {'W': (7, 2, 3, 3), 'N1': (2, 5, 6, 1), 'N2': (1, 3, 12, 1),
           'N3': (1, 3, 24, 0.5), 'R': (2, 5, 4, 1.5)}
    for st, (a_al, a_th, a_sl, a_fa) in amp.items():
        mask = stage_per_sample == st
        for ci, ch in enumerate(CHANNELS):
            if ch == 'A2':
                continue
            g = gains[ch]
            data[ci, mask] += (a_al * g['alpha'] * alpha[mask] + a_th * theta[mask]
                               + a_sl * g['slow'] * slow[mask] + a_fa * fast[mask])
    for ci, ch in enumerate(CHANNELS):
        background = _shaped_noise(rng, n, _power_law(2.0))
        if ch == 'A2':                                                # mastoid: quiet reference site
            data[ci] = 8.0 * background
        else:
            data[ci] += 10.0 * background + 15.0 * gains[ch]['spindle'] * spindles
    # amplifier noise floor (white, ~0.7 uV): without it the smooth 1/f signal changes so little from
    # one sample to the next that tool 6's flat-signal metric flags healthy channels
    data += 0.7 * rng.standard_normal(data.shape)
    return data * 1e-6


def make_events(stages, rng):
    """Event CSV rows (index, Name, Start, Duration) drawn inside the sleep epochs, sorted by onset."""
    sleep_epochs = np.flatnonzero(stages != 'W')
    rows = []
    for label, count in EVENT_LABELS:
        for ei in rng.choice(sleep_epochs, size=count, replace=False):
            start = ei * EPOCH_S + rng.uniform(0, EPOCH_S - 1)
            dur = rng.uniform(3, 15) if 'Arousal' in label else rng.uniform(10, 30)
            rows.append((label, round(float(start), 2), round(float(dur), 5)))
    rows.sort(key=lambda r: r[1])
    return rows


def write_recording(folder, file_id):
    """Write {file_id}.edf + _Hypnogram_remapped.txt + _event_xml.csv into `folder`."""
    seed, build, phys_overrides = DEFECTS[file_id]
    rng = np.random.default_rng(seed)
    stages = make_hypnogram(rng)
    info = mne.create_info(CHANNELS, SFREQ, ch_types='eeg')
    raw = build(mne.io.RawArray(make_signal(stages, rng), info, verbose=False))
    signals = []
    for ch in CHANNELS:
        phys = (phys_overrides or {}).get(ch, PHYS_RANGE_UV)
        data_uv = np.clip(raw.get_data(picks=[ch])[0] * 1e6, phys[0], phys[1])
        signals.append(edfio.EdfSignal(data_uv, sampling_frequency=SFREQ, label=ch,
                                       physical_dimension='uV', physical_range=phys))
    # anonymous by construction: edfio writes 'X X X X' patient / recording fields and no start date
    # beyond its 1985 default
    edfio.Edf(signals, data_record_duration=1).write(str(Path(folder) / f'{file_id}.edf'))
    (Path(folder) / f'{file_id}_Hypnogram_remapped.txt').write_text('\n'.join(stages) + '\n',
                                                                    encoding='utf-8')
    lines = [',Name,Start,Duration'] + [f'{i},{n},{s},{d}' for i, (n, s, d)
                                        in enumerate(make_events(stages, rng))]
    (Path(folder) / f'{file_id}_event_xml.csv').write_text('\n'.join(lines) + '\n', encoding='utf-8')
