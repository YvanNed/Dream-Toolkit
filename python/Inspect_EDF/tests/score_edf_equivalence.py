r"""Score_EDF -> tool 3bis: does the new tool reproduce what Score_EDF computed? (tools/plan_score_edf_fixes.md)

Not a pytest test: a before / after check run by hand, like outside_chain.py. Both tools are run in the SAME
environment (autoscore: the one Score_EDF was written for, YASA 0.7), so a difference comes from the code only.

    $py = "$env:LOCALAPPDATA\miniforge3\envs\autoscore\python.exe"
    & $py tests/score_edf_equivalence.py build            # simulated study: 3 nights of 2 h, full PSG montage
    & $py tests/score_edf_equivalence.py baseline-sim     # Score_EDF's notebook, UNMODIFIED, on the simulated nights
    & $py tests/score_edf_equivalence.py baseline-real    # Score_EDF's functions on night 100 (tools/test_data), timed
    & $py tests/score_edf_equivalence.py run-3bis-sim     # tool 3bis (step 1 + step 2) on the same nights
    & $py tests/score_edf_equivalence.py run-3bis-real
    & $py tests/score_edf_equivalence.py compare          # tables side by side, differences listed

Two datasets:
- simulated (3 nights): every channel Score_EDF looks for is present under a name its alias list knows, so its
  notebook runs untouched (YASA, GSSC, SleepFM, recalibrations A and B). Nights of exactly 2 h (a multiple of
  SleepFM's 5-min window): no "tail" pass, so SleepFM must match exactly.
- night 100 (real, tools/test_data, never versioned): F4 is absent, which Score_EDF's notebook cannot handle
  (KeyError 'F4'); its functions are called directly with the derivations present, like its notebook would.
  One night: no recalibration. Its SleepFM tail is where tool 3bis's shorter tail context shows (plan 3.2).

Channel semantics shared by the two tools in this check: Score_EDF's vEOG = the LEFT EOG (its aliases: LOC,
E1, EOG Left), hEOG = the RIGHT EOG. Its GSSC and YASA use hEOG alone: tool 3bis is set to "Right only".
"""
import json
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
SCORE_EDF = REPO_ROOT.parent / 'Score_EDF'
sys.path.insert(0, str(HERE))

ROOT = Path(os.path.realpath(os.path.join(tempfile.gettempdir(), 'dtk_scoreedf')))
SIM_IDS = ['simA', 'simB', 'simC']
SIM_SEEDS = {'simA': 201, 'simB': 202, 'simC': 203}
SIM_EPOCHS = 240                       # 2 h: a multiple of SleepFM's 5-min window (no tail pass)
EEG_FS = 256
REAL_EDF = REPO_ROOT / 'tools' / 'test_data' / '100.edf'
REAL_HYPNO = REPO_ROOT / 'tools' / 'test_data' / '100_Hypnogram_remapped.txt'
SLEEPFM_ROOT = Path(os.environ.get('SLEEPFM_ROOT', Path.home() / 'sleepfm-clinical'))

# Night 100 in Score_EDF's standard names (what a user picks in its channel-mapping dropdowns)
REAL_MAPPING = {'C3': 'C3', 'C4': 'C4', 'F3': 'F3', 'F4': None, 'O1': 'O1', 'O2': 'O2', 'M1': 'A1', 'M2': 'A2',
                'vEOG': 'EOG G', 'hEOG': 'EOG D', 'EMG': 'EMG 1', 'ECG': 'ECG', 'Thoracic': 'Thoracic',
                'Abdominal': 'Abdominal', 'Flux': 'Flux', 'SpO2': 'SpO2', 'Pulse': 'Pulse'}
REAL_GSSC_EEG = ('C3-M2', 'C4-M1', 'F3-M2')                     # no F4
REAL_SLEEPFM_EEG = ['F3-M2', 'C3-M2', 'C4-M1', 'O1-M2', 'O2-M1']   # EEG_DERIVATIONS order, without F4-M1


# =====================================================================================================
# Simulated nights
# =====================================================================================================
def _hypnogram(rng, n_epochs):
    """W at lights-off, then short cycles (N1 -> N2 -> N3 -> N2 -> R), a few W at the end."""
    stages = ['W'] * int(rng.integers(8, 16))
    cycle = 0
    while len(stages) < n_epochs - 6:
        stages += (['N1'] * int(rng.integers(3, 7)) + ['N2'] * int(rng.integers(15, 25))
                   + ['N3'] * max(0, 30 - 12 * cycle + int(rng.integers(-4, 5)))
                   + ['N2'] * int(rng.integers(8, 15)) + ['R'] * (10 + 8 * cycle + int(rng.integers(-2, 3))))
        if rng.random() < 0.5:
            stages += ['W'] * int(rng.integers(1, 4))
        cycle += 1
    return np.array(stages[:n_epochs - 6] + ['W'] * 6)


def _shaped(rng, n, fs, shape):
    spec = np.fft.rfft(rng.standard_normal(n))
    spec *= shape(np.fft.rfftfreq(n, 1 / fs))
    spec[0] = 0
    x = np.fft.irfft(spec, n)
    return x / x.std()


def _peak(c, w):
    return lambda f: np.exp(-0.5 * ((f - c) / w) ** 2)


def _power_law(e):
    return lambda f: np.maximum(f, 0.3) ** (-e / 2.0)


def make_sim_night(seed):
    """(signals {label: (data, fs, dimension, physical range)}, stages). Monopolar EEG on a common recording
    reference (C3, C4, F3, F4, O1, O2 + mastoids M1, M2), left / right EOG, chin EMG, ECG, respiration,
    SpO2 and pulse, each following the hypnogram well enough for the three models to score it."""
    rng = np.random.default_rng(seed)
    stages = _hypnogram(rng, SIM_EPOCHS)
    fs = EEG_FS
    n_ep = 30 * fs
    n = len(stages) * n_ep
    per_sample = np.repeat(stages, n_ep)
    t = np.arange(n) / fs
    alpha, theta = _shaped(rng, n, fs, _peak(10, 0.8)), _shaped(rng, n, fs, _peak(5.5, 1.2))
    slow, fast = _shaped(rng, n, fs, _power_law(3.0)), _shaped(rng, n, fs, _peak(25, 6))
    spindles = np.zeros(n)
    for ei in np.flatnonzero(stages == 'N2'):
        for _ in range(int(rng.integers(1, 4))):
            m = int(rng.uniform(0.7, 1.5) * fs)
            s = ei * n_ep + int(rng.integers(0, n_ep - m))
            spindles[s:s + m] += np.sin(2 * np.pi * rng.uniform(12, 14) * np.arange(m) / fs) * np.hanning(m)
    amp = {'W': (7, 2, 3, 3), 'N1': (2, 5, 6, 1), 'N2': (1, 3, 12, 1), 'N3': (1, 3, 24, 0.5), 'R': (2, 5, 4, 1.5)}
    gains = {'F3': (0.4, 1.3, 0.7), 'F4': (0.4, 1.3, 0.7), 'C3': (0.7, 1.0, 1.0), 'C4': (0.7, 1.0, 1.0),
             'O1': (1.5, 0.8, 0.5), 'O2': (1.5, 0.8, 0.5)}           # (alpha, slow, spindle)
    common_ref = 6.0 * _shaped(rng, n, fs, _power_law(2.0))          # picked up by every monopolar lead
    sig = {}
    for ch, (g_al, g_sl, g_sp) in gains.items():
        x = np.zeros(n)
        for st, (a_al, a_th, a_sl, a_fa) in amp.items():
            m = per_sample == st
            x[m] = a_al * g_al * alpha[m] + a_th * theta[m] + a_sl * g_sl * slow[m] + a_fa * fast[m]
        x += 9.0 * _shaped(rng, n, fs, _power_law(2.0)) + 15.0 * g_sp * spindles + common_ref
        sig[ch] = x
    for ch in ('M1', 'M2'):                                          # mastoids: reference + a little EEG
        sig[ch] = common_ref + 3.0 * _shaped(rng, n, fs, _power_law(2.0)) + 0.2 * slow * 10
    # Eye movements: slow rolling in W / N1, rapid saccades in R, opposite polarity left / right.
    eye = np.zeros(n)
    rolling = _shaped(rng, n, fs, _peak(0.3, 0.15))
    for st, a in (('W', 40.0), ('N1', 60.0)):
        eye[per_sample == st] += a * rolling[per_sample == st]
    for ei in np.flatnonzero(stages == 'R'):
        for _ in range(int(rng.integers(2, 6))):
            s = ei * n_ep + int(rng.integers(0, n_ep - fs))
            step = rng.choice([-1, 1]) * rng.uniform(80, 150)
            eye[s:s + fs // 2] += step * np.hanning(fs // 2)
    leak = 0.3 * sig['F3']
    sig['LOC'] = eye + leak + 2.0 * rng.standard_normal(n)
    sig['ROC'] = -eye + leak + 2.0 * rng.standard_normal(n)
    # Chin EMG: tone W > N1 > N2 > N3 >> R (atonia)
    emg_amp = {'W': 20.0, 'N1': 12.0, 'N2': 8.0, 'N3': 6.0, 'R': 1.5}
    emg_noise = _shaped(rng, n, fs, lambda f: ((f > 10) & (f < 100)).astype(float) + 0.01)
    sig['Chin'] = np.array([emg_amp[s] for s in per_sample]) * emg_noise + 0.5 * rng.standard_normal(n)
    # ECG: QRS-like spikes, ~60 bpm with variability
    ecg = np.zeros(n)
    beat = 0.0
    while beat < t[-1] - 1:
        i = int(beat * fs)
        ecg[i:i + 12] += 1000.0 * np.hanning(12)
        beat += 1.0 + 0.05 * rng.standard_normal()
    sig['ECG'] = ecg + 10.0 * rng.standard_normal(n)
    eeg_like = {k: (v, fs, 'uV', (-1000.0, 1000.0)) for k, v in sig.items() if k not in ('Chin', 'ECG')}
    eeg_like['Chin'] = (sig['Chin'], fs, 'uV', (-500.0, 500.0))
    eeg_like['ECG'] = (sig['ECG'], fs, 'uV', (-3000.0, 3000.0))
    # Respiration at 32 Hz: regular breathing in NREM, irregular in W / R
    rfs = 32
    rn = len(stages) * 30 * rfs
    rt = np.arange(rn) / rfs
    r_stage = np.repeat(stages, 30 * rfs)
    jitter = np.cumsum(np.where(np.isin(r_stage, ['W', 'R']), 0.02, 0.003) * rng.standard_normal(rn))
    breath = np.sin(2 * np.pi * 0.25 * rt + jitter)
    eeg_like['Thoracic'] = (300.0 * breath + 5.0 * rng.standard_normal(rn), rfs, 'uV', (-1000.0, 1000.0))
    eeg_like['Abdominal'] = (250.0 * np.sin(2 * np.pi * 0.25 * rt + jitter + 0.2) + 5.0 * rng.standard_normal(rn),
                             rfs, 'uV', (-1000.0, 1000.0))
    eeg_like['Airflow'] = (np.gradient(breath) * 2000.0 + 5.0 * rng.standard_normal(rn), rfs, 'uV', (-1000.0, 1000.0))
    sfs, sn = 1, len(stages) * 30
    eeg_like['SpO2'] = (96.0 + 0.5 * rng.standard_normal(sn), sfs, '%', (0.0, 100.0))
    eeg_like['Pulse'] = (60.0 + 3.0 * rng.standard_normal(sn), sfs, 'bpm', (0.0, 250.0))
    return eeg_like, stages


def write_sim_night(folder, file_id):
    import edfio
    signals, stages = make_sim_night(SIM_SEEDS[file_id])
    out = []
    for label, (data, fs, dim, rng_) in signals.items():
        out.append(edfio.EdfSignal(np.clip(data, *rng_), sampling_frequency=fs, label=label,
                                   physical_dimension=dim, physical_range=rng_))
    edfio.Edf(out, data_record_duration=1).write(str(Path(folder) / f'{file_id}.edf'))
    (Path(folder) / f'{file_id}_Hypnogram_remapped.txt').write_text('\n'.join(stages) + '\n', encoding='utf-8')


def tool2_entry(eeg, mastoids, eog_left, eog_right, emg, ecg):
    """A participant entry of tool 2's remap_reref_persubject.json (identity remap, linked mastoids)."""
    return {'config': 'config. 1', 'remap': {ch: ch for ch in eeg + mastoids}, 'ref_channels': list(mastoids),
            'context_channels': {'eog_left': eog_left, 'eog_right': eog_right, 'emg': emg, 'ecg': ecg}}


def build():
    data = ROOT / 'sim' / 'study'
    if data.exists():
        shutil.rmtree(data)
    (data / 'config_param').mkdir(parents=True)
    cfg = {}
    for fid in SIM_IDS:
        write_sim_night(data, fid)
        cfg[fid] = tool2_entry(['C3', 'C4', 'F3', 'F4', 'O1', 'O2'], ['M1', 'M2'], 'LOC', 'ROC', 'Chin', 'ECG')
    (data / 'config_param' / 'remap_reref_persubject.json').write_text(json.dumps(cfg, indent=2), encoding='utf-8')
    # night 100: a study folder holding a copy of the EDF + its tool-3 hypnogram
    real = ROOT / 'real' / 'study'
    if real.exists():
        shutil.rmtree(real)
    (real / 'config_param').mkdir(parents=True)
    # Night 100's only 512-Hz channel is DHR (a derived heart rate; every other channel is at 256 Hz or less).
    # Score_EDF reads all 58 channels, so MNE brings them all to 512 Hz: its EEG is upsampled, tool 3bis's (which
    # reads only the channels it uses) is not. Without DHR both tools get the same 256-Hz signals.
    import edfio
    edf = edfio.read_edf(str(REAL_EDF))
    edf.drop_signals(['DHR'])
    edf.write(str(real / REAL_EDF.name))
    shutil.copy2(REAL_HYPNO, real / REAL_HYPNO.name)
    entry = tool2_entry(['C3', 'C4', 'F3', 'O1', 'O2'], ['A1', 'A2'], 'EOG G', 'EOG D', 'EMG 1', 'ECG')
    (real / 'config_param' / 'remap_reref_persubject.json').write_text(json.dumps({'100': entry}, indent=2),
                                                                       encoding='utf-8')
    print(f'built {data} ({len(SIM_IDS)} nights) and {real} (night 100)')


# =====================================================================================================
# Score_EDF (reference)
# =====================================================================================================
def _score_edf_copy(study, dest):
    """Score_EDF scans a folder for EDF + hypnograms and writes autoscore_outputs/ beside it: give it a clean
    copy holding only the EDF files and the expert hypnograms."""
    if dest.exists():
        shutil.rmtree(dest)
    (dest / 'data').mkdir(parents=True)
    for f in study.iterdir():
        if f.suffix.lower() == '.edf' or f.name.endswith('_Hypnogram_remapped.txt'):
            shutil.copy2(f, dest / 'data' / f.name)
    return dest / 'data'


def baseline_sim():
    """Score_EDF's own notebook, untouched, driven like a user: pick the folders, Scan, Run (every stager)."""
    import nbdriver
    nbdriver.CELL_TIMEOUT_S = 4 * 3600          # Score_EDF scores the 3 nights in one click
    os.environ['SLEEPFM_ROOT'] = str(SLEEPFM_ROOT)
    data = _score_edf_copy(ROOT / 'sim' / 'study', ROOT / 'sim' / 'score_edf')
    t0 = time.time()
    with nbdriver.NotebookSession(SCORE_EDF / 'score_edf_voila.ipynb', cwd=SCORE_EDF) as nb:
        nb.pick('chooser_edf', data)
        nb.pick('chooser_scoring', data)
        nb.click('scan_button')
        ticked = nb.value("[k for k, cb in STATE['checkboxes'].items() if cb.value]")
        nb.click('run_button')
        text = nb.widget_text()
    (ROOT / 'sim' / 'score_edf' / 'run_log.txt').write_text(text, encoding='utf-8')
    print(f'Score_EDF on the simulated nights: stagers {ticked}, {time.time() - t0:.0f} s')
    print(text[-2500:])


class _Timer:
    """Wall time and process memory of each step (psutil: resident set now, peak so far)."""
    def __init__(self):
        import psutil
        self.p = psutil.Process()
        self.rows = []

    def step(self, name, fn):
        t0 = time.perf_counter()
        out = fn()
        mem = self.p.memory_info()
        self.rows.append({'step': name, 'seconds': round(time.perf_counter() - t0, 1),
                          'rss_gb': round(mem.rss / 1e9, 2),
                          'peak_gb': round(getattr(mem, 'peak_wset', mem.rss) / 1e9, 2)})
        print(self.rows[-1], flush=True)
        return out


def baseline_real():
    """Score_EDF's functions on night 100, in the order its notebook calls them, with the derivations present
    (no F4). Timed step by step: the reference for the speed audit (plan phase 0)."""
    os.environ['SLEEPFM_ROOT'] = str(SLEEPFM_ROOT)
    sys.path.insert(0, str(SCORE_EDF))
    from scoredf.io_score import detect_and_read_scoring, window_from_scoring
    from scoredf.channels import channel_mapping
    from scoredf.plots import (compute_agreement, plot_agreement_summary, plot_confusion_matrices,
                               plot_hypnodensities, plot_hypnogram_grid)
    from sleephd.alignment import build_epoch_table
    from sleephd.signals import infer_sleepfm_montage
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    data = _score_edf_copy(ROOT / 'real' / 'study', ROOT / 'real' / 'score_edf')
    out_dir = data.parent / 'autoscore_outputs'
    cache = out_dir / 'cache'
    cache.mkdir(parents=True, exist_ok=True)
    edf = data / '100.edf'
    sid = '100'
    timer = _Timer()
    hyp, ep_s = detect_and_read_scoring(data / REAL_HYPNO.name)
    window = window_from_scoring(hyp, edf, ep_s)
    five = {}
    with channel_mapping(REAL_MAPPING):
        # Imported here, as Score_EDF's notebook does: stagers.py binds signals._load_edf_channels at its first
        # import, and only an import made inside channel_mapping gets the renaming version of it.
        from sleephd.stagers import infer_gssc, infer_yasa
        five['YASA'] = timer.step('YASA (read + score)', lambda: infer_yasa(
            edf, window, cache_path=cache / f'{sid}_yasa.csv', force=True))
        five['GSSC'] = timer.step('GSSC (read + 7 combinations)', lambda: infer_gssc(
            edf, window, eeg=REAL_GSSC_EEG, cache_path=cache / f'{sid}_gssc.csv', force=True))
        five['SleepFM'] = timer.step('SleepFM (read + montage + primary + tail)', lambda: infer_sleepfm_montage(
            edf, cache, eeg_derivations=REAL_SLEEPFM_EEG, cache_tag='aasm_all_eeg', include_eog=True,
            device='cpu', cover_tail=True, force=True))
    frames = []
    for name in ['YASA', 'GSSC', 'SleepFM']:
        table = build_epoch_table(window, None, hyp, hypnodensity_5s=five[name], epoch_seconds=ep_s)
        table['minutes_from_start'] = (table['start_seconds'] - window.lights_off_seconds) / 60
        frames.append(table.assign(subject=sid, stager=name))
    epochs_df = pd.concat(frames, ignore_index=True)
    epochs_df.to_csv(out_dir / 'epochs_all_stagers.csv', index=False)
    scored = epochs_df.loc[epochs_df['manual_stage'] != 'Unknown']
    agreement = compute_agreement(scored, ['YASA', 'GSSC', 'SleepFM'], [sid])
    agreement.to_csv(out_dir / 'agreement_summary.csv')
    for name, fig in [('agreement_by_stager', plot_agreement_summary(agreement, ['YASA', 'GSSC', 'SleepFM'], [sid])),
                      ('confusion_by_stager', plot_confusion_matrices(scored, ['YASA', 'GSSC', 'SleepFM'])),
                      ('hypnograms_by_stager', plot_hypnogram_grid(epochs_df, ['YASA', 'GSSC', 'SleepFM'], [sid])),
                      (f'hypnodensities_{sid}', plot_hypnodensities(epochs_df, ['YASA', 'GSSC', 'SleepFM'], sid))]:
        fig.savefig(out_dir / f'{name}.png', dpi=130, bbox_inches='tight')
        plt.close(fig)
    pd.DataFrame(timer.rows).to_csv(ROOT / 'real' / 'score_edf' / 'timings.tsv', sep='\t', index=False)
    print(agreement.round(3).to_string())


# =====================================================================================================
# Tool 3bis (new)
# =====================================================================================================
T3BIS = 'tools/3bis_autoscoring_voila.ipynb'


def _clean_outputs(study):
    """Remove what an earlier 3bis run wrote, so this run starts from the inputs only."""
    for name in ('derivatives', 'autoscore_outputs', 'autoscore_outputs_strict'):
        shutil.rmtree(study / name, ignore_errors=True)
    for d in study.glob('reports_autoscoring_*'):
        shutil.rmtree(d, ignore_errors=True)
    # Only the models' own hypnograms (on Windows a glob ignores case: '*_hypnogram_*' would take the
    # expert '*_Hypnogram_remapped.txt' too)
    for algo in ('gssc', 'yasa', 'sleepfm'):
        for f in study.glob(f'*_hypnogram_{algo}.txt'):
            f.unlink()


def run_3bis(dataset, gssc_eeg):
    """Tool 3bis driven like a user, with the settings that reproduce Score_EDF's inputs: the three models,
    contralateral mastoid reference, GSSC on `gssc_eeg` + the right EOG alone (Score_EDF's hEOG), YASA on C4 +
    right EOG + chin EMG. Section 4 is run twice: first with Score_EDF's recalibrations (N3 fixed: must match it
    exactly, saved as autoscore_outputs_strict/), then with the default (N3 fitted when present)."""
    import nbdriver
    nbdriver.CELL_TIMEOUT_S = 4 * 3600
    os.environ['SLEEPFM_ROOT'] = str(SLEEPFM_ROOT)
    study = ROOT / dataset / 'study'
    _clean_outputs(study)
    t0 = time.time()
    with nbdriver.NotebookSession(T3BIS) as nb:
        nb.pick('fc_data', study)
        nb.set('cb_yasa', True)
        nb.set('cb_sleepfm', True)
        nb.click('btn_scan')
        nb.set('dd_ref', 'Contralateral mastoid (as the models were trained)')
        nb.set('tb_mode', 'Several electrodes (consensus)')
        nb.run(f'for _ch, _cb in eeg_checks.items(): _cb.value = _ch in {list(gssc_eeg)!r}')
        nb.set('dd_eog', 'Right only')
        nb.set('dd_yasa_eeg', 'C4')
        nb.set('dd_yasa_eog', 'Right')
        nb.click('btn_run')
        step1 = time.time() - t0
        print(nb.value('_dtk_widget_text()')[-3000:], flush=True)
        nb.run('RECAL_FREE_N3 = False')
        nb.click('btn_compare_find')
        ticked = nb.value("[k for k, cb in cmp_checks.items() if cb.value]")
        print('find:', nb.value('cmp_info.value'), ticked, flush=True)
        nb.click('btn_compare')
        if not (study / 'autoscore_outputs').is_dir():
            raise RuntimeError('Section 4 wrote nothing:\n' + nb.widget_text()[-3000:])
        shutil.copytree(study / 'autoscore_outputs', study / 'autoscore_outputs_strict')
        nb.run("RECAL_FREE_N3 = 'auto'")
        nb.click('btn_compare')
        text = nb.widget_text()
    (ROOT / dataset / 'run_3bis_log.txt').write_text(text, encoding='utf-8')
    print(f'tool 3bis on {dataset}: step 1 in {step1:.0f} s, total {time.time() - t0:.0f} s; compared {ticked}')
    print(text[-2500:])


def run_3bis_sim():
    run_3bis('sim', ['C3', 'C4', 'F3', 'F4'])


def run_3bis_real():
    run_3bis('real', ['C3', 'C4', 'F3'])


# =====================================================================================================
# Comparison
# =====================================================================================================
KEY = ['subject', 'stager', 'epoch']
PROBS = ['Wake', 'N1', 'N2', 'N3', 'REM']


def _epochs(path):
    d = pd.read_csv(path, dtype={'subject': str})
    return d.sort_values(KEY).reset_index(drop=True)


def compare_tables(ref_dir, new_dir, label, lines, expect_close=()):
    """Score_EDF's epochs_all_stagers.csv / agreement_summary.csv against tool 3bis's, model by model."""
    ref, new = _epochs(ref_dir / 'epochs_all_stagers.csv'), _epochs(new_dir / 'epochs_all_stagers.csv')
    lines.append(f'\n## {label}')
    same_cols = list(ref.columns) == list(new.columns)
    lines.append(f'columns identical: {same_cols}' + ('' if same_cols else f'\n  Score_EDF {list(ref.columns)}'
                                                                           f'\n  3bis      {list(new.columns)}'))
    for stager in sorted(set(ref['stager']) | set(new['stager'])):
        r, n = ref[ref['stager'] == stager], new[new['stager'] == stager]
        if len(r) != len(n) or not (r[KEY].values == n[KEY].values).all():
            lines.append(f'- {stager}: DIFFERENT epochs ({len(r)} vs {len(n)} rows)')
            continue
        dp = np.nanmax(np.abs(r[PROBS].to_numpy(float) - n[PROBS].to_numpy(float))) if len(r) else 0.0
        stage_diff = int((r['predicted_stage'].values != n['predicted_stage'].values).sum())
        manual_diff = int((r['manual_stage'].values != n['manual_stage'].values).sum())
        verdict = 'IDENTICAL' if dp == 0 and stage_diff == 0 else (
            'same stages, probabilities within 5e-5 (step 1 stores them with 4 decimals)'
            if stage_diff == 0 and dp <= 5e-5 else 'DIFFERENT')
        if verdict == 'DIFFERENT' and stager in expect_close:
            verdict = 'different (intended, see below)'
        lines.append(f'- {stager}: {len(r)} epochs, predicted stage differs on {stage_diff}, expert stage differs '
                     f'on {manual_diff}, max |probability difference| = {dp:.2e} -> {verdict}')
    ra = pd.read_csv(ref_dir / 'agreement_summary.csv', index_col=0)
    na = pd.read_csv(new_dir / 'agreement_summary.csv', index_col=0)
    common = [c for c in ra.columns if c in na.columns]
    diff = (ra[common] - na.reindex(ra.index)[common]).abs().max().max()
    lines.append(f'agreement_summary.csv: same rows {list(ra.index) == list(na.index)}, same columns '
                 f'{list(ra.columns) == list(na.columns)}, max |difference| = {diff:.2e}')
    lines.append('```\n' + pd.concat({'Score_EDF': ra[['accuracy', 'kappa_pooled']],
                                     '3bis': na[['accuracy', 'kappa_pooled']]}, axis=1).round(4).to_string() + '\n```')


def compare_tail(lines):
    """Night 100: SleepFM's 5-s probabilities. Score_EDF's tail used the whole night as context, tool 3bis's 30
    minutes (plan 3.2): the primary steps must be identical, the few tail steps close."""
    cache = ROOT / 'real' / 'score_edf' / 'autoscore_outputs' / 'cache'
    prim = pd.read_csv(cache / '100_aasm_all_eeg_sleepfm_hypnodensity_5s.csv')
    new = pd.read_csv(next((ROOT / 'real' / 'study').rglob('100_hypnodensity_5s.tsv')), sep='\t')
    cols_new = ['p_W', 'p_N1', 'p_N2', 'p_N3', 'p_R']
    n = len(prim)
    dp = np.abs(prim[PROBS].to_numpy() - new.loc[:n - 1, cols_new].to_numpy()).max()
    tail = new[new['source'] == 'tail']
    lines.append(f'\n## Night 100, SleepFM 5-s probabilities')
    lines.append(f'- primary pass: {n} steps, max |difference| = {dp:.2e}')
    lines.append(f'- tail: {len(tail)} steps added by 3bis (from {tail["start_s"].min():.0f} s), 30 min of context')


def compare():
    lines = ['# Score_EDF vs tool 3bis', f'({time.strftime("%Y-%m-%d %H:%M")}, autoscore environment)']
    compare_tables(ROOT / 'sim' / 'score_edf' / 'autoscore_outputs', ROOT / 'sim' / 'study' / 'autoscore_outputs_strict',
                   'Simulated nights, Score_EDF settings (N3 fixed in the recalibrations): must be identical', lines)
    compare_tables(ROOT / 'sim' / 'score_edf' / 'autoscore_outputs', ROOT / 'sim' / 'study' / 'autoscore_outputs',
                   'Simulated nights, tool 3bis defaults (N3 fitted in the recalibrations, plan 2.7)', lines,
                   expect_close=('SleepFM recal A', 'SleepFM recal B'))
    compare_tables(ROOT / 'real' / 'score_edf' / 'autoscore_outputs', ROOT / 'real' / 'study' / 'autoscore_outputs',
                   'Night 100 (no F4; 30-min SleepFM tail context, plan 3.2)', lines, expect_close=('SleepFM',))
    compare_tail(lines)
    text = '\n'.join(lines)
    (ROOT / 'compare_report.md').write_text(text, encoding='utf-8')
    print(text)


COMMANDS = {'build': build, 'baseline-sim': baseline_sim, 'baseline-real': baseline_real,
            'run-3bis-sim': run_3bis_sim, 'run-3bis-real': run_3bis_real, 'compare': compare}

if __name__ == '__main__':
    for cmd in sys.argv[1:] or ['build']:
        t0 = time.time()
        COMMANDS[cmd]()
        print(f'[{cmd}] done in {time.time() - t0:.0f} s')
