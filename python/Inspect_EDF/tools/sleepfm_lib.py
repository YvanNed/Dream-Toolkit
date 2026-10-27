"""
SleepFM-Clinical for tool 3bis (automatic sleep scoring): find or download the official model, prepare its input
montage, run it. Ported from python/Score_EDF (sleephd/sleepfm.py, sleephd/signals.py), see
tools/plan_score_edf_fixes.md (part F).

SleepFM (Thapa et al. 2026, Nature Medicine) is not a Python package: it is a folder of code + pretrained weights
published on GitHub (licence CC BY-NC 4.0, non-commercial use). Only the 9 files the scoring needs are used, at
a fixed version (SLEEPFM_COMMIT). They are looked for in the folder named by the SLEEPFM_ROOT environment
variable when it is set, else in a cache folder of the user, where download_sleepfm() fetches them from the
official repository the first time (each file checked against its SHA-256 fingerprint).

What differs from Score_EDF (intended, plan part F): the BiLSTM features used by recalibration B are captured
during the scoring pass itself (Score_EDF ran SleepFM a second time for them), and the "tail" pass that covers the
last minutes of the night uses TAIL_CONTEXT_CHUNKS of context instead of the whole night (Score_EDF re-ran the
encoder on almost the whole recording for those few minutes). tail_context_chunks=None gives Score_EDF's tail.
"""
import hashlib
import json
import os
import sys
import tempfile
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.signal import butter, sosfiltfilt

SLEEPFM_COMMIT = '2bcbae04c3592f61352addb7ac3d4193f0a3ca25'
SLEEPFM_URL = f'https://raw.githubusercontent.com/zou-group/sleepfm-clinical/{SLEEPFM_COMMIT}/'
# (path in the repository, SHA-256, size in bytes): the code SleepFM's scoring imports, its channel
# dictionary and its two pretrained models (encoder + sleep-staging head)
SLEEPFM_FILES = [
    ('sleepfm/models/__init__.py', 'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855', 0),
    ('sleepfm/models/dataset.py', '31ecb7c839f6c3de4180a2bb3094d6a3a75a86f24ad63fc62743b6f3e7db858d', 35151),
    ('sleepfm/models/models.py', 'b4f66873f863c0205bba14d02ba27d54de623b038c674144f6ff04b742f81e9c', 15906),
    ('sleepfm/utils.py', '8103b1bb76d68617d3af1df6d2fe9cb9b12ffc2ed2a2532d212f646e8b7ec5a6', 2956),
    ('sleepfm/configs/channel_groups.json', '03faba1ed1fd296ba345f98484f1677d856497ae4abc9812035296ed5de8ba43', 16342),
    ('sleepfm/checkpoints/model_base/best.pt', 'ffc9fc10233ebc4d0aae71abce87070db51b8bac6e4b16a5c1b4401a5f73f799', 19444020),
    ('sleepfm/checkpoints/model_base/config.json', '8ca2c62066ccdb19de18b44794c6050c1fa69c6dc2c62ea2767095c2fc3c2922', 884),
    ('sleepfm/checkpoints/model_sleep_staging/best.pth', 'a742faf4dd9a1f5ebcccc464894081480217e683974750eb3f10e6c4d1b0fde9', 9351915),
    ('sleepfm/checkpoints/model_sleep_staging/config.json', 'ec814f93f479eb26bdb6fdfddcd6be56577a779f9c7b5b8375f4749d7b20d0ed', 873),
]
SLEEPFM_TOTAL_BYTES = sum(size for _, _, size in SLEEPFM_FILES)

TARGET_HZ = 128                      # SleepFM works at 128 Hz
STEP_S = 5                           # one prediction every 5 s
CHUNK_SAMPLES = 5 * 60 * TARGET_HZ   # its encoder reads 5-min chunks (an incomplete last chunk is dropped)
# Filters applied before SleepFM's own preprocessing (Score_EDF's choice, AASM ranges): without them 50-Hz mains
# and other > 35 Hz content pass SleepFM's 64-Hz anti-alias filter and the model gives almost no REM.
BANDPASS = {'BAS': (0.3, 35.0), 'EMG': (10.0, 100.0)}
TAIL_CONTEXT_CHUNKS = 6              # context of the tail pass: 6 chunks = 30 min (plan 3.2, Q5)
MODALITIES = ['BAS', 'RESP', 'EKG', 'EMG']
# Respiratory channels, found by name (Score_EDF's alias lists) -> name given to SleepFM
RESP_ALIASES = {
    'Thoracic': ['Thoracic', 'Thor', 'Chest', 'Thor Effort', 'Thorax', 'Effort THO', 'RIP Thoracic', 'THO',
                 'Chest Effort'],
    'Abdominal': ['Abdominal', 'Abdo', 'Abd', 'Abdomen', 'Abd Effort', 'Effort ABD', 'RIP Abdominal', 'ABD',
                  'Abdom Effort'],
    'Airflow': ['Flux', 'Airflow', 'Flow', 'Nasal', 'NasalFlow', 'Nasal Pressure', 'Nasal Flow', 'Thermistor',
                'AIRFLOW', 'Airflow-0'],
    'SpO2': ['SpO2', 'SaO2', 'Oximetry', 'SAO2', 'SPO2', 'Ox Status'],
    'Pulse': ['Pulse', 'Pleth', 'PPG', 'Pulse Rate', 'HR', 'PULSE'],
}

# torch, h5py and the packages SleepFM's own code imports (einops, loguru, yaml): without them SleepFM cannot run,
# but the module still imports (the tool greys SleepFM out and says how to install them).
SLEEPFM_DEPS_OK = False
SLEEPFM_DEPS_MSG = ''
try:
    import torch
    import h5py
    import einops  # noqa: F401  (imported by SleepFM's code)
    import loguru  # noqa: F401
    import yaml  # noqa: F401
    SLEEPFM_DEPS_OK = True
except Exception as _e:
    SLEEPFM_DEPS_MSG = (f'SleepFM needs the packages of the autoscoring environment ({_e}). Run this tool in the '
                        '<code>autoscore</code> environment, or add them to the current one: '
                        '<code>conda install -c conda-forge h5py einops loguru pyyaml</code>.')


# ---------------------------------------------------------------------------------------------------------
# Locate / download the model
# ---------------------------------------------------------------------------------------------------------
def cache_root():
    """Where download_sleepfm() puts the model: a cache folder of the user (short path, outside OneDrive)."""
    base = os.environ.get('LOCALAPPDATA') or str(Path.home() / '.cache')
    return Path(base) / 'DreamToolkit' / f'sleepfm-clinical-{SLEEPFM_COMMIT[:8]}'


def _complete(root):
    return all((Path(root) / rel).is_file() for rel, _, _ in SLEEPFM_FILES)


def find_sleepfm():
    """(folder, how it was found) when the model is on disk, else (None, why not). SLEEPFM_ROOT comes first."""
    env = os.environ.get('SLEEPFM_ROOT')
    if env:
        if _complete(env):
            return Path(env), f'SLEEPFM_ROOT ({env})'
        if _complete(cache_root()):
            return cache_root(), f'download cache ({cache_root()}); SLEEPFM_ROOT={env} is incomplete'
        return None, f'SLEEPFM_ROOT={env} does not hold the SleepFM files'
    if _complete(cache_root()):
        return cache_root(), f'download cache ({cache_root()})'
    return None, 'not downloaded yet'


def _sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def download_sleepfm(progress=None):
    """Download the SleepFM files from the official repository into cache_root(), checking each against its
    SHA-256 fingerprint (a file already there and correct is kept). progress(done_bytes, total_bytes, file) is
    called while downloading. Returns the folder. Raises with a readable message on a network or checksum error."""
    root = cache_root()
    done = 0
    for rel, sha, size in SLEEPFM_FILES:
        dest = root / rel
        if dest.is_file() and dest.stat().st_size == size and _sha256(dest) == sha:
            done += size
            if progress:
                progress(done, SLEEPFM_TOTAL_BYTES, rel)
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        part = dest.with_name(dest.name + '.part')      # renamed only once complete and verified
        h = hashlib.sha256()
        try:
            with urllib.request.urlopen(SLEEPFM_URL + rel, timeout=60) as resp, open(part, 'wb') as f:
                while True:
                    block = resp.read(1 << 20)
                    if not block:
                        break
                    f.write(block)
                    h.update(block)
                    done += len(block)
                    if progress:
                        progress(done, SLEEPFM_TOTAL_BYTES, rel)
        except Exception as e:
            part.unlink(missing_ok=True)
            raise RuntimeError(f'download of {rel} failed ({e}): check the internet connection') from None
        if h.hexdigest() != sha:
            part.unlink(missing_ok=True)
            raise RuntimeError(f'{rel}: downloaded file does not match its fingerprint (corrupted or changed)')
        os.replace(part, dest)
    return root


# ---------------------------------------------------------------------------------------------------------
# Input montage (Score_EDF's prepare_sleepfm_montage, minus the EDF reading)
# ---------------------------------------------------------------------------------------------------------
def bandpass(signal, sfreq, low, high):
    """Zero-phase 4th-order Butterworth band-pass; `high` is capped below Nyquist."""
    high = min(high, 0.99 * sfreq / 2.0)
    sos = butter(4, [low, high], btype='band', fs=sfreq, output='sos')
    return sosfiltfilt(sos, signal)


def anti_alias_resample_standardize(signal, sfreq):
    """SleepFM's own preprocessing: low-pass at 64 Hz when faster than 128 Hz, linear interpolation to 128 Hz,
    z-score over the whole night."""
    if sfreq > TARGET_HZ:
        sos = butter(4, (TARGET_HZ / 2.0) / (sfreq / 2.0), btype='low', output='sos')
        signal = sosfiltfilt(sos, signal)
    duration = len(signal) / sfreq
    old_time = np.linspace(0, duration, len(signal), endpoint=False)
    new_time = np.linspace(0, duration, int(duration * TARGET_HZ), endpoint=False)
    resampled = np.interp(new_time, old_time, signal)
    std = resampled.std()
    return (resampled - resampled.mean()) / std if std > 0 else resampled - resampled.mean()


def channel_groups(root):
    with open(Path(root) / 'sleepfm' / 'configs' / 'channel_groups.json', encoding='utf-8') as f:
        return json.load(f)


def modality_of(name, groups):
    return next((m for m in MODALITIES if name in groups.get(m, [])), None)


def write_montage(path, bas, emg, others, sfreq, filters=BANDPASS):
    """The SleepFM input file (HDF5, one 128-Hz z-scored dataset per channel). bas = {name: signal} (EEG
    derivations + EOG, band-passed together), emg = chin EMG signal or None, others = {name: signal} (ECG,
    respiration), all in volts on the same time grid at `sfreq`. Same steps and order as Score_EDF."""
    selected = dict(bas)
    if filters:
        selected = {name: bandpass(values, sfreq, *filters['BAS']) for name, values in selected.items()}
    if emg is not None:
        selected['EMG'] = bandpass(emg, sfreq, *filters['EMG']) if filters else emg
    selected.update(others)
    processed = {name: anti_alias_resample_standardize(values, sfreq) for name, values in selected.items()}
    if len({len(v) for v in processed.values()}) != 1:
        raise ValueError('SleepFM montage channels have inconsistent lengths')
    if not all(np.isfinite(v).all() for v in processed.values()):
        raise ValueError('SleepFM montage contains non-finite values (flat or corrupted channel?)')
    with h5py.File(path, 'w') as handle:
        for name, values in processed.items():
            handle.create_dataset(name, data=values.astype('float32'), compression='gzip')
        handle.attrs['sampling_frequency_hz'] = TARGET_HZ
        handle.attrs['bandpass_hz'] = json.dumps(filters or {})
    return path


# ---------------------------------------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------------------------------------
_MODELS = {}


def _load_json(path):
    with open(path, encoding='utf-8') as f:
        return json.load(f)


def _state_dict_for_model(checkpoint, model):
    """The checkpoints were saved from a DataParallel model or not: add / strip the 'module.' prefix to match."""
    state = checkpoint.get('state_dict', checkpoint) if isinstance(checkpoint, dict) else checkpoint
    expects_module = next(iter(model.state_dict().keys())).startswith('module.')
    has_module = next(iter(state)).startswith('module.')
    if has_module and not expects_module:
        state = {key.removeprefix('module.'): value for key, value in state.items()}
    elif expects_module and not has_module:
        state = {f'module.{key}': value for key, value in state.items()}
    return state


def load_models(root, device='cpu'):
    """SleepFM's encoder and sleep-staging head with their pretrained weights (loaded once per session)."""
    key = (str(root), device)
    if key in _MODELS:
        return _MODELS[key]
    package = str(Path(root) / 'sleepfm')
    if package not in sys.path:
        sys.path.insert(0, package)
    from models.dataset import SetTransformerDataset, collate_fn
    from models.models import SetTransformer, SleepEventLSTMClassifier
    ckpt = Path(root) / 'sleepfm' / 'checkpoints'
    base_config = _load_json(ckpt / 'model_base' / 'config.json')
    base_config['num_workers'] = 1
    base = SetTransformer(base_config['in_channels'], base_config['patch_size'], base_config['embed_dim'],
                          base_config['num_heads'], base_config['num_layers'],
                          pooling_head=base_config['pooling_head'], dropout=0.0).to(device)
    # Official weights from the pinned repository (checked against their fingerprint): weights_only=False is
    # what these checkpoints need since torch 2.6 made weights_only=True the default (plan 4.1).
    base.load_state_dict(_state_dict_for_model(
        torch.load(ckpt / 'model_base' / 'best.pt', map_location=device, weights_only=False), base))
    base.eval()
    head = SleepEventLSTMClassifier(**_load_json(ckpt / 'model_sleep_staging' / 'config.json')['model_params']).to(device)
    head.load_state_dict(_state_dict_for_model(
        torch.load(ckpt / 'model_sleep_staging' / 'best.pth', map_location=device, weights_only=False), head))
    head.eval()
    _MODELS[key] = {'base': base, 'head': head, 'base_config': base_config, 'groups': channel_groups(root),
                    'dataset_class': SetTransformerDataset, 'collate_fn': collate_fn, 'device': device}
    return _MODELS[key]


def run_pass(hdf5_path, models):
    """One SleepFM pass over a prepared file: (probabilities (steps, 5) in the order W, N1, N2, N3, R, the
    staging head's BiLSTM output (steps, 128) = the input of its final layer, kept for recalibration B)."""
    from torch.utils.data import DataLoader
    device = models['device']
    cfg = models['base_config']
    dataset = models['dataset_class'](cfg, models['groups'], hdf5_paths=[str(hdf5_path)], split='test')
    if len(dataset) == 0:
        raise ValueError("recording shorter than SleepFM's 5-minute encoder window")
    loader = DataLoader(dataset, batch_size=1, shuffle=False, num_workers=0, collate_fn=models['collate_fn'])
    chunks = defaultdict(list)
    with torch.inference_mode():
        for batch_data, mask_list, _, _, chunk_starts in loader:
            for modality, values, mask in zip(cfg['modality_types'], batch_data, mask_list):
                emb = models['base'](values.to(device, dtype=torch.float32),
                                     mask.to(device, dtype=torch.bool))[1].cpu().numpy()
                for start, e in zip(chunk_starts, emb):
                    chunks[modality].append((int(start), e))
    seqs = [np.concatenate([a for _, a in sorted(chunks[m], key=lambda it: it[0])], axis=0)
            for m in cfg['modality_types']]
    steps = min(map(len, seqs))
    x = torch.from_numpy(np.stack([s[:steps] for s in seqs], axis=0)).unsqueeze(0).to(device, dtype=torch.float32)
    mask = torch.zeros((1, len(seqs), steps), dtype=torch.bool, device=device)
    captured = {}
    hook = models['head'].lstm.register_forward_hook(
        lambda module, inputs, output: captured.__setitem__('lstm', output[0]))
    try:
        with torch.inference_mode():
            logits, _ = models['head'](x, mask)
    finally:
        hook.remove()
    return torch.softmax(logits, dim=-1)[0].cpu().numpy(), captured['lstm'][0].cpu().numpy()


def head_parameters(models):
    """Weight and bias of the head's final linear layer (recalibration B starts from them)."""
    fc = models['head'].fc
    return fc.weight.detach().cpu().numpy(), fc.bias.detach().cpu().numpy()


def score_montage(hdf5_path, models, tail_context_chunks=TAIL_CONTEXT_CHUNKS):
    """SleepFM on a prepared file, the last incomplete 5-min chunk covered by a second pass whose chunks end at
    the end of the recording (only its predictions after the first pass are kept, the first one clipped to
    where the first pass stops). Returns (starts_s, stops_s, probs, features, source per step)."""
    probs, feats = run_pass(hdf5_path, models)
    starts = np.arange(len(probs)) * float(STEP_S)
    stops = starts + STEP_S
    source = np.array(['primary'] * len(probs), dtype=object)
    with h5py.File(hdf5_path, 'r') as handle:
        n_samples = min(len(handle[name]) for name in handle)
    offset = n_samples % CHUNK_SAMPLES
    if offset < STEP_S * TARGET_HZ:
        return starts, stops, probs, feats, source
    # Score_EDF started the tail pass at `offset` (whole night again); a few chunks of context are enough to
    # predict the last minutes. Same alignment: the chunks end exactly at the end of the recording.
    tail_start = offset if tail_context_chunks is None else max(offset, n_samples - tail_context_chunks * CHUNK_SAMPLES)
    tail_path = Path(hdf5_path).with_name(Path(hdf5_path).stem + '_tail.hdf5')
    try:
        with h5py.File(hdf5_path, 'r') as src, h5py.File(tail_path, 'w') as dst:
            for name in src:
                dst.create_dataset(name, data=src[name][tail_start:n_samples], compression='gzip')
            dst.attrs.update(dict(src.attrs))
        t_probs, t_feats = run_pass(tail_path, models)
    finally:
        Path(tail_path).unlink(missing_ok=True)
    primary_end = float(stops.max())
    t_starts = tail_start / TARGET_HZ + np.arange(len(t_probs)) * float(STEP_S)
    t_stops = t_starts + STEP_S
    keep = t_stops > primary_end
    t_starts = np.maximum(t_starts[keep], primary_end)
    return (np.concatenate([starts, t_starts]), np.concatenate([stops, t_stops[keep]]),
            np.vstack([probs, t_probs[keep]]), np.vstack([feats, t_feats[keep]]),
            np.concatenate([source, np.array(['tail'] * int(keep.sum()), dtype=object)]))


def to_epochs(starts, stops, probs, n_epochs, epoch_s=30.0):
    """5-s probabilities -> one row per 30-s epoch from the recording start (overlap-weighted mean, as Score_EDF
    aggregates them). Epochs with no 5-s step get NaN."""
    out = np.full((n_epochs, probs.shape[1]), np.nan)
    for i in range(n_epochs):
        overlap = np.maximum(0.0, np.minimum(stops, (i + 1) * epoch_s) - np.maximum(starts, i * epoch_s))
        total = overlap.sum()
        if total > 0:
            out[i] = (probs * overlap[:, None]).sum(axis=0) / total
    return out


def score_signals(bas, emg, others, sfreq, root, device='cpu', tail_context_chunks=TAIL_CONTEXT_CHUNKS,
                  workdir=None):
    """Whole SleepFM scoring of one night from its signals (see write_montage). The prepared 128-Hz file is a
    temporary file, deleted at the end. Returns a dict: starts_s, stops_s, probs (5-s), features, source,
    fc_weight, fc_bias, channels (name -> modality)."""
    models = load_models(root, device)
    groups = models['groups']
    channels = {name: modality_of(name, groups) for name in list(bas) + (['EMG'] if emg is not None else [])
                + list(others)}
    unknown = [n for n, m in channels.items() if m is None]
    if unknown:
        raise ValueError(f'channel name(s) unknown to SleepFM: {unknown}')
    missing = [m for m in MODALITIES if m not in channels.values()]
    if missing:
        raise ValueError(f'SleepFM needs at least one channel per modality; none for {missing}')
    tmp = Path(tempfile.mkdtemp(prefix='sleepfm_', dir=workdir))
    try:
        path = write_montage(tmp / 'montage_128hz.hdf5', bas, emg, others, sfreq)
        starts, stops, probs, feats, source = score_montage(path, models, tail_context_chunks)
    finally:
        for f in tmp.glob('*'):
            f.unlink(missing_ok=True)
        tmp.rmdir()
    weight, bias = head_parameters(models)
    return {'starts_s': starts, 'stops_s': stops, 'probs': probs, 'features': feats, 'source': source,
            'fc_weight': weight, 'fc_bias': bias, 'channels': channels}
