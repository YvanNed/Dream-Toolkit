"""GSSC and YASA sleep staging with each toolbox's documented inputs.

Both adapters return 30-s probability rows (`start_seconds`, `stop_seconds`
from the EDF start, then Wake/N1/N2/N3/REM), the format consumed by
`sleephd.alignment.build_epoch_table`, so every downstream analysis is shared
with SleepFM. Only the stager and its own documented preprocessing differ:

- The whole recording is passed, as both toolboxes are designed for complete
  PSGs. It is cropped only at the start (by LIGHTS OFF modulo 30 s) so that the
  toolbox epochs coincide with the LIGHTS OFF-anchored analysis grid.
- No SleepFM filtering is applied. GSSC filters to 0.3-30 Hz itself; YASA states
  that signals must not be filtered or transformed beforehand.

GSSC (Hanna & Flöel, v0.0.9, README): mainly trained on C3, C4, F3, F4 and on
HEOG (left - right EOG); VEOG is not recommended. By default every EEG x EOG
combination, including EEG-only and EOG-only, is inferred and each epoch takes
the combination with the lowest entropy.

YASA (Vallat & Walker 2021, v0.7): one central EEG referenced to the mastoid
(C4-M1 preferred), one EOG (preferably LOC-M2), one chin EMG. These EDFs hold
only hardware-bipolar `vEOG`/`hEOG`, so the horizontal `hEOG` is used as the
closest available EOG. Age/sex are unknown, so the classifier without
demographics is used.
"""
from __future__ import annotations

from contextlib import redirect_stdout
import io
import logging
from pathlib import Path
from typing import Sequence
import warnings

import numpy as np
import pandas as pd

from .io import NapWindow
from .signals import EEG_DERIVATIONS, _load_edf_channels

PROBABILITY_COLUMNS = ["Wake", "N1", "N2", "N3", "REM"]
GSSC_EEG = ("C3-M2", "C4-M1", "F3-M2", "F4-M1")
GSSC_EOG = ("hEOG",)
YASA_EEG, YASA_EOG, YASA_EMG = "C4-M1", "hEOG", "EMG"


def _grid_offset(window: NapWindow, epoch_seconds: float = 30.0) -> float:
    """Seconds to drop at the recording start so toolbox epochs align with LIGHTS OFF."""
    return window.lights_off_seconds % epoch_seconds


def _raw_array(edf_path: str | Path, channels: dict[str, str], start_seconds: float):
    """MNE Raw with referenced derivations (volts) and per-channel types.

    `channels` maps a derivation or EDF channel name to its MNE type.
    """
    import mne

    signals, sfreq = _load_edf_channels(edf_path)
    data = []
    for name in channels:
        if name in EEG_DERIVATIONS:
            active, reference = EEG_DERIVATIONS[name]
            data.append(signals[active] - signals[reference])
        else:
            data.append(signals[name])
    info = mne.create_info(list(channels), sfreq, ch_types=list(channels.values()))
    raw = mne.io.RawArray(np.stack(data), info, verbose="error")
    first = int(round(start_seconds * sfreq))
    return raw.crop(tmin=first / sfreq, verbose="error") if first else raw


def _rows(probabilities: np.ndarray, start_seconds: np.ndarray) -> pd.DataFrame:
    frame = pd.DataFrame(np.asarray(probabilities, float), columns=PROBABILITY_COLUMNS)
    frame.insert(0, "stop_seconds", start_seconds + 30.0)
    frame.insert(0, "start_seconds", start_seconds)
    return frame


def infer_gssc(
    edf_path: str | Path,
    window: NapWindow,
    eeg: Sequence[str] = GSSC_EEG,
    eog: Sequence[str] = GSSC_EOG,
    cache_path: str | Path | None = None,
    force: bool = False,
) -> pd.DataFrame:
    """GSSC hypnodensity with its default channel combinations and 0.3-30 Hz filter."""
    if cache_path is not None and Path(cache_path).exists() and not force:
        return pd.read_csv(cache_path)
    from gssc.infer import EEGInfer
    import gssc.infer as _gssc_mod
    import torch
    from scipy.special import softmax

    offset = _grid_offset(window)
    raw = _raw_array(edf_path, {**{c: "eeg" for c in eeg}, **{c: "eog" for c in eog}}, offset)

    # Patch loudest_vote to capture the raw logits before they are argmax'd.
    _captured = {}
    _orig_vote = _gssc_mod.loudest_vote
    def _capturing_vote(logits):
        _captured["logits"] = logits
        return _orig_vote(logits)
    _gssc_mod.loudest_vote = _capturing_vote

    with warnings.catch_warnings(), redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        model = EEGInfer(use_cuda=False)
        _, times = model.mne_infer(raw, eeg=list(eeg), eog=list(eog), filter=True)

    _gssc_mod.loudest_vote = _orig_vote

    raw_logits = _captured.get("logits")
    if raw_logits is not None and np.ndim(raw_logits) == 3:
        logits_t = torch.FloatTensor(np.array(raw_logits))
        entrs = torch.zeros(logits_t.shape[:2])
        loss_func = torch.nn.NLLLoss(reduction="none")
        for idx in range(len(logits_t)):
            targs = torch.argmax(logits_t[idx], dim=-1)
            entrs[idx] = loss_func(logits_t[idx], targs)
        min_inds = np.argmin(entrs.numpy(), axis=0)
        min_logits = np.array(raw_logits)[min_inds, np.arange(raw_logits.shape[1])]
        probabilities = softmax(min_logits, axis=-1)
    else:
        labels = _
        n = len(labels)
        probabilities = np.zeros((n, 5))
        probabilities[np.arange(n), labels.astype(int)] = 1.0

    result = _rows(probabilities, offset + np.asarray(times, float))
    if cache_path is not None:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(cache_path, index=False)
    return result


def infer_yasa(
    edf_path: str | Path,
    window: NapWindow,
    eeg: str = YASA_EEG,
    eog: str | None = YASA_EOG,
    emg: str | None = YASA_EMG,
    cache_path: str | Path | None = None,
    force: bool = False,
) -> pd.DataFrame:
    """YASA hypnodensity from unfiltered EEG/EOG/EMG (its classifier filters internally)."""
    if cache_path is not None and Path(cache_path).exists() and not force:
        return pd.read_csv(cache_path)
    import yasa

    offset = _grid_offset(window)
    channels = {eeg: "eeg", **({eog: "eog"} if eog else {}), **({emg: "emg"} if emg else {})}
    raw = _raw_array(edf_path, channels, offset)
    level = yasa.logger.level
    yasa.logger.setLevel(logging.ERROR)
    try:
        staging = yasa.SleepStaging(raw, eeg_name=eeg, eog_name=eog, emg_name=emg)
        staging.predict()
        probabilities = staging.predict_proba()
    finally:
        yasa.logger.setLevel(level)
    probabilities = probabilities.rename(columns={"W": "Wake", "WAKE": "Wake", "R": "REM"})
    result = _rows(probabilities[PROBABILITY_COLUMNS].to_numpy(), offset + 30.0 * np.arange(len(probabilities)))
    if cache_path is not None:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(cache_path, index=False)
    return result
