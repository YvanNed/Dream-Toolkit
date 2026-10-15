from __future__ import annotations

import json
from pathlib import Path
from typing import Sequence

import h5py
import mne
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt, spectrogram

from .io import NapWindow


EEG_DERIVATIONS: dict[str, tuple[str, str]] = {
    "F3-M2": ("F3", "M2"),
    "F4-M1": ("F4", "M1"),
    "C3-M2": ("C3", "M2"),
    "C4-M1": ("C4", "M1"),
    "O1-M2": ("O1", "M2"),
    "O2-M1": ("O2", "M1"),
}

EEG_BANDS: dict[str, tuple[float, float]] = {
    "Delta": (0.5, 4.0),
    "Theta": (4.0, 8.0),
    "Alpha": (8.0, 12.0),
    "Beta": (12.0, 30.0),
    "Gamma": (30.0, 45.0),
}

SLEEPFM_TARGET_HZ = 128

# Fixed SleepFM input filters (AASM-recommended; also the Compumedics display
# filters in the XML): EEG and EOG 0.3-35 Hz, chin EMG 10-100 Hz. Applied at the
# native rate before SleepFM's own anti-aliasing, 128-Hz interpolation and
# z-scoring. Without them, 50-Hz mains and other >35-Hz content pass SleepFM's
# 64-Hz anti-alias filter and the model assigns almost no probability to REM.
SLEEPFM_BANDPASS: dict[str, tuple[float, float]] = {"BAS": (0.3, 35.0), "EMG": (10.0, 100.0)}


def lagged_cross_correlation(
    reference: Sequence[float],
    response: Sequence[float],
    max_lag: int,
) -> pd.DataFrame:
    """Pearson correlation of reference[t] with response[t + lag].

    Positive lags mean that the response follows the reference. Each lag uses
    only pairs that overlap; no circular wrapping or interpolation is applied.
    """
    x = np.asarray(reference, dtype=float)
    y = np.asarray(response, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y):
        raise ValueError("reference and response must be one-dimensional and equally long")
    if max_lag < 0 or max_lag >= len(x):
        raise ValueError("max_lag must be between zero and len(reference) - 1")

    rows = []
    for lag in range(-max_lag, max_lag + 1):
        start = max(0, -lag)
        stop = min(len(x), len(x) - lag)
        x_overlap = x[start:stop]
        y_overlap = y[start + lag : stop + lag]
        valid = np.isfinite(x_overlap) & np.isfinite(y_overlap)
        correlation = (
            np.corrcoef(x_overlap[valid], y_overlap[valid])[0, 1]
            if valid.sum() >= 3
            and np.std(x_overlap[valid]) > 0
            and np.std(y_overlap[valid]) > 0
            else np.nan
        )
        rows.append({"lag_epochs": lag, "correlation": correlation, "n_pairs": int(valid.sum())})
    return pd.DataFrame(rows)


def _load_edf_channels(edf_path: str | Path) -> tuple[dict[str, np.ndarray], float]:
    """Load EDF channels in SI units, matching SleepFM's current MNE reader."""
    raw = mne.io.read_raw_edf(edf_path, preload=True, verbose="ERROR")
    return {
        name: raw.get_data(picks=[name])[0].astype(float, copy=False)
        for name in raw.ch_names
    }, float(raw.info["sfreq"])


def referenced_eeg(
    edf_path: str | Path,
    derivations: Sequence[str] = tuple(EEG_DERIVATIONS),
) -> tuple[dict[str, np.ndarray], float]:
    """Return requested EEG derivations in volts on MNE's common time grid."""
    signals, sfreq = _load_edf_channels(edf_path)
    unknown = sorted(set(derivations) - set(EEG_DERIVATIONS))
    if unknown:
        raise ValueError(f"Unknown EEG derivations: {unknown}")
    output = {}
    for name in derivations:
        active, reference = EEG_DERIVATIONS[name]
        output[name] = signals[active] - signals[reference]
    return output, sfreq


def bandpass(signal: np.ndarray, sfreq: float, low: float, high: float) -> np.ndarray:
    """Zero-phase 4th-order Butterworth band-pass; ``high`` is capped below Nyquist."""
    high = min(high, 0.99 * sfreq / 2.0)
    sos = butter(4, [low, high], btype="band", fs=sfreq, output="sos")
    return sosfiltfilt(sos, signal)


def _anti_alias_resample_standardize(signal: np.ndarray, sfreq: float) -> np.ndarray:
    """Reproduce SleepFM's 128-Hz interpolation and global standardization."""
    if sfreq > SLEEPFM_TARGET_HZ:
        nyquist = sfreq / 2.0
        sos = butter(4, (SLEEPFM_TARGET_HZ / 2.0) / nyquist, btype="low", output="sos")
        signal = sosfiltfilt(sos, signal)
    duration = len(signal) / sfreq
    old_time = np.linspace(0, duration, len(signal), endpoint=False)
    new_time = np.linspace(
        0, duration, int(duration * SLEEPFM_TARGET_HZ), endpoint=False
    )
    resampled = np.interp(new_time, old_time, signal)
    std = resampled.std()
    return (resampled - resampled.mean()) / std if std > 0 else resampled - resampled.mean()


def prepare_sleepfm_montage(
    edf_path: str | Path,
    output_path: str | Path,
    eeg_derivations: Sequence[str] = tuple(EEG_DERIVATIONS),
    include_eog: bool = True,
    force: bool = False,
    filters: dict[str, tuple[float, float]] | None = SLEEPFM_BANDPASS,
) -> Path:
    """Create a SleepFM HDF5 with explicit Compumedics channel semantics.

    EEG references follow the companion XML. EOG names are mapped to names
    recognized by SleepFM; Pleth is deliberately excluded from EMG. ``filters``
    band-passes the BAS (EEG + EOG) and EMG channels first; ``None`` gives the
    unfiltered official preprocessing.
    """
    output_path = Path(output_path)
    if output_path.exists() and not force:
        with h5py.File(output_path, "r") as handle:
            cached = handle.attrs.get("bandpass_hz", "{}")  # older caches were unfiltered
        if json.loads(cached) != json.loads(json.dumps(filters or {})):
            raise ValueError(
                f"{output_path} was prepared with filters {cached}; use force=True or another cache_tag"
            )
        return output_path

    signals, sfreq = _load_edf_channels(edf_path)
    selected: dict[str, np.ndarray] = {}
    for name in eeg_derivations:
        if name not in EEG_DERIVATIONS:
            raise ValueError(f"Unknown EEG derivation: {name}")
        active, reference = EEG_DERIVATIONS[name]
        selected[name.replace("-", ":")] = signals[active] - signals[reference]

    if include_eog:
        selected["EOG1"] = signals["vEOG"]
        selected["EOG2"] = signals["hEOG"]

    emg = signals["EMG"]
    if filters:
        selected = {name: bandpass(values, sfreq, *filters["BAS"]) for name, values in selected.items()}
        emg = bandpass(emg, sfreq, *filters["EMG"])
    selected.update(
        {
            "ECG": signals["ECG"],
            "EMG": emg,
            "Thoracic": signals["Thoracic"],
            "Abdominal": signals["Abdominal"],
            "Airflow": signals["Flux"],
            "SpO2": signals["SpO2"],
            "Pulse": signals["Pulse"],
        }
    )
    processed = {
        name: _anti_alias_resample_standardize(values, sfreq)
        for name, values in selected.items()
    }
    if len({len(values) for values in processed.values()}) != 1:
        raise ValueError("SleepFM montage channels have inconsistent lengths")
    if not all(np.isfinite(values).all() for values in processed.values()):
        raise ValueError("SleepFM montage contains non-finite values")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_path, "w") as handle:
        for name, values in processed.items():
            handle.create_dataset(
                name, data=values.astype("float32"), compression="gzip"
            )
        handle.attrs["eeg_derivations"] = json.dumps(list(eeg_derivations))
        handle.attrs["include_eog"] = include_eog
        handle.attrs["sampling_frequency_hz"] = SLEEPFM_TARGET_HZ
        handle.attrs["bandpass_hz"] = json.dumps(filters or {})
    return output_path


def infer_sleepfm_montage(
    edf_path: str | Path,
    output_dir: str | Path,
    eeg_derivations: Sequence[str] = tuple(EEG_DERIVATIONS),
    cache_tag: str = "aasm_all_eeg",
    include_eog: bool = True,
    device: str = "cpu",
    force: bool = False,
    cover_tail: bool = False,
) -> pd.DataFrame:
    """Prepare an explicit montage and run the official SleepFM checkpoints.

    Official SleepFM drops the final incomplete 5-minute encoder chunk. With
    ``cover_tail=True`` a second pass is run on the same standardized HDF5,
    cropped at its start so that its chunks end at the end of the recording;
    only its predictions after the primary pass are appended (``source`` =
    ``"tail"``), with the first appended step clipped to the primary end.
    """
    from .sleepfm import infer_hypnodensity

    output_dir = Path(output_dir)
    hdf5_path = output_dir / f"{Path(edf_path).stem}_{cache_tag}_128hz.hdf5"
    prepare_sleepfm_montage(
        edf_path,
        hdf5_path,
        eeg_derivations=eeg_derivations,
        include_eog=include_eog,
        force=force,
    )
    primary = infer_hypnodensity(
        edf_path,
        output_dir=output_dir,
        device=device,
        force=force,
        prepared_hdf5_path=hdf5_path,
        cache_tag=cache_tag,
    )
    if not cover_tail:
        return primary

    primary = primary.assign(source="primary")
    primary_end = float(primary["stop_seconds"].max())
    with h5py.File(hdf5_path, "r") as handle:
        n_samples = min(len(handle[name]) for name in handle)
    # Samples before the last whole 5-minute chunk that ends at the recording end.
    offset = n_samples % (5 * 60 * SLEEPFM_TARGET_HZ)
    if offset < 5 * SLEEPFM_TARGET_HZ:
        return primary

    tail_path = output_dir / f"{Path(edf_path).stem}_{cache_tag}_tail_128hz.hdf5"
    if force or not tail_path.exists():
        with h5py.File(hdf5_path, "r") as source, h5py.File(tail_path, "w") as target:
            for name in source:
                target.create_dataset(name, data=source[name][offset:n_samples], compression="gzip")
            target.attrs.update(dict(source.attrs))
            target.attrs["start_offset_samples"] = offset
    tail = infer_hypnodensity(
        edf_path,
        output_dir=output_dir,
        device=device,
        force=force,
        prepared_hdf5_path=tail_path,
        cache_tag=f"{cache_tag}_tail",
    )
    shift = offset / SLEEPFM_TARGET_HZ
    tail = tail.assign(
        start_seconds=tail["start_seconds"] + shift,
        stop_seconds=tail["stop_seconds"] + shift,
        source="tail",
    )
    tail = tail.loc[tail["stop_seconds"] > primary_end].copy()
    tail["start_seconds"] = tail["start_seconds"].clip(lower=primary_end)
    return pd.concat([primary, tail], ignore_index=True)


def normalized_bandpower_by_epoch(
    edf_path: str | Path,
    epochs: pd.DataFrame,
    derivations: Sequence[str] = tuple(EEG_DERIVATIONS),
    bands: dict[str, tuple[float, float]] = EEG_BANDS,
) -> pd.DataFrame:
    """Estimate relative Welch bandpower for every epoch and EEG derivation."""
    eeg, sfreq = referenced_eeg(edf_path, derivations)
    names = list(eeg)
    data = np.stack([eeg[name] for name in names]) * 1e6  # microvolts
    rows: list[dict[str, float | int | str]] = []
    for epoch in epochs.itertuples():
        start = max(0, int(round(epoch.start_seconds * sfreq)))
        stop = min(data.shape[1], int(round(epoch.stop_seconds * sfreq)))
        segment = data[:, start:stop]
        n_per_seg = min(int(round(4 * sfreq)), segment.shape[1])
        if n_per_seg < int(sfreq):
            continue
        psd, frequencies = mne.time_frequency.psd_array_welch(
            segment,
            sfreq=sfreq,
            fmin=min(low for low, _ in bands.values()),
            fmax=max(high for _, high in bands.values()),
            n_fft=n_per_seg,
            n_per_seg=n_per_seg,
            n_overlap=n_per_seg // 2,
            average="mean",
            verbose=False,
        )
        for channel_index, channel in enumerate(names):
            row: dict[str, float | int | str] = {
                "epoch": int(epoch.epoch),
                "channel": channel,
            }
            band_power = {}
            final_high = max(high for _, high in bands.values())
            for band, (low, high) in bands.items():
                # Welch bins are equally spaced. Half-open adjacent intervals
                # partition 0.5--45 Hz without double-counting boundaries.
                mask = (frequencies >= low) & (
                    frequencies <= high if high == final_high else frequencies < high
                )
                band_power[band] = psd[channel_index, mask].sum()
            total = sum(band_power.values())
            for band, power in band_power.items():
                row[band] = float(power / total)
            rows.append(row)
    return pd.DataFrame(rows)


def nap_spectrogram(
    edf_path: str | Path,
    window: NapWindow,
    derivations: Sequence[str] = tuple(EEG_DERIVATIONS),
    segment_seconds: float = 5.0,
    overlap_fraction: float = 0.5,
    fmin: float = 0.5,
    fmax: float = 35.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return channel-mean nap PSD as dB microvolt squared/Hz."""
    eeg, sfreq = referenced_eeg(edf_path, derivations)
    start = int(round(window.lights_off_seconds * sfreq))
    stop = int(round(window.reveil_seconds * sfreq))
    data = np.stack([values[start:stop] for values in eeg.values()]) * 1e6
    n_per_seg = int(round(segment_seconds * sfreq))
    frequencies, times, psd = spectrogram(
        data,
        fs=sfreq,
        window="hann",
        nperseg=n_per_seg,
        noverlap=int(round(n_per_seg * overlap_fraction)),
        detrend="constant",
        scaling="density",
        mode="psd",
        axis=-1,
    )
    keep = (frequencies >= fmin) & (frequencies <= fmax)
    mean_psd = psd[:, keep, :].mean(axis=0)
    return times / 60.0, frequencies[keep], 10 * np.log10(mean_psd + np.finfo(float).tiny)
