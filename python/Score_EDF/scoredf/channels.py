"""Channel alias detection, requirement checking, and transparent renaming."""
from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import Sequence

import mne


ALIASES: dict[str, list[str]] = {
    # EEG electrodes (raw names needed for derivations)
    "C3": ["C3", "EEG C3", "EEG C3-REF", "C3-REF", "EEG C3-LER", "EEG C3-LE"],
    "C4": ["C4", "EEG C4", "EEG C4-REF", "C4-REF", "EEG C4-LER", "EEG C4-LE"],
    "F3": ["F3", "EEG F3", "EEG F3-REF", "F3-REF", "EEG F3-LER", "EEG F3-LE"],
    "F4": ["F4", "EEG F4", "EEG F4-REF", "F4-REF", "EEG F4-LER", "EEG F4-LE"],
    "O1": ["O1", "EEG O1", "EEG O1-REF", "O1-REF", "EEG O1-LER", "EEG O1-LE"],
    "O2": ["O2", "EEG O2", "EEG O2-REF", "O2-REF", "EEG O2-LER", "EEG O2-LE"],
    "M1": ["M1", "A1", "EEG M1", "EEG A1", "M1-REF", "A1-REF", "EEG M1-REF", "EEG A1-REF"],
    "M2": ["M2", "A2", "EEG M2", "EEG A2", "M2-REF", "A2-REF", "EEG M2-REF", "EEG A2-REF"],
    # EOG
    "vEOG": ["vEOG", "EOG LOC", "LOC", "EOG1", "E1", "EOG E1", "EOG Left",
             "EOG LOC-REF", "EOG LOC-M2", "EOG1-REF", "LEOG", "EOG-L", "LOC-M2",
             "EOG LOC-LER"],
    "hEOG": ["hEOG", "EOG ROC", "ROC", "EOG2", "E2", "EOG E2", "EOG Right",
             "EOG ROC-REF", "EOG ROC-M1", "EOG2-REF", "REOG", "EOG-R", "ROC-M1",
             "EOG ROC-LER"],
    # EMG
    "EMG": ["EMG", "EMG Chin", "EMG1", "Chin", "EMG chin", "EMG Chin1",
            "EMG Chin1-Chin2", "Chin1-Chin2", "EMG1-EMG2", "EMG-REF",
            "EMG Chin-REF"],
    # ECG
    "ECG": ["ECG", "EKG", "ECG1", "ECG I", "ECG-REF", "EKG-REF", "ECG II",
            "ECG1-ECG2"],
    # Respiratory (for SleepFM)
    "Thoracic": ["Thoracic", "Thor", "Chest", "Thor Effort", "Thorax",
                 "Effort THO", "RIP Thoracic", "THO", "Chest Effort"],
    "Abdominal": ["Abdominal", "Abdo", "Abd", "Abdomen", "Abd Effort",
                  "Effort ABD", "RIP Abdominal", "ABD", "Abdom Effort"],
    "Flux": ["Flux", "Airflow", "Flow", "Nasal", "NasalFlow", "Nasal Pressure",
             "Nasal Flow", "Thermistor", "AIRFLOW", "Airflow-0"],
    "SpO2": ["SpO2", "SaO2", "Oximetry", "SAO2", "SPO2", "Ox Status"],
    "Pulse": ["Pulse", "Pleth", "PPG", "Pulse Rate", "HR", "PULSE"],
}

STAGER_REQUIREMENTS: dict[str, dict[str, list[str]]] = {
    "YASA": {
        "required": ["C4", "M1"],
        "optional": ["hEOG", "EMG"],
    },
    "GSSC": {
        "required": ["hEOG"],
        "optional": [],
        "any_of": [["C3", "M2"], ["C4", "M1"], ["F3", "M2"], ["F4", "M1"]],
    },
    "SleepFM": {
        "required": [
            "C3", "C4", "F3", "F4", "O1", "O2", "M1", "M2",
            "vEOG", "hEOG", "EMG", "ECG",
            "Thoracic", "Abdominal", "Flux", "SpO2", "Pulse",
        ],
        "optional": [],
    },
}
# Recal A/B use SleepFM output, so same channel requirements
STAGER_REQUIREMENTS["SleepFM recal A"] = STAGER_REQUIREMENTS["SleepFM"]
STAGER_REQUIREMENTS["SleepFM recal B"] = STAGER_REQUIREMENTS["SleepFM"]


def detect_channels(edf_path: str | Path) -> dict[str, str | None]:
    """Auto-detect standard channel roles from an EDF's channel names."""
    raw = mne.io.read_raw_edf(str(edf_path), preload=False, verbose="ERROR")
    edf_names = raw.ch_names
    del raw

    lower_to_original = {}
    for name in edf_names:
        lower_to_original[name.strip().lower()] = name

    mapping: dict[str, str | None] = {}
    used: set[str] = set()

    for standard, aliases in ALIASES.items():
        found = None
        for alias in aliases:
            key = alias.strip().lower()
            if key in lower_to_original and lower_to_original[key] not in used:
                found = lower_to_original[key]
                used.add(found)
                break
        mapping[standard] = found

    return mapping


def get_edf_channel_names(edf_path: str | Path) -> list[str]:
    """Return all channel names from an EDF file."""
    raw = mne.io.read_raw_edf(str(edf_path), preload=False, verbose="ERROR")
    names = list(raw.ch_names)
    del raw
    return names


def check_stager_requirements(
    mapping: dict[str, str | None],
    stager: str,
) -> tuple[bool, list[str]]:
    """Check whether the channel mapping satisfies a stager's requirements."""
    reqs = STAGER_REQUIREMENTS.get(stager)
    if reqs is None:
        return True, []

    missing = [ch for ch in reqs.get("required", []) if mapping.get(ch) is None]

    any_of = reqs.get("any_of", [])
    if any_of:
        has_any = any(
            all(mapping.get(ch) is not None for ch in group)
            for group in any_of
        )
        if not has_any:
            missing.append(f"at least one EEG derivation pair from {any_of}")

    return len(missing) == 0, missing


# ---------------------------------------------------------------------------
# Channel renaming via monkeypatch
# ---------------------------------------------------------------------------
_original_load = None
_current_mapping: dict[str, str] | None = None


def _patched_load_edf_channels(edf_path):
    """Wrapper that renames channel dict keys to standard names."""
    signals, sfreq = _original_load(edf_path)
    if _current_mapping is None:
        return signals, sfreq

    reverse = {edf_name: standard for standard, edf_name in _current_mapping.items() if edf_name is not None}
    renamed = {}
    for name, data in signals.items():
        standard = reverse.get(name, name)
        renamed[standard] = data
    return renamed, sfreq


@contextmanager
def channel_mapping(mapping: dict[str, str | None]):
    """Context manager that patches channel loading to use the given mapping."""
    global _original_load, _current_mapping
    import sleephd.signals as _signals

    _original_load = _signals._load_edf_channels
    _current_mapping = {k: v for k, v in mapping.items() if v is not None}
    _signals._load_edf_channels = _patched_load_edf_channels
    try:
        yield
    finally:
        _signals._load_edf_channels = _original_load
        _original_load = None
        _current_mapping = None
