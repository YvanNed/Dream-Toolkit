"""
Generate synthetic EDF test files with controlled defects.

Derives test fixtures from a clean baseline recording (73.edf) by injecting
one controlled defect per output file. A "combined" file mixes defects on
multiple channels for integration testing.

Each defect file is given a unique participant-like ID (731–738) so that
2_select&remap_channels_edf does not group them all under participant "73".
The mapping is defined in DEFECT_IDS below.

The channels kept in each output are exactly those listed for the base
participant in config_param/remap_reref_persubject.json — so the test files
remain compatible with the standard remap/load workflow.

This script also writes a config entry for each generated file into
config_param/remap_reref_persubject.json (copied from the "73" entry),
so the quality_overview notebook can process them without manual edits.

A manifest TSV (test_data_manifest.tsv) records each generated file with
ground truth (defect type, channel, parameters, expected detection metric).

Usage:
    python tools/generate_test_data.py              # idempotent: skip existing
    python tools/generate_test_data.py --force      # regenerate all
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import edfio
import mne
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TEST_DATA = ROOT / "tools" / "test_data"
CONFIG = TEST_DATA / "config_param" / "remap_reref_persubject.json"
BASE = "73"
SEED = 42
LINE_FREQ = 50.0

# Unique participant-like ID for each defect type.
# Using 731–738 avoids grouping all defect files under participant "73"
# in 2_select&remap_channels_edf (which extracts the ID from the filename stem).
DEFECT_IDS = {
    "clipping":        "731",
    "combined":        "732",
    "dc_drift":        "733",
    "dead_channel":    "734",
    "flat_segment":    "735",
    "line_noise":      "736",
    "movement_bursts": "737",
    "quantization":    "738",
    "not_anonymized":  "739",
}


def load_base_channels() -> list[str]:
    """Return the source channels to keep for the base file, per remap config."""
    cfg = json.loads(CONFIG.read_text(encoding="utf-8"))
    return list(cfg[BASE]["remap"].keys())


def load_base_raw(channels: list[str]) -> mne.io.BaseRaw:
    fp = TEST_DATA / f"{BASE}.edf"
    # include= avoids loading non-EEG channels (e.g. DHR at 512 Hz), which would
    # cause MNE to upsample all channels to 512 Hz and corrupt sfreq/filter headers.
    raw = mne.io.read_raw_edf(fp, preload=True, include=channels, verbose="ERROR")
    return raw


def read_original_phys_bounds(channels: list[str]) -> dict[str, tuple[float, float]]:
    """Read physical_min/max (in µV) directly from the 73.edf binary header."""
    with open(TEST_DATA / f"{BASE}.edf", "rb") as f:
        hdr = f.read(256)
        n = int(hdr[252:256].decode("ascii").strip())
        per_ch = f.read(n * (16 + 80 + 8 + 8 + 8 + 8 + 8 + 80 + 8 + 32))
    labels_raw  = per_ch[0 : 16 * n]
    phmin_start = 16 * n + 80 * n + 8 * n      # after label + transducer + units
    phmax_start = phmin_start + 8 * n
    phmin_raw   = per_ch[phmin_start : phmin_start + 8 * n]
    phmax_raw   = per_ch[phmax_start : phmax_start + 8 * n]
    bounds = {}
    for i in range(n):
        label = labels_raw[16 * i : 16 * (i + 1)].decode("ascii").strip()
        if label in channels:
            bounds[label] = (
                float(phmin_raw[8 * i : 8 * (i + 1)].decode("ascii").strip()),
                float(phmax_raw[8 * i : 8 * (i + 1)].decode("ascii").strip()),
            )
    return bounds


def copy_companions(file_id: str) -> None:
    """Copy hypnogram, XML, and event CSV next to the new EDF."""
    companions = [
        (f"{BASE}_Hypnogram_Export.txt",    f"{file_id}_Hypnogram_Export.txt"),
        (f"{BASE}.edf.XML",                  f"{file_id}.edf.XML"),
        (f"{BASE}_event_xml.csv",            f"{file_id}_event_xml.csv"),
    ]
    for src_name, dst_name in companions:
        src = TEST_DATA / src_name
        dst = TEST_DATA / dst_name
        if src.exists():
            shutil.copy(src, dst)


def export_edf(raw: mne.io.BaseRaw, file_id: str,
               phys_overrides: dict[str, tuple[float, float]] | None = None) -> Path:
    """Export raw to EDF using physical bounds from 73.edf header.

    phys_overrides: per-channel overrides {ch_name: (phys_min_uV, phys_max_uV)}.
    Use when a defect changes the effective signal range (e.g. clipping to ±75 µV).
    """
    out = TEST_DATA / f"{file_id}.edf"
    base_bounds = read_original_phys_bounds(raw.ch_names)
    signals = []
    for ch in raw.ch_names:
        data_uV = raw.get_data(picks=[ch])[0] * 1e6
        phys_range = (phys_overrides or {}).get(ch, base_bounds[ch])
        # clip to physical range so no sample falls outside the declared EDF bounds
        data_clipped = np.clip(data_uV, phys_range[0], phys_range[1])
        signals.append(edfio.EdfSignal(
            data_clipped,
            sampling_frequency=raw.info["sfreq"],
            label=ch,
            physical_dimension="uV",
            physical_range=phys_range,
        ))
    edfio.Edf(signals, data_record_duration=1).write(str(out))
    copy_companions(file_id)
    return out


def update_config(file_id: str) -> None:
    """Add a config entry for file_id (copied from BASE) if not already present."""
    cfg = json.loads(CONFIG.read_text(encoding="utf-8"))
    if file_id not in cfg:
        cfg[file_id] = cfg[BASE]
        CONFIG.write_text(
            json.dumps(cfg, indent=2, ensure_ascii=False), encoding="utf-8"
        )


# Header patching -----------------------------------------------------------

def patch_patient_id(edf_path: Path, patient_id_str: str) -> None:
    """Overwrite the 80-byte Local Patient ID field in an EDF binary header.

    Bytes 8–87 (0-indexed) hold this field. Signal data is never touched.
    The string is ASCII-encoded, truncated to 80 chars, and right-padded with
    spaces — exactly the EDF spec layout.
    """
    field = patient_id_str.encode("ascii", errors="replace")[:80].ljust(80)
    with open(edf_path, "r+b") as f:
        f.seek(8)
        f.write(field)


# Defect injection ----------------------------------------------------------

def inject_flat_segment(raw, channel, start_s, duration_s,
                        residual_uv=0.01, seed=SEED):
    raw = raw.copy()
    sf = raw.info["sfreq"]
    a = int(start_s * sf)
    b = a + int(duration_s * sf)
    idx = raw.ch_names.index(channel)
    rng = np.random.default_rng(seed)
    raw._data[idx, a:b] = rng.normal(0, residual_uv * 1e-6, b - a)
    return raw


def inject_dead_channel(raw, channel, scale=0.01):
    raw = raw.copy()
    idx = raw.ch_names.index(channel)
    raw._data[idx, :] *= scale
    return raw


def inject_clipping(raw, channel, threshold_uv=75.0):
    raw = raw.copy()
    idx = raw.ch_names.index(channel)
    thr = threshold_uv * 1e-6
    np.clip(raw._data[idx], -thr, thr, out=raw._data[idx])
    return raw


def inject_dc_drift(raw, channel, shift_uv=50.0, center_s=None, width_s=1800):
    """Sigmoidal DC step (mimics sweat/polarization-induced baseline shift)."""
    raw = raw.copy()
    sf = raw.info["sfreq"]
    idx = raw.ch_names.index(channel)
    n = raw.n_times
    if center_s is None:
        center_s = raw.times[-1] / 2
    t = np.arange(n) / sf
    sig = 1.0 / (1.0 + np.exp(-6.0 * (t - center_s) / width_s))
    raw._data[idx] += shift_uv * 1e-6 * sig
    return raw


def inject_movement_bursts(raw, channel, n_bursts=10, amp_uv=300.0,
                           dur_s_range=(2.0, 5.0), seed=SEED):
    raw = raw.copy()
    sf = raw.info["sfreq"]
    idx = raw.ch_names.index(channel)
    n = raw.n_times
    rng = np.random.default_rng(seed)
    positions = rng.uniform(0.05 * n, 0.95 * n, n_bursts).astype(int)
    for pos in positions:
        dur = rng.uniform(*dur_s_range)
        length = int(dur * sf)
        burst = rng.normal(0, 1, length) * np.hanning(length)
        burst = burst / np.max(np.abs(burst)) * amp_uv * 1e-6
        end = min(pos + length, n)
        raw._data[idx, pos:end] += burst[: end - pos]
    return raw


def inject_line_noise(raw, channel, amp_uv_pp=50.0, freq=LINE_FREQ):
    raw = raw.copy()
    sf = raw.info["sfreq"]
    idx = raw.ch_names.index(channel)
    n = raw.n_times
    t = np.arange(n) / sf
    sine = np.sin(2 * np.pi * freq * t) * (amp_uv_pp / 2) * 1e-6
    raw._data[idx] += sine
    return raw


def inject_quantization(raw, channel, n_levels=16, signal_range_uv=400.0):
    """Coarse quantization → comb-like multimodal histogram."""
    raw = raw.copy()
    idx = raw.ch_names.index(channel)
    half = (signal_range_uv / 2) * 1e-6
    step = (2 * half) / n_levels
    raw._data[idx] = np.round(raw._data[idx] / step) * step
    return raw


def inject_combined(raw, mid_s):
    raw = inject_flat_segment(raw, "Fp1", start_s=mid_s - 900, duration_s=1800)
    raw = inject_dc_drift(raw, "C3", shift_uv=50.0, center_s=mid_s, width_s=1800)
    raw = inject_line_noise(raw, "O1", amp_uv_pp=50.0)
    return raw


# Driver --------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="Regenerate even if output files already exist")
    args = parser.parse_args()

    channels = load_base_channels()
    print(f"Base: {BASE}.edf — keeping channels: {channels}")

    raw = load_base_raw(channels)
    total_s = raw.times[-1]
    mid_s = total_s / 2
    print(f"  duration: {total_s/3600:.2f} h, sfreq: {raw.info['sfreq']} Hz")

    rows = []

    def already_done(file_id):
        return (TEST_DATA / f"{file_id}.edf").exists() and not args.force

    def add_row(file_id, channel, defect, params, expected_metric):
        rows.append({
            "file": f"{file_id}.edf",
            "base_file": f"{BASE}.edf",
            "defect": defect,
            "channel": channel,
            "params": json.dumps(params),
            "expected_metric": expected_metric,
        })

    # Each entry: (suffix, channel, defect_id, expected_metric, params, build_fn,
    #              phys_overrides, header_patch)
    # phys_overrides: {ch_name: (phys_min_uV, phys_max_uV)} or None (use original bounds).
    # header_patch: callable(Path) applied after export, or None. Used to inject
    #               header-only defects without touching signal data.
    _NOT_ANON_PID = "01016 M 15-MAR-1980 DUPONT_JEAN"
    plan = [
        ("flat_segment", "Fp1", "flat_segment_30min", "flat_pct_high",
         {"start_s": mid_s - 900, "duration_s": 1800, "residual_uv": 0.01},
         lambda r: inject_flat_segment(r, "Fp1", start_s=mid_s - 900, duration_s=1800),
         None, None),
        ("dead_channel", "C3", "dead_channel", "low_variance_high_kurtosis",
         {"scale": 0.01},
         lambda r: inject_dead_channel(r, "C3", scale=0.01),
         None, None),
        ("clipping", "Fp1", "clipping", "physical_bounds_pct_high",
         {"threshold_uv": 75.0},
         lambda r: inject_clipping(r, "Fp1", threshold_uv=75.0),
         {"Fp1": (-75.0, 75.0)}, None),   # Fp1 is clipped to ±75 µV
        ("dc_drift", "C3", "dc_drift_bimodal", "bimodal_distribution",
         {"shift_uv": 50.0, "center_s": mid_s, "width_s": 1800},
         lambda r: inject_dc_drift(r, "C3", shift_uv=50.0, center_s=mid_s, width_s=1800),
         None, None),
        ("movement_bursts", "Fp1", "movement_bursts", "peak_to_peak_rejection_step2",
         {"n_bursts": 10, "amp_uv": 300.0, "dur_s_range": [2, 5]},
         lambda r: inject_movement_bursts(r, "Fp1", n_bursts=10, amp_uv=300.0),
         None, None),
        ("line_noise", "O1", "line_noise_50Hz", "spectral_line_50hz",
         {"amp_uv_pp": 50.0, "freq": LINE_FREQ},
         lambda r: inject_line_noise(r, "O1", amp_uv_pp=50.0, freq=LINE_FREQ),
         None, None),
        ("quantization", "C3", "quantization_16levels", "multimodal_distribution",
         {"n_levels": 16, "signal_range_uv": 400.0},
         lambda r: inject_quantization(r, "C3", n_levels=16, signal_range_uv=400.0),
         None, None),
        ("combined", "Fp1+C3+O1", "combined_defects", "multiple_metrics_flagged",
         {"Fp1": "flat_segment", "C3": "dc_drift", "O1": "line_noise"},
         lambda r: inject_combined(r, mid_s),
         None, None),
        # Header-only defect: clean signal, non-anonymized patient_id in EDF header.
        # Expected: anonymization_check flags "header NOT anonymized (file name looks clean)"
        # (the file name 739_not_anonymized.edf does not contain the patient name DUPONT_JEAN).
        ("not_anonymized", "n/a", "non_anonymized_header", "anon_warning_header_not_anonymized",
         {"patient_id": _NOT_ANON_PID},
         lambda r: r.copy(),
         None, lambda p: patch_patient_id(p, _NOT_ANON_PID)),
    ]

    for i, (suffix, channel, defect_id, metric, params, build, phys_overrides, header_patch) in enumerate(plan, 1):
        file_id = f"{DEFECT_IDS[suffix]}_{suffix}"
        tag = f"[{i}/{len(plan)}] {file_id}"
        if already_done(file_id):
            print(f"{tag} — skip (exists)")
        else:
            print(f"{tag} — generating ...")
            modified = build(raw)
            out_path = export_edf(modified, file_id, phys_overrides)
            if header_patch is not None:
                header_patch(out_path)
                print(f"  header patched: patient_id set to {params.get('patient_id', '')!r}")
        update_config(file_id)
        add_row(file_id, channel, defect_id, params, metric)

    manifest_path = TEST_DATA / "test_data_manifest.tsv"
    pd.DataFrame(rows).to_csv(manifest_path, sep="\t", index=False)
    print(f"\nManifest written: {manifest_path}")
    print("Config JSON updated with entries for all generated files.")


if __name__ == "__main__":
    main()
