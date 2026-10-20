"""Scoring file readers, EDF-scoring file matching, and analysis window detection."""
from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

from sleephd.io import NapWindow, STAGE_TO_NUMBER, STAGE_NAMES, read_edf_start, read_hypnogram


XML_STAGE_MAP: dict[int, tuple[float, str]] = {
    0: (0, "W"),
    1: (1, "1"),
    2: (2, "2"),
    3: (3, "3"),
    5: (4, "R"),
    9: (np.nan, "?"),
}

_CSV_STAGE_COLUMNS = {"stage", "sleep_stage", "score", "hypnogram", "sleep stage"}


def read_scoring_xml(xml_path: str | Path) -> tuple[pd.DataFrame, float]:
    """Read sleep stages and epoch length from a Compumedics .edf.XML file."""
    root = ET.parse(xml_path).getroot()

    epoch_length_el = root.findtext("EpochLength")
    epoch_seconds = float(epoch_length_el) if epoch_length_el is not None else 30.0

    stages_el = root.find("SleepStages")
    if stages_el is None:
        raise ValueError(f"No <SleepStages> element in {xml_path}")
    codes = [int(el.text) for el in stages_el.findall("SleepStage")]
    if not codes:
        raise ValueError(f"<SleepStages> is empty in {xml_path}")

    rows = []
    for i, code in enumerate(codes):
        number, token = XML_STAGE_MAP.get(code, (np.nan, "?"))
        rows.append({
            "source_epoch": i,
            "stage_token": token,
            "stage_number": number,
            "stage": STAGE_NAMES.get(int(number), "Unknown") if pd.notna(number) else "Unknown",
        })
    df = pd.DataFrame(rows)
    df["stage_number"] = df["stage_number"].astype("Float64")
    return df, epoch_seconds


def read_scoring_txt(txt_path: str | Path) -> tuple[pd.DataFrame, float]:
    """Read a Compumedics one-token-per-line hypnogram export."""
    return read_hypnogram(txt_path), 30.0


def read_scoring_csv(csv_path: str | Path) -> tuple[pd.DataFrame, float]:
    """Read a CSV with a stage/score column."""
    raw = pd.read_csv(csv_path)
    raw.columns = [c.strip() for c in raw.columns]

    stage_col = None
    for col in raw.columns:
        if col.lower() in _CSV_STAGE_COLUMNS:
            stage_col = col
            break
    if stage_col is None:
        raise ValueError(
            f"Cannot find a stage column in {csv_path}. "
            f"Expected one of {sorted(_CSV_STAGE_COLUMNS)}; found {list(raw.columns)}"
        )

    tokens = raw[stage_col].astype(str).str.strip().str.upper()
    numbers = [STAGE_TO_NUMBER.get(t, np.nan) for t in tokens]

    df = pd.DataFrame({
        "source_epoch": np.arange(len(tokens), dtype=int),
        "stage_token": tokens.values,
        "stage_number": pd.array(numbers, dtype="Float64"),
        "stage": [STAGE_NAMES.get(int(x), "Unknown") if pd.notna(x) else "Unknown" for x in numbers],
    })
    return df, 30.0


def detect_and_read_scoring(path: str | Path) -> tuple[pd.DataFrame, float]:
    """Dispatch to the right reader based on file extension."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".xml":
        return read_scoring_xml(path)
    if suffix == ".txt":
        return read_scoring_txt(path)
    if suffix == ".csv":
        return read_scoring_csv(path)
    raise ValueError(f"Unsupported scoring file format: {path.name} (expected .xml, .txt, or .csv)")


def _extract_subject_night(name: str) -> tuple[str, str] | None:
    """Extract (subject_id, night_id) from a filename.

    subject_id = first 4 characters, night_id = 2 characters after the first '_'.
    Returns None if the name is too short or has no '_'.
    """
    stem = Path(name).stem
    # Strip a trailing .edf for double-extension names like "foo.edf.XML"
    if stem.lower().endswith(".edf"):
        stem = stem[:-4]
    if len(stem) < 4 or "_" not in stem:
        return None
    subject = stem[:4]
    after_us = stem.split("_", 1)[1]
    if len(after_us) < 2:
        return None
    night = after_us[:2]
    return (subject.lower(), night.lower())


def match_files(
    edf_folder: str | Path,
    scoring_folder: str | Path,
) -> tuple[list[dict], list[Path], list[Path]]:
    """Match EDF files with scoring files by subject ID and night ID.

    Subject ID = first 4 characters of the filename.
    Night ID   = 2 characters immediately after the first '_'.

    Returns (matched, unmatched_edfs, unmatched_scorings).
    Each matched entry: {"edf": Path, "scoring": Path, "stem": str}.
    """
    edf_folder = Path(edf_folder)
    scoring_folder = Path(scoring_folder)

    edfs = sorted(
        f for f in edf_folder.rglob("*.edf")
        if not f.name.startswith("._")
    )
    scorings = sorted(
        f for f in scoring_folder.rglob("*")
        if f.suffix.lower() in (".xml", ".txt", ".csv") and not f.name.startswith("._")
    )

    scoring_by_key: dict[tuple[str, str], Path] = {}
    for sc in scorings:
        key = _extract_subject_night(sc.name)
        if key is not None:
            scoring_by_key[key] = sc

    matched: list[dict] = []
    used_scorings: set[Path] = set()

    for edf in edfs:
        key = _extract_subject_night(edf.name)
        if key is None:
            continue
        sc_path = scoring_by_key.get(key)
        if sc_path is not None and sc_path not in used_scorings:
            matched.append({"edf": edf, "scoring": sc_path, "stem": edf.stem})
            used_scorings.add(sc_path)

    unmatched_edfs = [e for e in edfs if not any(m["edf"] == e for m in matched)]
    unmatched_scorings = [s for s in scorings if s not in used_scorings]
    return matched, unmatched_edfs, unmatched_scorings


def window_from_scoring(
    hypnogram: pd.DataFrame,
    edf_path: str | Path,
    epoch_seconds: float = 30.0,
) -> NapWindow:
    """Derive the analysis window from the first and last scored epochs."""
    scored = hypnogram.loc[hypnogram["stage_number"].notna(), "source_epoch"]
    if scored.empty:
        raise ValueError("No scored epochs found in the hypnogram")

    first_epoch = int(scored.min())
    last_epoch = int(scored.max())
    lights_off = first_epoch * epoch_seconds
    reveil = (last_epoch + 1) * epoch_seconds

    edf_start = read_edf_start(edf_path)
    lo_dt = edf_start + pd.Timedelta(seconds=lights_off)
    rv_dt = edf_start + pd.Timedelta(seconds=reveil)

    return NapWindow(
        lights_off_seconds=lights_off,
        reveil_seconds=reveil,
        lights_off_clock=lo_dt.strftime("%H:%M:%S"),
        reveil_clock=rv_dt.strftime("%H:%M:%S"),
    )
