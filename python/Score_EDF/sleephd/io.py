from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd


STAGE_TO_NUMBER = {
    "W": 0,
    "WAKE": 0,
    "0": 0,
    "1": 1,
    "N1": 1,
    "S1": 1,
    "2": 2,
    "N2": 2,
    "S2": 2,
    "3": 3,
    "N3": 3,
    "S3": 3,
    "SWS": 3,
    "4": 3,
    "S4": 3,
    "R": 4,
    "REM": 4,
    "5": 4,
}
STAGE_NAMES = {0: "Wake", 1: "N1", 2: "N2", 3: "N3", 4: "REM"}


@dataclass(frozen=True)
class NapWindow:
    """Nap markers expressed in seconds from the EDF recording start."""

    lights_off_seconds: float
    reveil_seconds: float
    lights_off_clock: str
    reveil_clock: str

    @property
    def duration_seconds(self) -> float:
        return self.reveil_seconds - self.lights_off_seconds


def read_edf_start(edf_path: str | Path) -> datetime:
    """Read the EDF start timestamp without loading its signal data."""
    raw = Path(edf_path).read_bytes()[:256]
    if len(raw) < 256:
        raise ValueError(f"Invalid EDF header: {edf_path}")
    date_text = raw[168:176].decode("latin-1").strip()
    time_text = raw[176:184].decode("latin-1").strip()
    return datetime.strptime(f"{date_text} {time_text}", "%d.%m.%y %H.%M.%S")


def read_epoch_length(xml_path: str | Path) -> float:
    root = ET.parse(xml_path).getroot()
    value = root.findtext("EpochLength")
    if value is None:
        raise ValueError(f"EpochLength is missing from {xml_path}")
    return float(value)


def _clock_after_recording_start(clock: str, recording_start: datetime) -> datetime:
    parsed = datetime.strptime(clock.strip(), "%H:%M:%S").time()
    candidate = datetime.combine(recording_start.date(), parsed)
    if candidate < recording_start:
        candidate += timedelta(days=1)
    return candidate


def read_study_log(log_path: str | Path, edf_path: str | Path) -> tuple[pd.DataFrame, NapWindow]:
    """Read Compumedics markers and locate LIGHTS OFF / REVEIL."""
    log = pd.read_csv(
        log_path,
        header=None,
        names=["clock", "export_epoch", "code", "event"],
        dtype={"clock": str, "event": str},
    )
    log["event"] = log["event"].str.strip()
    start = read_edf_start(edf_path)
    log["timestamp"] = log["clock"].map(lambda value: _clock_after_recording_start(value, start))
    log["seconds_from_edf_start"] = (log["timestamp"] - start).dt.total_seconds()

    def marker_seconds(marker: str) -> tuple[float, str]:
        rows = log.loc[log["event"].str.upper() == marker]
        if len(rows) != 1:
            raise ValueError(f"Expected exactly one {marker!r} marker in {log_path}; found {len(rows)}")
        row = rows.iloc[0]
        return float(row["seconds_from_edf_start"]), str(row["clock"])

    lights_off, lights_clock = marker_seconds("LIGHTS OFF")
    reveil, reveil_clock = marker_seconds("REVEIL")
    if reveil <= lights_off:
        reveil += 24 * 60 * 60
    return log, NapWindow(lights_off, reveil, lights_clock, reveil_clock)


def read_hypnogram(path: str | Path) -> pd.DataFrame:
    """Read one Compumedics sleep-stage token per scoring epoch."""
    tokens = [line.strip().upper() for line in Path(path).read_text(encoding="utf-8-sig").splitlines()]
    numbers = [STAGE_TO_NUMBER.get(token, np.nan) for token in tokens]
    return pd.DataFrame(
        {
            "source_epoch": np.arange(len(tokens), dtype=int),
            "stage_token": tokens,
            "stage_number": pd.array(numbers, dtype="Float64"),
            "stage": [STAGE_NAMES.get(int(x)) if pd.notna(x) else "Unknown" for x in numbers],
        }
    )


def read_tet(path: str | Path) -> pd.DataFrame:
    """Read a WebPlotDigitizer-style TET workbook with a title row then X/Y headers."""
    frame = pd.read_excel(path, header=1, usecols="A:B")
    frame.columns = [str(column).strip().upper() for column in frame.columns]
    if not {"X", "Y"}.issubset(frame.columns):
        raise ValueError(f"Expected X and Y columns in {path}; found {list(frame.columns)}")
    frame = frame[["X", "Y"]].apply(pd.to_numeric, errors="coerce").dropna()
    frame = frame.groupby("X", as_index=False, sort=True)["Y"].mean()
    if len(frame) < 2:
        raise ValueError(f"TET needs at least two numeric points: {path}")
    return frame

