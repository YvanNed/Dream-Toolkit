from __future__ import annotations

import math

import numpy as np
import pandas as pd

from .io import NapWindow, STAGE_NAMES


def make_nap_epochs(window: NapWindow, epoch_seconds: float = 30.0) -> pd.DataFrame:
    """Create LIGHTS OFF-anchored epochs, retaining a final partial epoch."""
    count = math.ceil(window.duration_seconds / epoch_seconds)
    starts = window.lights_off_seconds + np.arange(count) * epoch_seconds
    stops = np.minimum(starts + epoch_seconds, window.reveil_seconds)
    return pd.DataFrame(
        {
            "epoch": np.arange(count, dtype=int),
            "start_seconds": starts,
            "stop_seconds": stops,
            "duration_seconds": stops - starts,
            "minutes_from_lights_off": (starts - window.lights_off_seconds) / 60,
        }
    )


def _piecewise_linear_mean(x: np.ndarray, y: np.ndarray, left: float, right: float) -> float:
    if right <= left:
        return np.nan
    inside = (x > left) & (x < right)
    points = np.concatenate(([left], x[inside], [right]))
    values = np.interp(points, x, y)
    return float(np.trapz(values, points) / (right - left))


def resample_tet_to_epochs(
    tet: pd.DataFrame,
    epochs: pd.DataFrame,
    window: NapWindow,
) -> pd.DataFrame:
    """Map normalized TET X anchors to clock time and average Y within each epoch."""
    clean = tet.loc[:, ["X", "Y"]].copy()
    # Digitized strokes may overshoot calibration axes slightly. Values outside
    # X0/X1 are not part of the nap; small Y overshoots are bounded by Y0/Y1.
    clean = clean.loc[clean["X"].between(0.0, 1.0)].copy()
    clean["Y"] = clean["Y"].clip(0.0, 1.0)
    clean = clean.groupby("X", as_index=False, sort=True)["Y"].mean()
    x_seconds = window.lights_off_seconds + clean["X"].to_numpy() * window.duration_seconds
    y = clean["Y"].to_numpy()
    values = [
        _piecewise_linear_mean(x_seconds, y, row.start_seconds, row.stop_seconds)
        for row in epochs.itertuples()
    ]
    result = epochs.copy()
    result["tet"] = values
    return result


def align_manual_hypnogram(
    hypnogram: pd.DataFrame,
    epochs: pd.DataFrame,
    source_epoch_seconds: float = 30.0,
) -> pd.DataFrame:
    """Assign the manual stage having the largest overlap with each nap epoch."""
    source_starts = hypnogram["source_epoch"].to_numpy(dtype=float) * source_epoch_seconds
    source_stops = source_starts + source_epoch_seconds
    source_values = hypnogram["stage_number"].to_numpy(dtype=float, na_value=np.nan)
    aligned: list[float] = []
    coverage: list[float] = []
    for row in epochs.itertuples():
        overlap = np.maximum(0.0, np.minimum(source_stops, row.stop_seconds) - np.maximum(source_starts, row.start_seconds))
        valid = (~np.isnan(source_values)) & (overlap > 0)
        coverage.append(float(overlap[valid].sum() / row.duration_seconds))
        if not valid.any():
            aligned.append(np.nan)
            continue
        votes = {stage: overlap[valid & (source_values == stage)].sum() for stage in range(5)}
        aligned.append(float(max(votes, key=votes.get)))
    result = epochs.copy()
    result["manual_stage_number"] = aligned
    result["manual_stage"] = [STAGE_NAMES.get(int(value), "Unknown") if pd.notna(value) else "Unknown" for value in aligned]
    result["manual_coverage"] = coverage
    return result


def aggregate_hypnodensity(
    five_second: pd.DataFrame,
    epochs: pd.DataFrame,
    source_step_seconds: float = 5.0,
) -> pd.DataFrame:
    """Overlap-weight SleepFM 5-s probabilities into arbitrary nap epochs."""
    probability_columns = ["Wake", "N1", "N2", "N3", "REM"]
    source_starts = five_second["start_seconds"].to_numpy(float)
    source_stops = five_second.get("stop_seconds", five_second["start_seconds"] + source_step_seconds).to_numpy(float)
    probabilities = five_second[probability_columns].to_numpy(float)
    rows = []
    for epoch in epochs.itertuples():
        overlap = np.maximum(0.0, np.minimum(source_stops, epoch.stop_seconds) - np.maximum(source_starts, epoch.start_seconds))
        total = overlap.sum()
        if total == 0:
            rows.append([np.nan] * len(probability_columns))
        else:
            rows.append((probabilities * overlap[:, None]).sum(axis=0) / total)
    result = epochs.copy()
    result[probability_columns] = rows
    covered = result[probability_columns].notna().all(axis=1)
    stage_number = pd.Series(pd.NA, index=result.index, dtype="Int64")
    stage_number.loc[covered] = result.loc[covered, probability_columns].to_numpy().argmax(axis=1)
    result["predicted_stage_number"] = stage_number
    result["predicted_stage"] = stage_number.map(STAGE_NAMES).fillna("Unknown")
    return result


def build_epoch_table(
    window: NapWindow,
    tet: pd.DataFrame | None,
    hypnogram: pd.DataFrame,
    hypnodensity_5s: pd.DataFrame | None = None,
    epoch_seconds: float = 30.0,
) -> pd.DataFrame:
    epochs = make_nap_epochs(window, epoch_seconds)
    if tet is None:
        # Recordings without a TET trace keep all PSG-derived columns.
        table = epochs.copy()
        table["tet"] = np.nan
    else:
        table = resample_tet_to_epochs(tet, epochs, window)
    manual = align_manual_hypnogram(hypnogram, epochs, epoch_seconds)
    table[["manual_stage_number", "manual_stage", "manual_coverage"]] = manual[
        ["manual_stage_number", "manual_stage", "manual_coverage"]
    ]
    if hypnodensity_5s is not None:
        predicted = aggregate_hypnodensity(hypnodensity_5s, epochs)
        columns = ["Wake", "N1", "N2", "N3", "REM", "predicted_stage_number", "predicted_stage"]
        table[columns] = predicted[columns]
    return table
