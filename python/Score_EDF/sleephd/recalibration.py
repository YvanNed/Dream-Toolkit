"""Level-A recalibration of SleepFM hypnodensities (temperature + class biases).

The recalibrated 5-s probability is softmax(log p / T + b), where p is SleepFM's
cached 5-s softmax output. Because softmax is invariant to a per-row constant,
this equals softmax(logits / T + b): no SleepFM re-run is needed. Parameters are
log T and the biases of N1, N2 and REM; Wake is the reference (b = 0) and N3 is
fixed at b = 0 because the manual hypnograms contain no N3 (its bias would not
be identifiable).

Parameters are fitted by minimising the negative log-likelihood of the manual
stage under the 30-s probability, obtained from the 5-s probabilities with the
same overlap weights as `sleephd.alignment.aggregate_hypnodensity`. Used with
leave-one-recording-out cross-fitting, each recording's recalibrated
hypnodensity comes from parameters fitted on the other recordings only.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import log_softmax, softmax

from .alignment import align_manual_hypnogram, make_nap_epochs
from .io import NapWindow

PROBABILITY_COLUMNS = ["Wake", "N1", "N2", "N3", "REM"]
FREE_BIASES = [1, 2, 4]  # N1, N2, REM; Wake and N3 fixed at 0


@dataclass
class Recording:
    """Everything the recalibration needs from one recording."""

    subject: str
    five_second: pd.DataFrame   # start_seconds, stop_seconds, Wake..REM (SleepFM, cached)
    log_p: np.ndarray           # (steps, 5) log of the 5-s probabilities
    weights: np.ndarray         # (epochs, steps) overlap weights, rows sum to 1
    labels: np.ndarray          # (epochs,) manual stage 0-4, -1 when unscored


def prepare_recording(subject: str, five_second: pd.DataFrame, window: NapWindow,
                      hypnogram: pd.DataFrame, epoch_seconds: float = 30.0) -> Recording:
    epochs = make_nap_epochs(window, epoch_seconds)
    starts = five_second["start_seconds"].to_numpy(float)
    stops = five_second["stop_seconds"].to_numpy(float)
    overlap = np.maximum(0.0, np.minimum(stops[None], epochs["stop_seconds"].to_numpy()[:, None])
                         - np.maximum(starts[None], epochs["start_seconds"].to_numpy()[:, None]))
    weights = overlap / overlap.sum(axis=1, keepdims=True)
    manual = align_manual_hypnogram(hypnogram, epochs, epoch_seconds)["manual_stage_number"]
    labels = manual.fillna(-1).to_numpy(int)
    probabilities = five_second[PROBABILITY_COLUMNS].to_numpy(float)
    log_p = np.log(np.clip(probabilities, np.finfo(float).tiny, 1.0))
    return Recording(subject, five_second, log_p, weights, labels)


def _unpack(theta: np.ndarray) -> tuple[float, np.ndarray]:
    bias = np.zeros(len(PROBABILITY_COLUMNS))
    bias[FREE_BIASES] = theta[1:]
    return float(np.exp(theta[0])), bias


def apply(log_p: np.ndarray, temperature: float, bias: np.ndarray) -> np.ndarray:
    return softmax(log_p / temperature + bias, axis=1)


def _negative_log_likelihood(theta: np.ndarray, recordings: list[Recording]) -> float:
    temperature, bias = _unpack(theta)
    total, count = 0.0, 0
    for rec in recordings:
        scored = rec.labels >= 0
        p30 = rec.weights[scored] @ apply(rec.log_p, temperature, bias)
        total -= np.log(np.clip(p30[np.arange(scored.sum()), rec.labels[scored]], 1e-12, 1)).sum()
        count += scored.sum()
    return total / count


def fit(recordings: list[Recording]) -> dict[str, float]:
    """Maximum-likelihood temperature and biases (starting from the identity)."""
    result = minimize(_negative_log_likelihood, np.zeros(1 + len(FREE_BIASES)), args=(recordings,),
                      method="L-BFGS-B", bounds=[(-3, 3)] + [(-15, 15)] * len(FREE_BIASES))
    temperature, bias = _unpack(result.x)
    return {"temperature": temperature, **{f"bias_{PROBABILITY_COLUMNS[i]}": bias[i] for i in FREE_BIASES},
            "train_nll": float(result.fun), "converged": bool(result.success), "iterations": int(result.nit)}


def recalibrated_five_second(recording: Recording, parameters: dict[str, float]) -> pd.DataFrame:
    bias = np.zeros(len(PROBABILITY_COLUMNS))
    for i in FREE_BIASES:
        bias[i] = parameters[f"bias_{PROBABILITY_COLUMNS[i]}"]
    out = recording.five_second[["start_seconds", "stop_seconds"]].copy()
    out[PROBABILITY_COLUMNS] = apply(recording.log_p, parameters["temperature"], bias)
    return out


# --- Level B: re-fit SleepFM's final linear layer (fc) on frozen head features --------
FREE_ROWS = [0, 1, 2, 4]  # Wake, N1, N2, REM; the N3 row stays at its published values
LAMBDA_GRID = (1e-4, 1e-3, 1e-2, 1e-1, 1.0)


@dataclass
class HeadRecording:
    """Frozen BiLSTM output (input of fc) aligned to the nap epochs of one recording."""

    subject: str
    five_second: pd.DataFrame   # start_seconds, stop_seconds of each 5-s step
    features: np.ndarray        # (steps, 128)
    weights: np.ndarray         # (epochs, steps) overlap weights, rows sum to 1
    labels: np.ndarray          # (epochs,) manual stage 0-4, -1 when unscored


def prepare_head_recording(base: Recording, features: np.ndarray) -> HeadRecording:
    return HeadRecording(base.subject, base.five_second[["start_seconds", "stop_seconds"]].copy(),
                         features, base.weights, base.labels)


def _stack(recordings: list[HeadRecording]):
    import torch

    features = torch.from_numpy(np.concatenate([r.features for r in recordings])).double()
    rows, cols, values, labels, offset, row = [], [], [], [], 0, 0
    for r in recordings:
        scored = np.flatnonzero(r.labels >= 0)
        block = r.weights[scored]
        i, j = np.nonzero(block)
        rows.append(i + row)
        cols.append(j + offset)
        values.append(block[i, j])
        labels.append(r.labels[scored])
        row += len(scored)
        offset += len(r.features)
    weights = torch.sparse_coo_tensor(np.vstack([np.concatenate(rows), np.concatenate(cols)]),
                                      np.concatenate(values), (row, offset)).double()
    return features, weights, torch.from_numpy(np.concatenate(labels)).long()


def _fc_nll(weight, bias, features, weights, labels):
    import torch

    p5 = torch.softmax(features @ weight.T + bias, dim=1)
    p30 = torch.sparse.mm(weights, p5)
    return -torch.log(p30[torch.arange(len(labels)), labels].clamp_min(1e-12)).mean()


def fit_fc(recordings: list[HeadRecording], weight0: np.ndarray, bias0: np.ndarray, penalty: float) -> tuple[np.ndarray, np.ndarray]:
    """Penalised maximum likelihood of fc, starting from and shrunk toward the published fc."""
    import torch

    torch.set_num_threads(1)
    features, weights, labels = _stack(recordings)
    w0, b0 = torch.from_numpy(weight0).double(), torch.from_numpy(bias0).double()
    free = torch.zeros(5, dtype=torch.bool)
    free[FREE_ROWS] = True
    dw = torch.zeros_like(w0, requires_grad=True)
    db = torch.zeros_like(b0, requires_grad=True)
    optimizer = torch.optim.LBFGS([dw, db], max_iter=500, tolerance_grad=1e-7, line_search_fn="strong_wolfe")

    def closure():
        optimizer.zero_grad()
        weight = w0 + dw * free[:, None]
        bias = b0 + db * free
        loss = _fc_nll(weight, bias, features, weights, labels) + penalty * ((dw[free] ** 2).sum() + (db[free] ** 2).sum())
        loss.backward()
        return loss

    optimizer.step(closure)
    with torch.no_grad():
        return (w0 + dw * free[:, None]).numpy(), (b0 + db * free).numpy()


def heldout_nll(recordings: list[HeadRecording], weight: np.ndarray, bias: np.ndarray) -> tuple[float, int]:
    import torch

    features, weights, labels = _stack(recordings)
    with torch.no_grad():
        nll = _fc_nll(torch.from_numpy(weight).double(), torch.from_numpy(bias).double(), features, weights, labels)
    return float(nll) * len(labels), len(labels)


def fit_fc_nested(recordings: list[HeadRecording], weight0: np.ndarray, bias0: np.ndarray,
                  grid: tuple[float, ...] = LAMBDA_GRID, inner_folds: int = 5) -> dict:
    """Choose the penalty by grouped inner cross-validation (log-loss), then refit on all."""
    groups = np.arange(len(recordings)) % inner_folds
    inner = {}
    for penalty in grid:
        total, count = 0.0, 0
        for g in range(inner_folds):
            train = [r for r, k in zip(recordings, groups) if k != g]
            test = [r for r, k in zip(recordings, groups) if k == g]
            weight, bias = fit_fc(train, weight0, bias0, penalty)
            nll, n = heldout_nll(test, weight, bias)
            total, count = total + nll, count + n
        inner[penalty] = total / count
    best = min(inner, key=inner.get)
    weight, bias = fit_fc(recordings, weight0, bias0, best)
    return {"penalty": best, "inner_nll": inner, "weight": weight, "bias": bias}


def refit_five_second(recording: HeadRecording, weight: np.ndarray, bias: np.ndarray) -> pd.DataFrame:
    out = recording.five_second.copy()
    out[PROBABILITY_COLUMNS] = softmax(recording.features @ weight.T + bias, axis=1)
    return out
