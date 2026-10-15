"""Individual (within-recording) hypnodensities: catch22 + XGBoost, leave-one-epoch-out.

The feature extraction and classifier are vendored from
https://github.com/rubenherzog/ieeg-sleep-staging (commit 0ac620e):
`staging/catch22_features.py::get_catch22_from_epochs` and
`staging/single_night.py::{loo_train_and_predict, loo_xgboost_multiclass_par}`,
with their logic unchanged (only the timing print and the unused metrics are
dropped). Preprocessing follows that repo's configuration: 128-Hz sampling and a
0.1-45 Hz MNE band-pass on every EEG/EOG/EMG channel, 30-s epochs, 22 catch22
features per channel.

Differences from the repo, both deliberate:
- The Platt calibration in `hypnos_from_epochs` is not applied: it is fitted on
  the same epochs (and labels) it calibrates, so its output would no longer be
  out-of-sample.
- Epochs are the LIGHTS OFF-anchored nap grid. Epochs with a manual label are
  predicted by leave-one-epoch-out; unscored epochs and the final partial epoch
  are predicted by a model trained on all labeled epochs (their features use
  the 30-s window starting at the epoch onset).

The model is fitted to each recording's own manual hypnogram, so it cannot
stage recordings without one. Adjacent epochs are autocorrelated, so
leave-one-epoch-out agreement is optimistic; `buffer_epochs` also removes the
neighbours of the held-out epoch from training to quantify this.
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from .alignment import align_manual_hypnogram, make_nap_epochs
from .io import NapWindow
from .signals import EEG_DERIVATIONS, _load_edf_channels

PROBABILITY_COLUMNS = ["Wake", "N1", "N2", "N3", "REM"]
INDIVIDUAL_SFREQ = 128.0          # ieeg-sleep-staging: 4096 Hz / 32
INDIVIDUAL_BAND = (0.1, 45.0)     # ieeg-sleep-staging: EEG_F_LO, EEG_F_HI
INDIVIDUAL_CHANNELS = (*EEG_DERIVATIONS, "vEOG", "hEOG", "EMG")


# --- vendored from ieeg-sleep-staging --------------------------------------------------
def get_catch22_from_epochs(data, ncpus=4):
    """catch22 for each channel of epoched data (nepochs x nchan x nsamps).

    Returns a (nchan * 22) x nepochs feature matrix.
    """
    import pycatch22
    from joblib import Parallel, delayed

    nfeats = 22
    nepochs, nchan, nsamps = np.shape(data)
    feat_mat = np.zeros((nchan, nfeats, nepochs))
    for c in range(nchan):
        aux_out = Parallel(n_jobs=ncpus)(
            delayed(pycatch22.catch22_all)(data[e, c, :]) for e in range(nepochs))
        feat_mat[c, :, :] = np.array([np.array(a["values"]) for a in aux_out]).T
    return np.reshape(feat_mat, [nchan * nfeats, nepochs])


def loo_train_and_predict(x_train, y_train, x_test, y_test, params):
    """Train an XGBoost model and predict (worker of loo_xgboost_multiclass_par)."""
    import xgboost as xgb

    dtrain = xgb.DMatrix(x_train, label=y_train)
    dtest = xgb.DMatrix(x_test, label=y_test)
    try:
        classifier = xgb.train(params, dtrain)
        y_pred = classifier.predict(dtest)
    except Exception:
        y_pred = np.zeros((len(y_test), len(np.unique(y_train))))
    return y_pred


def xgboost_params(n_classes: int) -> dict:
    return {
        'objective': 'multi:softprob',
        'num_class': n_classes,
        'eval_metric': 'mlogloss',
        'learning_rate': 0.1,
        'max_depth': 6,
        'subsample': 0.8,
        'colsample_bytree': 0.8,
        'gamma': 0.1,
        'min_child_weight': 1,
        'n_jobs': -1
    }


def loo_xgboost_multiclass_par(x, y, ncpus=4):
    """Leave-one-out XGBoost probabilities, parallelized across folds."""
    from joblib import Parallel, delayed
    from sklearn.model_selection import LeaveOneOut

    params = xgboost_params(len(np.unique(y)))
    return Parallel(n_jobs=ncpus)(delayed(loo_train_and_predict)(
        x[train], y[train], x[test], y[test], params
    ) for train, test in LeaveOneOut().split(x, y))
# --------------------------------------------------------------------------------------


def _epoch_data(edf_path: str | Path, epochs: pd.DataFrame, channels: Sequence[str]) -> np.ndarray:
    """(epochs, channels, 30 s * 128 Hz) array, resampled then band-passed as in the repo."""
    import mne

    signals, sfreq = _load_edf_channels(edf_path)
    data = []
    for name in channels:
        if name in EEG_DERIVATIONS:
            active, reference = EEG_DERIVATIONS[name]
            data.append(signals[active] - signals[reference])
        else:
            data.append(signals[name])
    raw = mne.io.RawArray(np.stack(data), mne.create_info(list(channels), sfreq, "misc"), verbose="error")
    raw.resample(INDIVIDUAL_SFREQ, verbose="error")
    raw.filter(*INDIVIDUAL_BAND, picks="all", verbose="error")
    values = raw.get_data()
    n = int(30 * INDIVIDUAL_SFREQ)
    starts = np.round(epochs["start_seconds"].to_numpy() * INDIVIDUAL_SFREQ).astype(int)
    starts = np.minimum(starts, values.shape[1] - n)  # final epoch: last 30 s of data if needed
    return np.stack([values[:, s:s + n] for s in starts])


def infer_individual(
    edf_path: str | Path,
    window: NapWindow,
    hypnogram: pd.DataFrame,
    channels: Sequence[str] = INDIVIDUAL_CHANNELS,
    epoch_seconds: float = 30.0,
    buffer_epochs: int = 0,
    cache_path: str | Path | None = None,
    force: bool = False,
    ncpus: int = 1,
) -> pd.DataFrame:
    """Out-of-sample hypnodensity for every nap epoch, fitted to this recording's labels.

    Rows are the nap epochs (`start_seconds`, `stop_seconds`, Wake..REM) plus
    `held_out` (True where the epoch had a manual label and was predicted by
    leave-one-epoch-out). With `buffer_epochs = k`, the k epochs on each side
    of the held-out epoch are also excluded from its training set.
    """
    if cache_path is not None and Path(cache_path).exists() and not force:
        return pd.read_csv(cache_path)

    epochs = make_nap_epochs(window, epoch_seconds)
    manual = align_manual_hypnogram(hypnogram, epochs, epoch_seconds)
    labeled = manual["manual_stage_number"].notna().to_numpy() & np.isclose(epochs["duration_seconds"], epoch_seconds)
    features = get_catch22_from_epochs(_epoch_data(edf_path, epochs, channels), ncpus=ncpus).T
    stage = manual["manual_stage_number"].to_numpy()

    probabilities = np.full((len(epochs), len(PROBABILITY_COLUMNS)), np.nan)
    classes = np.unique(stage[labeled]).astype(int)
    if len(classes) == 1:  # a single manual stage: the classifier is trivial
        probabilities[:] = 0.0
        probabilities[:, classes[0]] = 1.0
    else:
        remap = {c: i for i, c in enumerate(classes)}
        index = np.flatnonzero(labeled)
        x, y = features[index], np.array([remap[int(c)] for c in stage[index]])
        if buffer_epochs == 0:
            held_out = loo_xgboost_multiclass_par(x, y, ncpus=ncpus)
        else:
            from joblib import Parallel, delayed
            params = xgboost_params(len(classes))
            splits = [(np.abs(index - index[i]) > buffer_epochs, i) for i in range(len(index))]
            held_out = Parallel(n_jobs=ncpus)(delayed(loo_train_and_predict)(
                x[keep], y[keep], x[[i]], y[[i]], params) for keep, i in splits)
        for i, prediction in zip(index, held_out):
            row = np.zeros(len(PROBABILITY_COLUMNS))
            prediction = np.asarray(prediction).reshape(-1)
            row[classes[:len(prediction)]] = prediction
            probabilities[i] = row
        others = np.flatnonzero(~labeled)
        if len(others):
            prediction = loo_train_and_predict(x, y, features[others], np.zeros(len(others)),
                                               xgboost_params(len(classes)))
            probabilities[others] = 0.0
            probabilities[np.ix_(others, classes)] = prediction
    # A training fold that lost a class yields fewer columns; renormalise defensively.
    probabilities /= probabilities.sum(axis=1, keepdims=True)

    result = pd.DataFrame(probabilities, columns=PROBABILITY_COLUMNS)
    result.insert(0, "stop_seconds", epochs["stop_seconds"].to_numpy())
    result.insert(0, "start_seconds", epochs["start_seconds"].to_numpy())
    result["held_out"] = labeled
    if cache_path is not None:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(cache_path, index=False)
    return result
