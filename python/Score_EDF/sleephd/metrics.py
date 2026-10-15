from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, cohen_kappa_score, confusion_matrix, f1_score

STAGES = ["Wake", "N1", "N2", "N3", "REM"]


def staging_metrics(epoch_table: pd.DataFrame) -> tuple[pd.Series, pd.DataFrame]:
    valid = epoch_table["manual_stage_number"].notna() & epoch_table["predicted_stage_number"].notna()
    truth = epoch_table.loc[valid, "manual_stage_number"].astype(int)
    predicted = epoch_table.loc[valid, "predicted_stage_number"].astype(int)
    metrics = pd.Series(
        {
            "n_scored_epochs": len(truth),
            "accuracy": accuracy_score(truth, predicted),
            "macro_f1": f1_score(truth, predicted, labels=range(5), average="macro", zero_division=0),
            "cohen_kappa": cohen_kappa_score(truth, predicted, labels=range(5)),
        },
        name="value",
    )
    matrix = confusion_matrix(truth, predicted, labels=range(5), normalize="true")
    matrix = np.nan_to_num(matrix)
    return metrics, pd.DataFrame(matrix, index=STAGES, columns=STAGES)
