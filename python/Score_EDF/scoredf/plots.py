"""Matplotlib figures for stager comparison results."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import seaborn as sns
from sklearn.metrics import cohen_kappa_score, f1_score


STAGES = ["Wake", "N1", "N2", "N3", "REM"]

STAGE_COLORS = {
    "Wake": "#7f7f7f",
    "N1": "#f4a261",
    "N2": "#2a9d8f",
    "N3": "#264653",
    "REM": "#9b5de5",
}

STAGER_COLORS = {
    "SleepFM": "#2a78d6",
    "SleepFM recal A": "#008300",
    "SleepFM recal B": "#7b3fb0",
    "GSSC": "#eb6834",
    "YASA": "#1baf7a",
}


def plot_hypnodensities(
    epochs_df: pd.DataFrame,
    stagers: list[str],
    subject: str,
) -> plt.Figure:
    """Stacked-area hypnodensity per stager + manual hypnogram for one recording."""
    n_panels = len(stagers) + 1  # +1 for manual hypnogram
    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(13, 1.4 * n_panels),
        sharex=True,
        gridspec_kw={"height_ratios": [1] * len(stagers) + [0.9]},
    )
    if n_panels == 1:
        axes = [axes]

    for ax, name in zip(axes, stagers):
        g = epochs_df.loc[
            (epochs_df["subject"] == subject) & (epochs_df["stager"] == name)
        ]
        if g.empty:
            ax.set_ylabel(name, rotation=0, ha="right", va="center", fontsize=9)
            ax.text(0.5, 0.5, "not available", transform=ax.transAxes, ha="center", va="center", color="grey")
            continue
        ax.stackplot(
            g["minutes_from_start"],
            *[g[s] for s in STAGES],
            colors=list(STAGE_COLORS.values()),
            step="post",
            alpha=0.9,
        )
        ax.set(ylim=(0, 1), yticks=[0, 1], xlim=(0, g["minutes_from_start"].max()))
        ax.set_ylabel(name, rotation=0, ha="right", va="center", fontsize=9)

    # Manual hypnogram panel
    g = epochs_df.loc[
        (epochs_df["subject"] == subject) & (epochs_df["stager"] == stagers[0])
    ]
    if not g.empty:
        ax_manual = axes[-1]
        ax_manual.step(
            g["minutes_from_start"],
            g["manual_stage_number"],
            where="post",
            color="black",
            lw=1.5,
        )
        present_stages = sorted(g["manual_stage_number"].dropna().unique())
        yticks = [int(s) for s in present_stages]
        ylabels = [STAGES[int(s)] for s in present_stages]
        ax_manual.set(
            yticks=yticks,
            yticklabels=ylabels,
            ylim=(min(yticks) - 0.3, max(yticks) + 0.3) if yticks else (0, 4),
        )
        ax_manual.set_ylabel("Manual", rotation=0, ha="right", va="center", fontsize=9)
        ax_manual.set_xlabel("Minutes from analysis start")

    axes[0].legend(
        handles=[Patch(color=c, label=s) for s, c in STAGE_COLORS.items()],
        ncol=5,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        frameon=False,
        fontsize=8,
    )
    fig.suptitle(f"{subject}: hypnodensity comparison", y=1.0, fontsize=11)
    fig.tight_layout()
    return fig


def plot_confusion_matrices(
    scored_df: pd.DataFrame,
    stagers: list[str],
) -> plt.Figure:
    """Confusion matrix heatmaps (pooled), one per stager."""
    n = len(stagers)
    cols = min(n, 4)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3.5 * rows))
    axes = np.atleast_2d(axes)

    labels = [s for s in STAGES if (scored_df["manual_stage"] == s).any()]

    for idx, name in enumerate(stagers):
        ax = axes.flat[idx]
        g = scored_df.loc[scored_df["stager"] == name]
        if g.empty:
            ax.set_title(name)
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", color="grey")
            continue
        counts = pd.crosstab(
            g["manual_stage"], g["predicted_stage"]
        ).reindex(index=labels, columns=STAGES, fill_value=0)
        row_sums = counts.sum(axis=1)
        pct = counts.div(row_sums.replace(0, 1), axis=0) * 100
        sns.heatmap(
            pct,
            annot=True,
            fmt=".1f",
            cmap="Blues",
            cbar=False,
            vmin=0,
            vmax=100,
            ax=ax,
        )
        ax.set(title=name, xlabel="Predicted", ylabel="Manual")

    for idx in range(len(stagers), axes.size):
        axes.flat[idx].axis("off")

    fig.suptitle("Confusion matrices % (pooled recordings)", y=1.01, fontsize=11)
    fig.tight_layout()
    return fig


def compute_agreement(
    scored_df: pd.DataFrame,
    stagers: list[str],
    subjects: list[str],
) -> pd.DataFrame:
    """Compute agreement metrics per stager (pooled and per recording)."""
    labels = [s for s in STAGES if (scored_df["manual_stage"] == s).any()]
    rows = []
    for name in stagers:
        g = scored_df.loc[scored_df["stager"] == name]
        if g.empty:
            continue
        row: dict = {
            "stager": name,
            "epochs": len(g),
            "accuracy": (g["manual_stage"] == g["predicted_stage"]).mean(),
            "kappa_pooled": cohen_kappa_score(g["manual_stage"], g["predicted_stage"]),
        }
        f1_vals = f1_score(
            g["manual_stage"],
            g["predicted_stage"],
            labels=labels,
            average=None,
            zero_division=0,
        )
        for stage, val in zip(labels, f1_vals):
            row[f"F1_{stage}"] = val
        for sid in subjects:
            sub = g.loc[g["subject"] == sid]
            if len(sub) >= 2:
                row[f"kappa_{sid}"] = cohen_kappa_score(sub["manual_stage"], sub["predicted_stage"])
        rows.append(row)
    return pd.DataFrame(rows).set_index("stager")


def plot_agreement_summary(
    agreement_df: pd.DataFrame,
    stagers: list[str],
    subjects: list[str],
) -> plt.Figure:
    """F1 per stage (left) and kappa per recording (right)."""
    labels = [s for s in STAGES if f"F1_{s}" in agreement_df.columns]
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.5), gridspec_kw={"width_ratios": [1.3, 1]})

    # F1 bars
    present_stagers = [s for s in stagers if s in agreement_df.index]
    f1_cols = [f"F1_{s}" for s in labels]
    f1 = agreement_df.loc[present_stagers, f1_cols].rename(columns=lambda c: c.removeprefix("F1_"))
    colors = [STAGER_COLORS.get(s, "#888888") for s in present_stagers]
    f1.T.plot.bar(ax=axes[0], color=colors, width=0.8, rot=0)
    axes[0].set(ylim=(0, 1), ylabel="F1 (pooled)", title="Agreement per stage")
    axes[0].legend(title="", ncol=min(4, len(present_stagers)), fontsize=7,
                   loc="lower center", bbox_to_anchor=(0.5, 1.05), frameon=False)

    # Kappa per recording
    markers = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "+"]
    kappa_cols = [c for c in agreement_df.columns if c.startswith("kappa_") and c != "kappa_pooled"]
    for i, col in enumerate(kappa_cols):
        sid = col.removeprefix("kappa_")
        marker = markers[i % len(markers)]
        vals = agreement_df.loc[present_stagers, col].values
        axes[1].plot(range(len(present_stagers)), vals, color="#8a8a8a", lw=0.8,
                     marker=marker, ms=7, label=sid)
    axes[1].scatter(
        range(len(present_stagers)),
        agreement_df.loc[present_stagers, "kappa_pooled"],
        marker="_", s=700, color="black", zorder=3, label="pooled",
    )
    axes[1].set(
        xticks=range(len(present_stagers)),
        ylim=(0, 1),
        ylabel="Cohen's kappa",
        title="Agreement per recording",
    )
    axes[1].set_xticklabels(present_stagers, rotation=30, ha="right", fontsize=8)
    axes[1].legend(title="", fontsize=7)

    fig.tight_layout()
    return fig


def plot_hypnogram_grid(
    epochs_df: pd.DataFrame,
    stagers: list[str],
    subjects: list[str],
) -> plt.Figure:
    """Color-coded epoch grid: manual + each stager, one subplot per recording."""
    stage_code = {s: i for i, s in enumerate(STAGES)}
    cmap = ListedColormap(list(STAGE_COLORS.values()))
    cmap.set_bad("#e0e0e0")

    fig, axes = plt.subplots(len(subjects), 1, figsize=(13, 2.6 * len(subjects)))
    axes_arr = np.atleast_1d(axes)

    for ax, sid in zip(axes_arr, subjects):
        g = epochs_df.loc[epochs_df["subject"] == sid]
        first_stager = next((s for s in stagers if s in g["stager"].unique()), None)
        if first_stager is None:
            ax.text(0.5, 0.5, f"{sid}: no stager data", transform=ax.transAxes, ha="center")
            continue

        manual = (
            g.loc[g["stager"] == first_stager]
            .set_index("epoch")["manual_stage"]
            .map(stage_code)
        )
        present_stagers = [s for s in stagers if s in g["stager"].unique()]
        predicted = (
            g.assign(code=g["predicted_stage"].map(stage_code))
            .pivot(index="stager", columns="epoch", values="code")
            .reindex(present_stagers)
        )
        image = np.vstack([manual.to_numpy(float), predicted.to_numpy(float)])
        ax.imshow(
            np.ma.masked_invalid(image),
            aspect="auto",
            cmap=cmap,
            vmin=-0.5,
            vmax=4.5,
            interpolation="nearest",
            extent=(0, image.shape[1] / 2, image.shape[0], 0),
        )
        ax.set(
            yticks=np.arange(image.shape[0]) + 0.5,
            yticklabels=["Manual", *present_stagers],
            title=sid,
        )
        ax.axhline(1, color="white", lw=3)
        ax.grid(False)

    axes_arr[-1].set_xlabel("Minutes from analysis start")
    axes_arr[0].legend(
        handles=[Patch(color=c, label=s) for s, c in STAGE_COLORS.items()],
        ncol=5,
        loc="lower right",
        bbox_to_anchor=(1, 1.02),
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout()
    return fig
