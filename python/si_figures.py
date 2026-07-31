"""Supplementary figures S1 (calibration reliability) and S2 (tuning sensitivity)."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figstyle as fs

R = os.path.join(os.path.dirname(__file__), "..", "results")
F = os.path.join(os.path.dirname(__file__), "..", "figures")
fs.setup()
MC = {"GLM": fs.OI["grey"], "GAM": fs.OI["orange"], "RF": fs.OI["skyblue"],
      "SpatialRF": fs.OI["blue"], "stRF": fs.OI["vermillion"]}
ORDER = ["GLM", "GAM", "RF", "SpatialRF", "stRF"]


def figS1_calibration():
    rel = pd.read_csv(f"{R}/si_reliability.csv")
    cal = pd.read_csv(f"{R}/si_calibration.csv").set_index("method")
    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    ax.plot([0, 1], [0, 1], ls="--", color="0.5", lw=1, zorder=1)
    for m in ORDER:
        d = rel[rel.method == m].sort_values("pred_mean")
        if d.empty:
            continue
        ax.plot(d["pred_mean"], d["obs_freq"], "-o", color=MC[m], ms=4, lw=1.6, zorder=3,
                label=f"{m}  (Brier {cal.loc[m, 'brier_mean']:.3f})")
    ax.set_xlabel("Mean predicted probability")
    ax.set_ylabel("Observed frequency of occurrence")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=8.2, title="Method (lower Brier = better)",
              title_fontsize=8.4)
    ax.text(0.97, 0.06, "Perfect calibration:\npoints on the diagonal", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7.6, color="0.4")
    fs.save(fig, f"{F}/FigureS2_calibration.png")


def figS2_tuning():
    d = pd.read_csv(f"{R}/si_tuning.csv")
    methods = ["RF", "SpatialRF", "stRF"]
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    rng = np.random.default_rng(0)
    for i, m in enumerate(methods):
        cm = d[d.method == m]
        xj = i + (rng.random(len(cm)) - 0.5) * 0.28
        ax.scatter(xj, cm["auc"], s=34, color=MC[m], edgecolor="black", linewidth=0.4,
                   alpha=0.9, zorder=3)
        dfl = cm[(cm.leaf == 1) & (cm.mf == "sqrt")]["auc"]
        if len(dfl):
            ax.scatter([i], [dfl.iloc[0]], s=120, marker="D", facecolor="none",
                       edgecolor="black", linewidth=1.3, zorder=4)
        ax.plot([i - 0.22, i + 0.22], [cm["auc"].mean()] * 2, color="black", lw=1.4, zorder=4)
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(["Standard RF", "Spatial RF", "stRF"])
    ax.set_ylabel("Test AUC across 9 hyperparameter settings")
    ax.set_xlabel("Method")
    ax.set_xlim(-0.5, len(methods) - 0.5)
    handles = [Line2D([], [], marker="o", ls="", markerfacecolor="0.6", markeredgecolor="black",
                      markersize=7, label="One hyperparameter setting"),
               Line2D([], [], marker="D", ls="", markerfacecolor="none", markeredgecolor="black",
                      markersize=9, label="Default setting"),
               Line2D([], [], color="black", lw=1.4, label="Mean across settings")]
    ax.legend(handles=handles, loc="lower right", fontsize=8)
    fs.save(fig, f"{F}/FigureS1_tuning.png")


if __name__ == "__main__":
    figS1_calibration()
    figS2_tuning()
    print("SI FIGURES DONE")
