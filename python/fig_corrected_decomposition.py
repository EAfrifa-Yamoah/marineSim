"""
fig_corrected_decomposition.py
Drop-in replacement for Figure 2.

(a) Shapley decomposition of the AUC gain over the aspatial GLM into the four
    feature blocks plus algorithmic uplift, order free.
(b) Placebo test. At one timestep the temporal features are constant by
    construction and must contribute nothing. The published contrast
    (stRF - Spatial RF) reports +0.047 AUC there; the corrected increment
    reports approximately zero.
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = "../results"
OK = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
      "vermillion": "#D55E00", "purple": "#CC79A7", "grey": "#999999"}

d = json.load(open(f"{R}/decomposition_corrected.json"))
df = pd.read_csv(f"{R}/ablation_decomposed.csv")

fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.0))

# ---- (a) Shapley decomposition -------------------------------------------
labels = ["Algorithmic\nflexibility", "Coordinates", "Distance\nfields",
          "Spatial\nlags", "Temporal\nfeatures"]
keys = ["ML", "COORD", "EDF", "LAG", "TEMP"]
cols = [OK["grey"], OK["blue"], OK["blue"], OK["blue"], OK["orange"]]
est = [d["shapley"][k][0] for k in keys]
lo = [d["shapley"][k][0] - d["shapley"][k][1] for k in keys]
hi = [d["shapley"][k][2] - d["shapley"][k][0] for k in keys]
total = d["shapley_total"]

b = ax[0].bar(range(5), est, yerr=[lo, hi], capsize=3, color=cols,
              edgecolor="black", linewidth=0.6, width=0.66)
for i, v in enumerate(est):
    ax[0].text(i, v + hi[i] + 0.0035, f"{100*v/total:.0f}%",
               ha="center", va="bottom", fontsize=9)
ax[0].set_xticks(range(5))
ax[0].set_xticklabels(labels, fontsize=8.5)
ax[0].set_ylabel("Shapley contribution to test AUC")
ax[0].set_ylim(0, max(est) + 0.022)
ax[0].axhline(0, color="black", lw=0.8)
ax[0].text(0.02, 0.95, "(a)", transform=ax[0].transAxes, fontweight="bold",
           va="top")
ax[0].grid(axis="y", alpha=0.25, lw=0.5)
ax[0].set_axisbelow(True)

# ---- (b) placebo test by temporal depth ----------------------------------
g = df.groupby("n_time")[["published_d_temporal", "temporal_last"]].agg(
    ["mean", "sem"])
T = g.index.values
w = 0.34
pub_m = g[("published_d_temporal", "mean")].values
pub_e = 1.96 * g[("published_d_temporal", "sem")].values
cor_m = g[("temporal_last", "mean")].values
cor_e = 1.96 * g[("temporal_last", "sem")].values

x = np.arange(len(T))
ax[1].bar(x - w/2, pub_m, w, yerr=pub_e, capsize=3, color=OK["vermillion"],
          edgecolor="black", linewidth=0.6,
          label="Published contrast (stRF − Spatial RF)")
ax[1].bar(x + w/2, cor_m, w, yerr=cor_e, capsize=3, color=OK["orange"],
          edgecolor="black", linewidth=0.6,
          label="Corrected increment (temporal block only)")
ax[1].axhline(0, color="black", lw=0.8)
ax[1].set_xticks(x)
ax[1].set_xticklabels([f"{t}" for t in T])
ax[1].set_xlabel("Number of timesteps")
ax[1].set_ylabel("Attributed temporal gain in test AUC")
ax[1].text(0.02, 0.95, "(b)", transform=ax[1].transAxes, fontweight="bold",
           va="top")
ax[1].grid(axis="y", alpha=0.25, lw=0.5)
ax[1].set_axisbelow(True)
ax[1].annotate("temporal features carry\nno information at T = 1",
               xy=(0 + w/2, cor_m[0]), xytext=(0.42, 0.62),
               textcoords="axes fraction", fontsize=8.5,
               arrowprops=dict(arrowstyle="->", lw=0.8, color="black"))

handles, lab = ax[1].get_legend_handles_labels()
fig.legend(handles, lab, loc="lower center", ncol=2, frameon=False,
           fontsize=9, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=[0, 0.06, 1, 1])
fig.savefig(f"{R}/Figure2_corrected_decomposition.png", dpi=300,
            bbox_inches="tight")
print("wrote Figure2_corrected_decomposition.png")
