"""Figure 3: the five method core benchmark from results/core_five_methods.csv."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
R = os.environ.get("MS_RESULTS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")) + "/"
FIG = os.environ.get("MS_FIGURES", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")) + "/"

import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})

c = pd.read_csv(R + "core_five_methods.csv")
c["field"] = c.world_cfg + "_" + c.realisation.astype(str)
M = [("GLM", "GLM", "#999999"), ("GAM", "GAM", "#009E73"), ("RF", "Standard RF", "#E69F00"),
     ("SRF", "Spatial RF", "#56B4E9"), ("STRF", "stRF", "#D55E00")]
ns = sorted(c.n.unique())

fig, ax = plt.subplots(1, 2, figsize=(10.4, 4.2), gridspec_kw=dict(width_ratios=[1.15, 1]))
fig.subplots_adjust(left=0.07, right=0.985, top=0.95, bottom=0.21, wspace=0.24)

# (a) pooled AUC vs n, bands = 2.5-97.5 percentiles of field-realisation means
a = ax[0]
for col, lab, colr in M:
    fm = c.groupby(["field", "n"])[col].mean().reset_index()
    mean = fm.groupby("n")[col].mean().reindex(ns)
    lo = fm.groupby("n")[col].quantile(0.025).reindex(ns)
    hi = fm.groupby("n")[col].quantile(0.975).reindex(ns)
    a.fill_between(ns, lo, hi, color=colr, alpha=0.15, lw=0)
    a.plot(ns, mean, "-o", color=colr, ms=4.5, lw=1.8)
a.axhline(0.5, ls=":", color="k", lw=1)
a.set_xscale("log"); a.set_xticks(ns); a.set_xticklabels([str(n) for n in ns]); a.minorticks_off()
a.set_xlabel("Sites per timestep"); a.set_ylabel("Mean test AUC")
a.text(-0.10, 1.02, "(a)", transform=a.transAxes, fontweight="bold")

# (b) change relative to random sampling, by method
b = ax[1]
dm = c.groupby(["design", "n"])[[m[0] for m in M]].mean()
for col, lab, colr in M:
    for des, ls, mk in [("clustered", "-", "o"), ("stratified", "--", "s")]:
        d = (dm.loc[des, col] - dm.loc["random", col]).reindex(ns)
        b.plot(ns, d, ls=ls, marker=mk, color=colr, ms=4.5, lw=1.6)
b.axhline(0, color="k", lw=1)
b.set_xscale("log"); b.set_xticks(ns); b.set_xticklabels([str(n) for n in ns]); b.minorticks_off()
b.set_xlabel("Sites per timestep"); b.set_ylabel("Change in mean AUC relative to random")
b.text(-0.13, 1.02, "(b)", transform=b.transAxes, fontweight="bold")

h1 = [Line2D([], [], color=colr, lw=2, marker="o", ms=4.5, label=lab) for _, lab, colr in M]
h2 = [Line2D([], [], color="k", lw=1.6, ls="-", marker="o", ms=4.5, label="Clustered vs random"),
      Line2D([], [], color="k", lw=1.6, ls="--", marker="s", ms=4.5, label="Stratified vs random")]
fig.legend(handles=h1 + h2, loc="lower center", ncol=7, frameon=False, bbox_to_anchor=(0.5, 0.015),
           columnspacing=1.4, handlelength=2.4)
os.makedirs(FIG + "main", exist_ok=True); fig.savefig(FIG + "main/Figure3.png", dpi=300)

# numbers for the caption/text check
print(c.groupby("n")[[m[0] for m in M]].mean().round(3))
print((dm.loc["clustered"] - dm.loc["random"]).round(3))
print((dm.loc["stratified"] - dm.loc["random"]).round(3))
