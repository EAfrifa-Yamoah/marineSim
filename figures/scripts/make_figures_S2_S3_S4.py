import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
R = os.environ.get("MS_RESULTS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")) + "/"
FIG = os.environ.get("MS_FIGURES", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")) + "/"

"""Render SI Figures S2, S3, S4 from R-generated CSVs (manuscript style: Okabe-Ito, no embedded titles, shared legends)."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

try:
    import figstyle as fs; fs.setup()
except Exception:
    fs = None
S = R + "si/"; OUT = FIG + "si/"; os.makedirs(OUT, exist_ok=True)
G = 50; CELL = 8.0

def grid_of(p):  # column-major cell index -> matrix (row = y, col = x)
    return np.asarray(p).reshape(G, G, order="F")

# ------------------------------------------------------------------ Figure S2
w = pd.read_csv(S + "figS2_worlds.csv"); des = pd.read_csv(S + "figS2_designs.csv"); rep = pd.read_csv(S + "figS2_repworld.csv")
fig, ax = plt.subplots(4, 3, figsize=(8.4, 10.6)); fig.subplots_adjust(left=0.05, right=0.90, top=0.97, bottom=0.05, hspace=0.22, wspace=0.10)
ext = [0, G * CELL, 0, G * CELL]
for r, stat in enumerate([True, False]):
    for c, rng in enumerate([30, 80, 200]):
        sub = w[(w.stationary == stat) & (w.range_km == rng)]
        im = ax[r, c].imshow(grid_of(sub.p.values), origin="lower", extent=ext, cmap="viridis", vmin=0, vmax=1)
        ax[r, c].text(0.03, 0.95, f"({'abcdef'[r*3+c]}) {'stationary' if stat else 'non stationary'}, {rng} km", transform=ax[r, c].transAxes, va="top", fontsize=8, color="white", fontweight="bold")
P = grid_of(rep.p.values)
for c, dsg in enumerate(["random", "clustered", "stratified"]):
    d = des[des.panel == f"design_{dsg}"]
    ax[2, c].imshow(P, origin="lower", extent=ext, cmap="Greys", vmin=0, vmax=1, alpha=0.6)
    ax[2, c].scatter(d.x[d.y_obs == 1], d.y[d.y_obs == 1], s=14, color="#D55E00", edgecolor="black", lw=0.3)
    ax[2, c].scatter(d.x[d.y_obs == 0], d.y[d.y_obs == 0], s=14, color="#56B4E9", edgecolor="black", lw=0.3)
    ax[2, c].text(0.03, 0.95, f"({'ghi'[c]}) {dsg}, n = 100", transform=ax[2, c].transAxes, va="top", fontsize=8, fontweight="bold")
for c, n in enumerate([30, 100, 500]):
    d = des[des.panel == f"n_{n}"]
    ax[3, c].imshow(P, origin="lower", extent=ext, cmap="Greys", vmin=0, vmax=1, alpha=0.6)
    ax[3, c].scatter(d.x[d.y_obs == 1], d.y[d.y_obs == 1], s=10, color="#D55E00", edgecolor="black", lw=0.3)
    ax[3, c].scatter(d.x[d.y_obs == 0], d.y[d.y_obs == 0], s=10, color="#56B4E9", edgecolor="black", lw=0.3)
    ax[3, c].text(0.03, 0.95, f"({'jkl'[c]}) random, n = {n}", transform=ax[3, c].transAxes, va="top", fontsize=8, fontweight="bold")
for a in ax.ravel():
    a.set_xticks([0, 200, 400]); a.set_yticks([0, 200, 400]); a.tick_params(labelsize=7)
for a in ax[:, 1:].ravel(): a.set_yticklabels([])
for a in ax[:3].ravel(): a.set_xticklabels([])
for a in ax[3]: a.set_xlabel("km", fontsize=8)
cax = fig.add_axes([0.915, 0.53, 0.018, 0.42]); cb = fig.colorbar(im, cax=cax); cb.set_label("Occurrence probability", fontsize=8); cb.ax.tick_params(labelsize=7)
fig.legend(handles=[Patch(fc="#D55E00", ec="black", label="presence"), Patch(fc="#56B4E9", ec="black", label="absence")],
           loc="center right", bbox_to_anchor=(0.995, 0.28), frameon=False, fontsize=8, title="Sampled cells", title_fontsize=8)
fig.savefig(OUT + "FigureS2.png", dpi=300); plt.close(fig); print("S2")

# ------------------------------------------------------------------ Figure S3
t = pd.read_csv(S + "figS3_tuning.csv")
agg = t.groupby(["method", "min_node_size", "mtry_rule"]).auc.mean().reset_index()
order = [(m, r) for m in [1, 3, 5] for r in ["sqrt", "half", "p8"]]
lab = {"sqrt": "√p", "half": "0.5p", "p8": "0.8p"}
MC = {"Standard RF": "#E69F00", "Spatial RF": "#56B4E9", "stRF": "#D55E00"}
fig, ax = plt.subplots(figsize=(7.2, 4.0)); fig.subplots_adjust(left=0.10, right=0.98, top=0.95, bottom=0.30)
x = np.arange(len(order))
for m in MC:
    y = [agg[(agg.method == m) & (agg.min_node_size == a) & (agg.mtry_rule == b)].auc.iloc[0] for a, b in order]
    ax.plot(x, y, "-o", color=MC[m], ms=5, lw=1.6, label=m)
ax.axvline(0, color="0.6", lw=0.8, ls="--"); ax.text(0.1, ax.get_ylim()[0] + 0.003, "default", fontsize=7, color="0.4")
ax.set_xticks(x); ax.set_xticklabels([f"node {a}\nmtry {lab[b]}" for a, b in order], fontsize=7.5)
ax.set_ylabel("Mean test AUC"); ax.set_xlabel("Hyperparameter setting (minimum node size, features per split)", fontsize=8.5)
fig.legend(handles=[Line2D([], [], color=MC[m], marker="o", ms=5, lw=1.6, label=m) for m in MC], loc="lower center", ncol=3, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, 0.01))
fig.savefig(OUT + "FigureS3.png", dpi=300); plt.close(fig); print("S3")

# ------------------------------------------------------------------ Figure S4
c = pd.read_csv(S + "figS4_calibration.csv")
MC5 = {"stRF": "#D55E00", "Spatial RF": "#56B4E9", "Standard RF": "#E69F00", "GLM": "#999999", "GAM": "#009E73"}
bins = np.linspace(0, 1, 11)
fig, ax = plt.subplots(1, 2, figsize=(8.6, 4.0), gridspec_kw=dict(width_ratios=[1.15, 1])); fig.subplots_adjust(left=0.08, right=0.98, top=0.94, bottom=0.26, wspace=0.28)
ax[0].plot([0, 1], [0, 1], "--", color="0.6", lw=1)
for m in MC5:
    s = c[c.method == m]; b = np.clip(np.digitize(s.p_hat, bins) - 1, 0, 9)
    xs, ys = [], []
    for k in range(10):
        sel = b == k
        if sel.sum() >= 30: xs.append(s.p_hat[sel].mean()); ys.append(s.y[sel].mean())
    ax[0].plot(xs, ys, "-o", color=MC5[m], ms=4.5, lw=1.5)
ax[0].set_xlabel("Predicted probability"); ax[0].set_ylabel("Observed frequency"); ax[0].set_xlim(0, 1); ax[0].set_ylim(0, 1)
ax[0].text(0.02, 0.97, "(a)", transform=ax[0].transAxes, va="top", fontweight="bold")
# (b) Brier by method with SD
tab5 = pd.read_csv(S + "tableS5_brier.csv")
ax[1].barh(range(len(tab5)), tab5.brier, xerr=tab5.sd, color=[MC5[m] for m in tab5.method], edgecolor="black", lw=0.5, capsize=3)
ax[1].set_yticks(range(len(tab5))); ax[1].set_yticklabels(tab5.method, fontsize=8.5); ax[1].invert_yaxis()
ax[1].set_xlabel("Brier score (lower is better)"); ax[1].text(0.02, 0.97, "(b)", transform=ax[1].transAxes, va="top", fontweight="bold")
fig.legend(handles=[Line2D([], [], color=MC5[m], marker="o", ms=4.5, lw=1.5, label=m) for m in MC5], loc="lower center", ncol=5, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, 0.01))
fig.savefig(OUT + "FigureS4.png", dpi=300); plt.close(fig); print("S4")
