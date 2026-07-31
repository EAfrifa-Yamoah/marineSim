"""Temporal figures: (T1) spatiotemporal reconstruction, (T2) honest forecast benchmark."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import figstyle as fs
fs.setup()

R, F = "../results", "../figures"
# colour-blind-safe qualitative palette for the five parks (distinct from method colours)
REG_ORDER = ["JBMP", "MMP", "CSMC", "SIMP", "NCMP"]
REG_LAB = {"JBMP": "Jurien Bay", "MMP": "Marmion", "CSMC": "Cockburn Sound",
           "SIMP": "Shoalwater", "NCMP": "Ngari Capes"}
REG_COL = {"JBMP": "#0072B2", "MMP": "#E69F00", "CSMC": "#009E73",
           "SIMP": "#CC79A7", "NCMP": "#5D3A9B"}


def figT1():
    traj = pd.read_csv(f"{R}/temporal_regional_traj.csv")
    fld = pd.read_csv(f"{R}/temporal_field.csv")
    fig = plt.figure(figsize=(9.2, 8.2))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.25], hspace=0.32, wspace=0.14,
                  left=0.085, right=0.965, top=0.95, bottom=0.09)

    # (a) regional trajectories: observed vs reconstruction
    axa = fig.add_subplot(gs[0, :])
    for rg in REG_ORDER:
        d = traj[traj["region"] == rg].sort_values("year")
        if d.empty:
            continue
        axa.plot(d["year"], d["obs"], "o", color=REG_COL[rg], ms=4, alpha=0.85)
        axa.plot(d["year"], d["pred"], "-", color=REG_COL[rg], lw=2, alpha=0.9)
    axa.set_xlabel("Year"); axa.set_ylabel("Mean shoot density 0.04 m$^{-2}$")
    axa.set_title("(a)  Regional trajectories: observed (points) vs stRF reconstruction "
                  "(lines, held-out sites)", fontsize=9.5, loc="left")
    handles = [Line2D([], [], color=REG_COL[r], lw=2, marker="o", ms=4, label=REG_LAB[r])
               for r in REG_ORDER]
    axa.legend(handles=handles, fontsize=7.6, ncol=5, loc="upper center",
               bbox_to_anchor=(0.5, -0.16), frameon=False)

    # (b)/(c) site x year heatmaps, sites ordered by latitude within park (north at top)
    order = (fld[["site", "region", "latitude"]].drop_duplicates("site")
             .assign(rank=lambda x: x["region"].map({r: i for i, r in enumerate(REG_ORDER)}))
             .sort_values(["rank", "latitude"], ascending=[True, False]))
    sites = order["site"].tolist()
    years = sorted(fld["year"].unique())
    obs = fld.pivot(index="site", columns="year", values="dens").reindex(sites)[years]
    pred = fld.pivot(index="site", columns="year", values="pred").reindex(sites)[years]
    vmax = float(np.nanpercentile(fld["dens"], 98))
    ext = [years[0] - 0.5, years[-1] + 0.5, len(sites) - 0.5, -0.5]

    for col, (M, ttl) in enumerate([(obs, "(b)  Observed (gaps = unsampled)"),
                                    (pred, "(c)  stRF reconstruction (all years)")]):
        ax = fig.add_subplot(gs[1, col])
        im = ax.imshow(M.to_numpy(), aspect="auto", cmap="YlGnBu", vmin=0, vmax=vmax,
                       extent=ext, interpolation="nearest")
        ax.set_title(ttl, fontsize=9.5, loc="left")
        ax.set_xlabel("Year")
        # park separators + labels
        bounds, lab_pos, cum = [], [], 0
        for rg in REG_ORDER:
            n = (order["region"] == rg).sum()
            if n == 0:
                continue
            lab_pos.append((rg, cum + n / 2)); cum += n; bounds.append(cum)
        for b in bounds[:-1]:
            ax.axhline(b - 0.5, color="white", lw=1.2)
        if col == 0:
            ax.set_yticks([p for _, p in lab_pos])
            ax.set_yticklabels([REG_LAB[r] for r, _ in lab_pos], fontsize=7)
        else:
            ax.set_yticks([])
        ax.set_facecolor("0.85")
    cax = fig.add_axes([0.40, 0.055, 0.22, 0.014])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("Shoot density 0.04 m$^{-2}$", fontsize=7.5); cb.ax.tick_params(labelsize=6.5)
    fs.save(fig, f"{F}/FigureT1_spatiotemporal_reconstruction.png")


def figT2():
    sk = pd.read_csv(f"{R}/temporal_skill.csv")
    f = pd.read_csv(f"{R}/temporal_forecast_raw.csv")
    ag = pd.read_csv(f"{R}/temporal_reconstruction.csv")   # has site-level obs
    siteclim = pd.read_csv(f"{R}/temporal_field.csv").groupby("site")["dens"].mean()
    f["obs_anom"] = f["obs"] - f["site"].map(siteclim)
    f["strf_anom"] = f["stRF"] - f["site"].map(siteclim)

    MC = {"stRF": "#D55E00", "persistence": "#0072B2", "climatology": "#666666"}
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 4.1))
    fig.subplots_adjust(left=0.08, right=0.985, top=0.9, bottom=0.14, wspace=0.26)

    # (a) skill vs horizon (Spearman)
    for m in ["stRF", "persistence", "climatology"]:
        d = sk[sk["method"] == m].sort_values("horizon")
        ax[0].plot(d["horizon"], d["spearman"], "-o", color=MC[m], ms=5, lw=1.8, label=m)
    ax[0].axvspan(-0.5, 0.5, color="0.9", zorder=0)
    ax[0].text(0.0, ax[0].get_ylim()[0], " interp.", fontsize=6.5, va="bottom", ha="center", color="0.4")
    ax[0].set_xlabel("Forecast horizon (years ahead)"); ax[0].set_ylabel("Spearman (predicted vs observed)")
    ax[0].set_title("(a)  Level skill: stRF ties the naive baselines", fontsize=9.5, loc="left")
    ax[0].legend(fontsize=8, loc="lower left")
    ax[0].set_xticks(range(0, 6))

    # (b) anomaly correlation vs horizon (skill at predicting CHANGE)
    from scipy.stats import pearsonr
    hs, rr = [], []
    for h in sorted(f["horizon"].unique()):
        s = f[f["horizon"] == h].dropna(subset=["obs_anom", "strf_anom"])
        if len(s) < 30:
            continue
        hs.append(h); rr.append(pearsonr(s["strf_anom"], s["obs_anom"])[0])
    ax[1].axhline(0, color="0.6", lw=1, ls="--")
    ax[1].plot(hs, rr, "-o", color="#D55E00", ms=5, lw=1.8)
    ax[1].fill_between(hs, rr, 0, where=[v < 0 for v in rr], color="#D55E00", alpha=0.12)
    ax[1].set_xlabel("Forecast horizon (years ahead)")
    ax[1].set_ylabel("Corr(predicted, observed) anomaly")
    ax[1].set_title("(b)  Change skill: stRF cannot predict departures", fontsize=9.5, loc="left")
    ax[1].set_ylim(-0.35, 0.35); ax[1].set_xticks(range(1, 6))
    ax[1].text(0.5, 0.97, "zero = no skill at predicting change", transform=ax[1].transAxes,
               fontsize=7, va="top", ha="center", color="0.4")
    fs.save(fig, f"{F}/FigureT2_forecast_skill.png")


if __name__ == "__main__":
    figT1(); figT2()
    print("temporal figures done")
