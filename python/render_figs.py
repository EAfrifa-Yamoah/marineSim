"""
render_figs.py — regenerate every figure whose content changed.

Style: Okabe-Ito colour-blind-safe palette; no titles embedded in figures;
one shared legend per figure; captions live in the document text.
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = "../results"
OK = dict(blue="#0072B2", orange="#E69F00", green="#009E73",
          verm="#D55E00", purple="#CC79A7", grey="#999999", sky="#56B4E9")
plt.rcParams.update({"font.family": "serif", "font.size": 9,
                     "axes.linewidth": 0.7})


# ---------------------------------------------------- Figure 2 (decomposition)
def figure2():
    d = json.load(open(f"{R}/decomposition_corrected.json"))
    df = pd.read_csv(f"{R}/ablation_decomposed.csv")
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 3.55))

    labels = ["Algorithmic\nflexibility", "Coordinates", "Distance\nfields",
              "Spatial\nlags", "Temporal\nfeatures"]
    keys = ["ML", "COORD", "EDF", "LAG", "TEMP"]
    cols = [OK["grey"], OK["blue"], OK["blue"], OK["blue"], OK["orange"]]
    est = [d["shapley"][k][0] for k in keys]
    lo = [d["shapley"][k][0] - d["shapley"][k][1] for k in keys]
    hi = [d["shapley"][k][2] - d["shapley"][k][0] for k in keys]
    total = d["shapley_total"]

    ax[0].bar(range(5), est, yerr=[lo, hi], capsize=2.5, color=cols,
              edgecolor="black", linewidth=0.5, width=0.68)
    for i, v in enumerate(est):
        ax[0].text(i, v + hi[i] + 0.003, f"{100*v/total:.0f}%",
                   ha="center", va="bottom", fontsize=7.5)
    ax[0].set_xticks(range(5))
    ax[0].set_xticklabels(labels, fontsize=6.8)
    ax[0].set_ylabel("Shapley contribution to test AUC", fontsize=8)
    ax[0].set_ylim(0, max(est) + 0.020)
    ax[0].tick_params(labelsize=7.5)
    ax[0].text(0.02, 0.96, "(a)", transform=ax[0].transAxes,
               fontweight="bold", va="top", fontsize=9)
    ax[0].grid(axis="y", alpha=0.25, lw=0.4); ax[0].set_axisbelow(True)

    g = df.groupby("n_time")[["published_d_temporal", "temporal_last"]].agg(
        ["mean", "sem"])
    T = g.index.values; x = np.arange(len(T)); w = 0.34
    ax[1].bar(x - w/2, g[("published_d_temporal", "mean")],
              w, yerr=1.96*g[("published_d_temporal", "sem")], capsize=2.5,
              color=OK["verm"], edgecolor="black", linewidth=0.5,
              label="Nested contrast (stRF − Spatial RF)")
    ax[1].bar(x + w/2, g[("temporal_last", "mean")],
              w, yerr=1.96*g[("temporal_last", "sem")], capsize=2.5,
              color=OK["orange"], edgecolor="black", linewidth=0.5,
              label="Isolated temporal increment")
    ax[1].axhline(0, color="black", lw=0.7)
    ax[1].set_xticks(x); ax[1].set_xticklabels([str(t) for t in T])
    ax[1].set_xlabel("Number of timesteps", fontsize=8)
    ax[1].set_ylabel("Attributed temporal gain in AUC", fontsize=8)
    ax[1].tick_params(labelsize=7.5)
    ax[1].text(0.02, 0.96, "(b)", transform=ax[1].transAxes,
               fontweight="bold", va="top", fontsize=9)
    ax[1].grid(axis="y", alpha=0.25, lw=0.4); ax[1].set_axisbelow(True)
    ax[1].annotate("no information\nat T = 1", xy=(0 + w/2, 0.002),
                   xytext=(0.30, 0.58), textcoords="axes fraction",
                   fontsize=6.8, ha="left",
                   arrowprops=dict(arrowstyle="->", lw=0.6))

    h, l = ax[1].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, frameon=False, fontsize=7.5,
               bbox_to_anchor=(0.5, -0.015))
    fig.tight_layout(rect=[0, 0.075, 1, 1])
    fig.savefig(f"{R}/FIG_decomposition.png", dpi=400, bbox_inches="tight")
    plt.close(fig)
    print("  FIG_decomposition.png")


# ------------------------------------------------- Figure 6 (case study)
def figure6():
    cv = pd.read_csv(f"{R}/case_cv_regimes.csv")
    loro = pd.read_csv(f"{R}/case_loro_corrected.csv")
    imp = pd.read_csv(f"{R}/case_importance.csv")
    drop = pd.read_csv(f"{R}/case_dropcol.csv")
    fig, ax = plt.subplots(1, 3, figsize=(9.0, 3.0))

    order = ["GLM", "RF", "SpatialRF", "stRF"]
    disp = ["GLM", "RF", "Spatial RF", "stRF"]
    rand = cv[cv.regime == "random 10-fold"].set_index("method")["Spearman"]
    site = cv[cv.regime == "site-grouped"].set_index("method")["Spearman"]
    lo = loro.groupby("method")["Spearman"].mean()
    x = np.arange(len(order)); w = 0.27
    ax[0].bar(x - w, [rand.get(m, np.nan) for m in order], w, color=OK["blue"],
              edgecolor="black", lw=0.5, label="Random 10-fold (interpolation)")
    ax[0].bar(x, [site.get(m, np.nan) for m in order], w, color=OK["sky"],
              edgecolor="black", lw=0.5, label="Site grouped (new site)")
    ax[0].bar(x + w, [lo.get(m, np.nan) for m in order], w, color=OK["orange"],
              edgecolor="black", lw=0.5, label="Leave one region out (new region)")
    ax[0].axhline(0, color="black", lw=0.7)
    ax[0].set_xticks(x); ax[0].set_xticklabels(disp, fontsize=7.5)
    ax[0].set_ylabel("Spearman rank correlation", fontsize=8)
    ax[0].tick_params(labelsize=7.5)
    ax[0].text(0.02, 0.96, "(a)", transform=ax[0].transAxes,
               fontweight="bold", va="top")
    ax[0].grid(axis="y", alpha=0.25, lw=0.4); ax[0].set_axisbelow(True)

    imp = imp.sort_values(imp.columns[1], ascending=True).tail(8)
    ax[1].barh(range(len(imp)), imp.iloc[:, 1], color=OK["blue"],
               edgecolor="black", lw=0.5)
    ax[1].set_yticks(range(len(imp)))
    ax[1].set_yticklabels(imp.iloc[:, 0], fontsize=6.8)
    ax[1].set_xlabel("Permutation importance", fontsize=8)
    ax[1].tick_params(labelsize=7.5)
    ax[1].text(0.02, 0.96, "(b)", transform=ax[1].transAxes,
               fontweight="bold", va="top")
    ax[1].grid(axis="x", alpha=0.25, lw=0.4); ax[1].set_axisbelow(True)

    drop = drop.sort_values("unique_loss", ascending=True).tail(8)
    ax[2].barh(range(len(drop)), drop["unique_loss"], color=OK["green"],
               edgecolor="black", lw=0.5)
    ax[2].set_yticks(range(len(drop)))
    ax[2].set_yticklabels(drop.iloc[:, 0], fontsize=6.8)
    ax[2].set_xlabel("Unique contribution (drop-column)", fontsize=8)
    ax[2].tick_params(labelsize=7.5)
    ax[2].text(0.02, 0.96, "(c)", transform=ax[2].transAxes,
               fontweight="bold", va="top")
    ax[2].grid(axis="x", alpha=0.25, lw=0.4); ax[2].set_axisbelow(True)

    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=7.5,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.09, 1, 1])
    fig.savefig(f"{R}/FIG_case_analysis.png", dpi=400, bbox_inches="tight")
    plt.close(fig)
    print("  FIG_case_analysis.png")


# ------------------------------------------- Figure S6 (forecast skill)
def figS_forecast():
    ro = pd.read_csv(f"{R}/case_rolling_origin.csv")
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.9))
    lv = ro.pivot_table(index="horizon", columns="method", values="rho_level")
    style = {"stRF": (OK["orange"], "o", "-"),
             "climatology": (OK["blue"], "s", "--"),
             "persistence": (OK["green"], "^", ":")}
    for m in ["stRF", "climatology", "persistence"]:
        c, mk, ls = style[m]
        ax[0].plot(lv.index, lv[m], marker=mk, color=c, ls=ls, lw=1.2,
                   ms=4.5, label=m if m != "stRF" else "stRF")
    ax[0].set_xlabel("Forecast horizon (years)", fontsize=8)
    ax[0].set_ylabel("Rank correlation, density level", fontsize=8)
    ax[0].tick_params(labelsize=7.5)
    ax[0].set_ylim(0.55, 0.82)
    ax[0].text(0.02, 0.96, "(a)", transform=ax[0].transAxes,
               fontweight="bold", va="top")
    ax[0].grid(alpha=0.25, lw=0.4); ax[0].set_axisbelow(True)

    an = ro[ro.method == "stRF"].groupby("horizon")["rho_anomaly"].mean()
    ax[1].bar(an.index.astype(str), an.values, color=OK["orange"],
              edgecolor="black", lw=0.5, width=0.6)
    ax[1].axhline(0, color="black", lw=0.7)
    ax[1].set_xlabel("Forecast horizon (years)", fontsize=8)
    ax[1].set_ylabel("Rank correlation, anomaly", fontsize=8)
    ax[1].set_ylim(-0.05, 0.25)
    ax[1].tick_params(labelsize=7.5)
    ax[1].text(0.02, 0.96, "(b)", transform=ax[1].transAxes,
               fontweight="bold", va="top")
    ax[1].grid(axis="y", alpha=0.25, lw=0.4); ax[1].set_axisbelow(True)

    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=7.5,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.10, 1, 1])
    fig.savefig(f"{R}/FIG_forecast_skill.png", dpi=400, bbox_inches="tight")
    plt.close(fig)
    print("  FIG_forecast_skill.png")


# --------------------------------------- Figure S7 (interval coverage)
def figS_coverage():
    iv = pd.read_csv(f"{R}/case_intervals.csv")
    g = iv.groupby(["regime", "interval"])[["coverage", "width"]].mean().reset_index()
    order = ["ensemble spread", "QRF", "split conformal"]
    disp = ["Ensemble\nspread", "Quantile\nregression forest", "Split\nconformal"]
    regimes = ["interpolation (random 10-fold)", "extrapolation (LORO)"]
    rlab = ["Interpolation", "Extrapolation (LORO)"]
    cols = [OK["blue"], OK["orange"]]
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.9))
    x = np.arange(3); w = 0.36
    for k, (rg, lab) in enumerate(zip(regimes, rlab)):
        s = g[g.regime == rg].set_index("interval")
        ax[0].bar(x + (k - 0.5) * w, [s.loc[o, "coverage"] for o in order], w,
                  color=cols[k], edgecolor="black", lw=0.5, label=lab)
        ax[1].bar(x + (k - 0.5) * w, [s.loc[o, "width"] for o in order], w,
                  color=cols[k], edgecolor="black", lw=0.5)
    ax[0].axhline(0.80, color="black", ls="--", lw=0.9)
    ax[0].text(2.42, 0.815, "nominal", fontsize=6.5, ha="right")
    ax[0].set_xticks(x); ax[0].set_xticklabels(disp, fontsize=6.8)
    ax[0].set_ylabel("Empirical coverage", fontsize=8)
    ax[0].set_ylim(0, 1.0); ax[0].tick_params(labelsize=7.5)
    ax[0].text(0.02, 0.96, "(a)", transform=ax[0].transAxes,
               fontweight="bold", va="top")
    ax[0].grid(axis="y", alpha=0.25, lw=0.4); ax[0].set_axisbelow(True)

    ax[1].set_xticks(x); ax[1].set_xticklabels(disp, fontsize=6.8)
    ax[1].set_ylabel("Mean interval width (shoots per 0.04 m²)", fontsize=8)
    ax[1].tick_params(labelsize=7.5)
    ax[1].text(0.02, 0.96, "(b)", transform=ax[1].transAxes,
               fontweight="bold", va="top")
    ax[1].grid(axis="y", alpha=0.25, lw=0.4); ax[1].set_axisbelow(True)

    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, frameon=False, fontsize=7.5,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.10, 1, 1])
    fig.savefig(f"{R}/FIG_interval_coverage.png", dpi=400, bbox_inches="tight")
    plt.close(fig)
    print("  FIG_interval_coverage.png")


if __name__ == "__main__":
    figure2()
    try:
        figure6()
    except Exception as e:
        print("  figure6 skipped:", e)
    figS_forecast()
    figS_coverage()
