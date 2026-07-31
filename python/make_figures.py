"""
make_figures.py - render the six manuscript figures (300 dpi PNG) into figures/.
"""
import os, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
import figstyle as fs

R = os.path.join(os.path.dirname(__file__), "..", "results")
F = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(F, exist_ok=True)
fs.setup()
D = json.load(open(f"{R}/derived.json"))
C = json.load(open(f"{R}/case_derived.json"))
CE = json.load(open(f"{R}/case_ext.json"))
METHODS = ["GLM", "GAM", "RF", "SpatialRF", "stRF"]


# ============================================================ Figure 1
def fig1_decomposition():
    dec = D["decomposition"]; sh = D["shares"]
    comps = [("ml", "Algorithmic\nflexibility\n(RF \u2212 GLM)"),
             ("spatial", "Spatial structure\n(SpatialRF \u2212 RF)"),
             ("temporal", "Temporal features\n(stRF \u2212 SpatialRF)")]
    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    y = np.arange(len(comps))[::-1]
    for yi, (key, lab) in zip(y, comps):
        m, lo, hi = dec[key]["mean"], dec[key]["lo"], dec[key]["hi"]
        ax.barh(yi, m, color=fs.COMP_COLOR[key], height=0.55,
                edgecolor="black", linewidth=0.8, zorder=3)
        ax.plot([lo, hi], [yi, yi], color="black", lw=1.4, zorder=4)
        for b in (lo, hi):
            ax.plot([b, b], [yi - 0.09, yi + 0.09], color="black", lw=1.4, zorder=4)
        ax.annotate(f"+{m:.3f} AUC   ({sh[key]*100:.0f}% of total)",
                    xy=(hi, yi), xytext=(8, 0), textcoords="offset points",
                    va="center", ha="left", fontsize=9.5, zorder=5)
    ax.set_yticks(y); ax.set_yticklabels([l for _, l in comps], fontsize=9.5)
    tot = dec["total"]["mean"]
    ax.axvline(0, color="black", lw=0.9)
    ax.set_xlabel("Gain in test AUC over the non-spatial GLM baseline")
    ax.set_xlim(-0.01, 0.205)
    _ = tot
    ax.grid(axis="y", visible=False)
    fs.save(fig, f"{F}/Figure1_decomposition.png")


# ============================================================ Figure 2
def fig2_auc_by_n():
    df = pd.read_csv(f"{R}/fig2_auc_by_n.csv")
    thr = D["thresholds"]
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ns = sorted(df["n"].unique())
    for m in METHODS:
        g = df[df.method == m].sort_values("n")
        se = g["sd"] / np.sqrt(36)
        ax.plot(g["n"], g["mean"], marker=fs.METHOD_MARKER[m], color=fs.METHOD_COLOR[m],
                lw=2.0 if m == "stRF" else 1.4, ms=7 if m == "stRF" else 5.5,
                label=m, zorder=4 if m == "stRF" else 3)
        ax.fill_between(g["n"], g["mean"] - se, g["mean"] + se,
                        color=fs.METHOD_COLOR[m], alpha=0.13, zorder=1)
    ax.axhline(0.75, color="0.4", ls="--", lw=1.0, zorder=2)
    ax.text(ns[0]*0.96, 0.756, "AUC = 0.75  (useful-skill threshold)",
            fontsize=8.3, color="0.35", va="bottom", ha="left")
    ax.axhline(0.5, color="0.6", ls=":", lw=0.8, zorder=2)
    ax.set_xscale("log"); ax.set_xticks(ns); ax.set_xticklabels(ns)
    ax.set_xlabel("Number of sampled sites (log scale)")
    ax.set_ylabel("Test AUC")
    ax.set_ylim(0.48, 0.95)
    # threshold annotations for stRF and aspatial RF (spaced apart)
    if thr.get("stRF"):
        ax.annotate(f"stRF reaches 0.75\nat n \u2248 {thr['stRF']}",
                    xy=(thr["stRF"], 0.75), xytext=(thr["stRF"]*1.05, 0.86),
                    fontsize=8.5, color=fs.METHOD_COLOR["stRF"], ha="left",
                    arrowprops=dict(arrowstyle="->", color=fs.METHOD_COLOR["stRF"], lw=1))
    if thr.get("RF"):
        ax.annotate(f"aspatial RF needs\nn \u2248 {thr['RF']}",
                    xy=(thr["RF"], 0.752), xytext=(thr["RF"]*0.40, 0.60),
                    fontsize=8.5, color=fs.METHOD_COLOR["RF"], ha="left",
                    arrowprops=dict(arrowstyle="->", color=fs.METHOD_COLOR["RF"], lw=1))
    ax.legend(loc="lower right", ncol=2, columnspacing=1.0, handletextpad=0.5)
    fs.save(fig, f"{F}/Figure2_auc_by_samplesize.png")


# ============================================================ Figure 3
def fig3_temporal_surface():
    df = pd.read_csv(f"{R}/fig3_temporal_surface.csv")
    piv = df.pivot(index="n", columns="n_time", values="d_temporal").sort_index(ascending=False)
    fig, ax = plt.subplots(figsize=(5.7, 4.4))
    im = ax.imshow(piv.values, cmap="YlGnBu", aspect="auto", vmin=0,
                   vmax=float(piv.values.max())*1.02)
    ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns)
    ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index)
    ax.set_xlabel("Number of time points"); ax.set_ylabel("Number of sampled sites")
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            ax.text(j, i, f"+{v:.3f}", ha="center", va="center",
                    color="white" if v > piv.values.max()*0.55 else "black", fontsize=9)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("\u0394AUC", fontsize=9)
    ax.set_xticks(np.arange(-.5, piv.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, piv.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.2); ax.grid(which="major", visible=False)
    ax.tick_params(which="minor", length=0)
    fs.save(fig, f"{F}/Figure3_temporal_surface.png")


# ============================================================ Figure 4
def fig4_design():
    df = pd.read_csv(f"{R}/fig4_design_n100.csv").set_index("design")
    designs = ["random", "clustered", "stratified"]
    dcol = {"random": fs.OI["skyblue"], "clustered": fs.OI["vermillion"],
            "stratified": fs.OI["green"]}
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    x = np.arange(len(METHODS)); w = 0.26
    for k, des in enumerate(designs):
        vals = [df.loc[des, m] for m in METHODS]
        ax.bar(x + (k-1)*w, vals, w, label=des.capitalize(), color=dcol[des],
               edgecolor="black", linewidth=0.6, zorder=3)
    ax.axhline(0.75, color="0.4", ls="--", lw=1.0, zorder=2)
    ax.set_xticks(x); ax.set_xticklabels(METHODS)
    ax.set_ylabel("Test AUC"); ax.set_ylim(0.5, 0.92)
    ax.set_xlabel("Method")
    ax.legend(loc="upper left", title="Design", ncol=3, columnspacing=1.0)
    ax.grid(axis="x", visible=False)
    fs.save(fig, f"{F}/Figure4_sampling_design.png")


# ============================================================ Figure 5
def fig5_case_study():
    grid = pd.read_csv(f"{R}/case_grid_unc.csv").sort_values("latitude")
    reg = pd.read_csv(f"{R}/case_regions.csv")
    code = {"Jurien Bay Marine Park": "JBMP", "Marmion Marine Park": "MMP",
            "Cockburn Sound": "CSMC", "Shoalwater Islands Marine Park": "SIMP",
            "Ngari Capes Marine Park": "NCMP"}
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    lat = grid["latitude"].values; pred = grid["pred_mean"].values
    plo = grid["pred_lo"].values; phi = grid["pred_hi"].values
    inside = grid["inside_aoa"].values.astype(bool)
    ymax = float(np.nanmax(phi))
    top = ymax * 1.30
    # outside-AOA bands
    outside = ~inside
    if outside.any():
        edges = np.where(np.diff(outside.astype(int)) != 0)[0]
        idx = np.r_[0, edges + 1, len(outside)]
        first = True
        for a, b in zip(idx[:-1], idx[1:]):
            if outside[a]:
                ax.axvspan(lat[a], lat[min(b, len(lat)-1)], color=fs.OI["vermillion"],
                           alpha=0.14, zorder=1,
                           label="Outside area of applicability" if first else None)
                first = False
    ax.fill_between(lat, plo, phi, color=fs.OI["blue"], alpha=0.18, zorder=3,
                    label="stRF 80% prediction band")
    ax.plot(lat, pred, color=fs.OI["blue"], lw=2.3, zorder=4,
            label="stRF predicted shoot density")

    # region centroids: dashed line + staggered code label with leader
    reg = reg.sort_values("lat", ascending=False)  # north -> south
    # explicit label latitudes (spread) and heights to avoid collisions
    place = {"JBMP": (-30.10, top*0.99), "MMP": (-31.55, top*0.99),
             "CSMC": (-32.05, top*0.88), "SIMP": (-32.62, top*0.99),
             "NCMP": (-33.58, top*0.99)}
    for _, rr in reg.iterrows():
        cd = code[rr["name"]]; lx, ly = place[cd]
        ax.axvline(rr["lat"], color="0.6", ls=(0, (4, 3)), lw=0.8, zorder=2)
        ax.scatter(rr["lat"], 0, s=55, color="black", zorder=6, clip_on=False)
        ax.annotate(cd, xy=(rr["lat"], ymax*1.02), xytext=(lx, ly),
                    ha="center", va="center", fontsize=8.6, fontweight="bold",
                    color="0.15", zorder=7,
                    arrowprops=dict(arrowstyle="-", color="0.55", lw=0.8,
                                    connectionstyle="arc3,rad=0.0"))
    ax.set_xlabel("Latitude (\u00b0S)")
    ax.set_ylabel("Predicted mean shoot density\n(shoots 0.04 m$^{-2}$)")
    ax.set_xlim(lat.max()+0.03, lat.min()-0.03)   # inverted: north (left) -> south (right)
    ax.set_ylim(0, top)
    # direction note
    ax.annotate("north", xy=(0.015, 1.005), xycoords="axes fraction", fontsize=8.5,
                color="0.4", ha="left", va="bottom")
    ax.annotate("south", xy=(0.985, 1.005), xycoords="axes fraction", fontsize=8.5,
                color="0.4", ha="right", va="bottom")
    # monitoring footprint table (compact, in open lower band)
    rows = [f"{code[r['name']]:5s} {int(r['nyears']):>2d} yr   {int(r['nsite']):>2d} sites"
            for _, r in reg.iterrows()]
    txt = "Monitoring footprint\n" + "\n".join(rows)
    ax.text(0.985, 0.04, txt, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7.8, family="DejaVu Sans Mono",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.7", lw=0.8))
    ax.legend(loc="lower left", fontsize=8.4)
    fs.save(fig, f"{F}/Figure5_seagrass_latitudinal.png")


# ============================================================ Figure 6
def _plab(ax, s):
    """panel label placed just outside the top-left of the axes (never over data)."""
    ax.text(-0.02, 1.045, s, transform=ax.transAxes, fontsize=12, fontweight="bold",
            ha="right", va="bottom")

def fig6_case_analysis():
    """(a) interpolation vs extrapolation skill; (b) permutation importance (model
    reliance); (c) drop-column unique contribution, exposing spatial redundancy."""
    ie = pd.read_csv(f"{R}/case_interp_extrap.csv")
    dc = pd.read_csv(f"{R}/case_dropcol.csv")
    meta = json.load(open(f"{R}/case_maps_meta.json"))
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(13.6, 4.5),
                                        gridspec_kw={"width_ratios": [1.0, 1.05, 1.05]})
    spatial = {"Spatial lags", "Distance fields", "Latitude", "Longitude"}
    fcol = lambda f: (fs.OI["blue"] if f in spatial else
                      (fs.OI["green"] if f == "Depth" else fs.OI["grey"]))

    # ----- (a) interpolation vs extrapolation -----
    methods = ["GLM", "RF", "SpatialRF", "stRF"]
    regimes = ["Interpolation (random CV)", "Extrapolation (leave-region-out)"]
    rcol = {regimes[0]: fs.OI["skyblue"], regimes[1]: fs.OI["vermillion"]}
    x = np.arange(len(methods)); w = 0.38
    for k, reg in enumerate(regimes):
        vals = [ie[(ie.method == m) & (ie.regime == reg)]["Spearman"].iloc[0] for m in methods]
        axA.bar(x + (k - 0.5) * w, vals, w, label=reg.split(" (")[0],
                color=rcol[reg], edgecolor="black", linewidth=0.6, zorder=3)
    axA.axhline(0, color="black", lw=0.8)
    axA.set_xticks(x); axA.set_xticklabels(methods)
    axA.set_ylabel("Spearman rank correlation\n(predicted vs observed)")
    axA.set_xlabel("Method"); axA.set_ylim(-0.12, 0.95)
    axA.legend(loc="upper left", fontsize=8.2); axA.grid(axis="x", visible=False)
    _plab(axA, "(a)")

    # shared feature order (by permutation reliance, ascending for horizontal bars)
    dc = dc.sort_values("perm_pct")
    yb = np.arange(len(dc)); cols = [fcol(f) for f in dc["feature"]]

    # ----- (b) permutation importance (reliance) -----
    axB.barh(yb, dc["perm_pct"], color=cols, edgecolor="black", linewidth=0.6, zorder=3)
    axB.set_yticks(yb); axB.set_yticklabels(dc["feature"], fontsize=9)
    axB.set_xlabel("Permutation importance (% of total)")
    for i, v in enumerate(dc["perm_pct"]):
        axB.text(v + 0.6, i, f"{v:.0f}%", va="center", fontsize=8)
    axB.set_xlim(0, dc["perm_pct"].max() * 1.20); axB.grid(axis="y", visible=False)
    handles = [Patch(facecolor=fs.OI["blue"], edgecolor="black", label="Spatial features"),
               Patch(facecolor=fs.OI["green"], edgecolor="black", label="Depth"),
               Patch(facecolor=fs.OI["grey"], edgecolor="black", label="Other env. / temporal")]
    axB.legend(handles=handles, loc="lower right", fontsize=7.6)
    _plab(axB, "(b)")
    axB.set_title("Model reliance", fontsize=9.5, color="0.3")

    # ----- (c) drop-column unique contribution -----
    loss = dc["unique_loss"].clip(lower=0)
    axC.barh(yb, loss, color=cols, edgecolor="black", linewidth=0.6, zorder=3)
    axC.set_yticks(yb); axC.set_yticklabels([])
    axC.set_xlabel("Unique contribution\n(Spearman loss when removed)")
    for i, v in enumerate(dc["unique_loss"]):
        axC.text(max(v, 0) + 0.0008, i, f"{v:+.3f}", va="center", fontsize=8)
    axC.set_xlim(0, max(loss.max(), 0.001) * 1.5); axC.grid(axis="y", visible=False)
    axC.text(0.97, 0.30, f"Spatial features\njointly: +{meta['all_spatial_loss']:.3f}",
             transform=axC.transAxes, ha="right", va="center", fontsize=8.2,
             color=fs.OI["blue"], fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=fs.OI["blue"], lw=1))
    _plab(axC, "(c)")
    axC.set_title("Unique information", fontsize=9.5, color="0.3")

    fig.subplots_adjust(wspace=0.30, top=0.90, bottom=0.16, left=0.07, right=0.985)
    fs.save(fig, f"{F}/Figure6_case_analysis.png")


def fig7_spatial_maps():
    """Three coastal maps: (a) observed shoot density at the monitoring sites;
    (b) stRF predicted density over a continuous coastal corridor; (c) the area of
    applicability, flagging the unmonitored latitudinal gaps."""
    site = pd.read_csv(f"{R}/case_sites.csv")
    grid = pd.read_csv(f"{R}/case_map_grid.csv")
    coast = json.load(open(f"{R}/coast_wa.json"))   # WA land polygons (geo-maps 100m)
    # region column already holds short codes (JBMP, MMP, ...); Okabe-Ito colours
    rc = {"JBMP": fs.OI["orange"], "MMP": fs.OI["skyblue"], "CSMC": fs.OI["blue"],
          "SIMP": fs.OI["green"], "NCMP": fs.OI["vermillion"]}
    site["code"] = site["region"]
    aspect = 1.0 / np.cos(np.radians(32))
    lon0, lon1 = 114.83, 115.97
    lat0, lat1 = site["latitude"].min() - 0.12, site["latitude"].max() + 0.12

    def draw_coast(ax):
        ax.add_collection(PolyCollection(coast, facecolors="#EAE6DA",
                          edgecolors="#8d8674", linewidths=0.45, zorder=0.6))

    fig, axes = plt.subplots(1, 3, figsize=(7.6, 7.4), sharey=True)
    for ax in axes:
        ax.set_xlim(lon0, lon1); ax.set_ylim(lat0, lat1)   # north-up (less-negative lat at top)
        ax.set_aspect(aspect); ax.set_xlabel("Longitude (\u00b0E)")
        ax.set_xticks([115.0, 115.4, 115.8]); draw_coast(ax)
    axes[0].set_ylabel("Latitude (\u00b0S)")

    # marker size to tile the corridor ribbon as a quasi-continuous surface
    ms = 26

    # ---- (a) observed sites ----
    axa = axes[0]
    smax = site["dens"].max()
    for cd, sub in site.groupby("code"):
        axa.scatter(sub["longitude"], sub["latitude"], s=18 + 90 * sub["dens"] / smax,
                    color=rc[cd], edgecolor="black", linewidth=0.5, alpha=0.9, zorder=4)
    # park code labels at cluster centroids (staggered to avoid collisions)
    yoff = {"NCMP": 0, "SIMP": -14, "CSMC": 14, "MMP": 0, "JBMP": 0}
    for cd, sub in site.groupby("code"):
        axa.annotate(cd, (sub["longitude"].mean(), sub["latitude"].mean()),
                     xytext=(-30, yoff.get(cd, 0)), textcoords="offset points", fontsize=8,
                     fontweight="bold", color="0.15", ha="right", va="center")
    # density size legend (proxy handles), placed in the empty centre-left
    for d in [15, 30]:
        axa.scatter([], [], s=18 + 90 * d / smax, color="0.6", edgecolor="black",
                    linewidth=0.5, label=f"{d} shoots")
    axa.legend(loc="center left", bbox_to_anchor=(0.0, 0.63), fontsize=7.0,
               title="Mean density", title_fontsize=7.4, framealpha=0.92,
               labelspacing=1.0, borderpad=0.8)
    _plab(axa, "(a)")
    axa.set_title("Observed shoot density\n(marker size \u221d density)", fontsize=9.5)

    # ---- (b) predicted surface ----
    # Colour limits are taken from the data (observed and predicted pooled) rather than
    # hard coded, so the surface and the observations are read on one scale. Observed
    # sites are drawn as circles filled on that same scale: where the model agrees they
    # blend into the surface, where it disagrees they stand out.
    axb = axes[1]
    vlo = float(np.floor(min(site["dens"].min(), grid["pred"].min())))
    vhi = float(np.ceil(max(site["dens"].quantile(0.98), grid["pred"].max())))
    scb = axb.scatter(grid["longitude"], grid["latitude"], c=grid["pred"], cmap="YlGnBu",
                      s=ms, marker="s", edgecolors="none", vmin=vlo, vmax=vhi, zorder=2)
    axb.scatter(site["longitude"], site["latitude"], c=site["dens"], cmap="YlGnBu",
                vmin=vlo, vmax=vhi, s=30, marker="o", edgecolor="black",
                linewidth=0.6, zorder=4)
    # Isobaths of the interpolated site depth (NOT bathymetry). The corridor depth is IDW
    # from the monitoring-site depths, spanning ~1-18 m, so contours are drawn only at 5,
    # 10 and 15 m. Triangles that span the inter-park gaps are masked so the contours do
    # not bridge unsampled water.
    from matplotlib.tri import Triangulation
    gx, gy, gz = (grid["longitude"].to_numpy(), grid["latitude"].to_numpy(),
                  grid["depth"].to_numpy())
    tri = Triangulation(gx, gy)
    sxr = np.cos(np.radians(np.mean(gy))); T = tri.triangles
    def _el(a, b): return np.sqrt((sxr * (gx[a] - gx[b])) ** 2 + (gy[a] - gy[b]) ** 2)
    emax = np.maximum.reduce([_el(T[:, 0], T[:, 1]), _el(T[:, 1], T[:, 2]), _el(T[:, 0], T[:, 2])])
    tri.set_mask(emax > np.percentile(emax, 92))
    iso = axb.tricontour(tri, gz, levels=[-15, -10, -5], colors="0.15",
                         linewidths=0.6, alpha=0.85, zorder=3)
    axb.clabel(iso, fmt={-5: "5 m", -10: "10 m", -15: "15 m"}, fontsize=5.2,
               inline=True, inline_spacing=2)
    obs_h = Line2D([], [], marker="o", ls="", markerfacecolor="0.75",
                   markeredgecolor="black", markersize=6, label="Observed site")
    iso_h = Line2D([], [], color="0.15", lw=0.8, label="Depth (m)")
    axb.legend(handles=[obs_h, iso_h], loc="upper right", bbox_to_anchor=(0.995, 0.995),
               fontsize=6.8, framealpha=0.92, borderpad=0.5, handletextpad=0.5)
    # Colour bar as a vertical inset sitting directly BELOW the "Observed site" legend,
    # over the land side of the panel where there is no data. An attached colour bar
    # (fig.colorbar(ax=axb)) would steal width from this panel alone and leave the three
    # panels unevenly spaced.
    cax = axb.inset_axes([0.845, 0.60, 0.05, 0.25])
    cb = fig.colorbar(scb, cax=cax, orientation="vertical")
    cb.ax.yaxis.set_ticks_position("left")
    cb.ax.tick_params(labelsize=6.2, pad=1.5, length=2)
    cb.ax.set_title("Shoots\n0.04 m$^{-2}$", fontsize=6.4, pad=3)
    cb.outline.set_linewidth(0.5)
    _plab(axb, "(b)")
    axb.set_title("stRF predicted density\n(circles = observed, same scale)", fontsize=9.5)

    # ---- (c) area of applicability ----
    axc = axes[2]
    ins = grid[grid["inside"] == 1]; out = grid[grid["inside"] == 0]
    axc.scatter(ins["longitude"], ins["latitude"], s=ms, marker="s",
                color=fs.OI["skyblue"], edgecolors="none", zorder=2, label="Inside AOA")
    axc.scatter(out["longitude"], out["latitude"], s=ms, marker="s",
                color=fs.OI["vermillion"], edgecolors="none", zorder=3,
                label="Outside AOA (extrapolation)")
    axc.scatter(site["longitude"], site["latitude"], s=7, color="black", zorder=4)
    axc.legend(loc="upper right", fontsize=7.0, framealpha=0.92)
    _plab(axc, "(c)")
    axc.set_title("Area of applicability\n(unmonitored gaps flagged)", fontsize=9.5)

    fig.subplots_adjust(wspace=0.10, top=0.92, bottom=0.09, left=0.10, right=0.98)
    fs.save(fig, f"{F}/Figure7_spatial_maps.png")


def fig8_design_overview():
    """Conceptual overview of the simulation design: the six ground-truth worlds
    (stationarity x autocorrelation range), the three sampling designs, and the
    sample-size gradient. Single realisations, shown for illustration. The worlds
    are rendered on a finer display grid than the benchmark (same physical extent
    and correlation ranges) purely so the continuous fields read smoothly."""
    import sim_engine as se
    from matplotlib.gridspec import GridSpec
    ranges = [30, 80, 200]
    rlabels = ["30 km (short)", "80 km (medium)", "200 km (long)"]

    g0, c0 = se.GRID, se.CELL_KM                     # finer grid for display only
    se.GRID = 160; se.CELL_KM = se.EXTENT_KM / se.GRID
    try:
        worlds = {}
        for i, st in enumerate([True, False]):
            for j, rk in enumerate(ranges):
                worlds[(i, j)] = se.World(stationary=st, range_km=rk, n_time=1,
                                          seed=100 + i * 10 + j)
        vmax = max(float(worlds[k].p[0].max()) for k in worlds)
        ext = [0, se.EXTENT_KM, 0, se.EXTENT_KM]

        fig = plt.figure(figsize=(9.6, 11.8))
        gs = GridSpec(4, 3, figure=fig, left=0.09, right=0.88, top=0.955, bottom=0.075,
                      hspace=0.22, wspace=0.10)

        im = None
        for i, stlab in enumerate(["Stationary", "Nonstationary"]):
            for j in range(3):
                ax = fig.add_subplot(gs[i, j])
                im = ax.imshow(worlds[(i, j)].p[0], origin="lower", extent=ext,
                               cmap="YlGnBu", vmin=0, vmax=vmax, aspect="equal",
                               interpolation="bilinear")
                ax.set_xticks([]); ax.set_yticks([])
                if i == 0:
                    ax.set_title(rlabels[j], fontsize=9.5)
                if j == 0:
                    ax.set_ylabel(stlab, fontsize=10.5, fontweight="bold")
        cax = fig.add_axes([0.90, 0.55, 0.016, 0.34])
        cb = fig.colorbar(im, cax=cax); cb.set_label("Occurrence probability", fontsize=8.5)
        cb.ax.tick_params(labelsize=7.5)

        wrep = worlds[(0, 1)]                 # stationary, 80 km, as the sampling backdrop
        bg, yv_all = wrep.p[0], wrep.y[0].ravel()

        def draw_samples(ax, idx, strata=None):
            ax.imshow(bg, origin="lower", extent=ext, cmap="Greys", vmin=0, vmax=vmax * 1.5,
                      alpha=0.28, aspect="equal", interpolation="bilinear")
            if strata is not None:              # show the lattice for stratified design
                for e in strata:
                    ax.axvline(e, color="0.45", lw=0.5, alpha=0.7, zorder=1.4)
                    ax.axhline(e, color="0.45", lw=0.5, alpha=0.7, zorder=1.4)
            xy = wrep.coords(idx); yv = yv_all[idx]; pres = yv == 1
            ax.scatter(xy[pres, 0], xy[pres, 1], s=15, c=fs.OI["vermillion"],
                       edgecolor="black", linewidth=0.3, zorder=3)
            ax.scatter(xy[~pres, 0], xy[~pres, 1], s=13, facecolor="white",
                       edgecolor=fs.OI["blue"], linewidth=0.5, zorder=2)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_xlim(0, se.EXTENT_KM); ax.set_ylim(0, se.EXTENT_KM)

        n_dsg = 100
        s_strata = int(np.ceil(np.sqrt(n_dsg)))
        strata_edges = np.linspace(0, se.EXTENT_KM, s_strata + 1)[1:-1]
        for j, dz in enumerate(["random", "clustered", "stratified"]):
            ax = fig.add_subplot(gs[2, j])
            draw_samples(ax, se.sample_sites(n_dsg, dz, np.random.default_rng(42)),
                         strata=strata_edges if dz == "stratified" else None)
            ax.set_title(dz.capitalize(), fontsize=9.5)
            if j == 0:
                ax.set_ylabel("Designs (n = 100)", fontsize=10.5, fontweight="bold")

        for j, nn in enumerate([30, 100, 500]):
            ax = fig.add_subplot(gs[3, j])
            draw_samples(ax, se.sample_sites(nn, "random", np.random.default_rng(7)))
            ax.set_title(f"n = {nn}", fontsize=9.5)
            if j == 0:
                ax.set_ylabel("Sample size\n(random design)", fontsize=10.5, fontweight="bold")

        pres_h = Line2D([], [], marker="o", ls="", markerfacecolor=fs.OI["vermillion"],
                        markeredgecolor="black", markersize=7, label="Sampled presence")
        abs_h = Line2D([], [], marker="o", ls="", markerfacecolor="white",
                       markeredgecolor=fs.OI["blue"], markersize=7, label="Sampled absence")
        strat_h = Line2D([], [], color="0.45", lw=0.8, label="Stratum boundaries")
        fig.legend(handles=[pres_h, abs_h, strat_h], loc="lower center", ncol=3,
                   bbox_to_anchor=(0.485, 0.015), fontsize=9, frameon=True)
        fs.save(fig, f"{F}/Figure8_design_overview.png")
    finally:
        se.GRID, se.CELL_KM = g0, c0                 # restore benchmark grid


if __name__ == "__main__":
    fig1_decomposition()
    fig2_auc_by_n()
    fig3_temporal_surface()
    fig4_design()
    fig5_case_study()
    fig6_case_analysis()
    fig7_spatial_maps()
    fig8_design_overview()
    print("ALL FIGURES DONE")
