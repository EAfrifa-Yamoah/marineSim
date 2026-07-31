"""Figure T1 (revised): latitude x time Hovmoller of stRF predicted density.

Fit stRF once on all site-years, then predict the coastal corridor for every year using
that year's ocean covariates (kd490 and the three heatwave metrics, interpolated from the
year's site values) and that year's density memory (previous-year regional mean). Collapse
the corridor to latitude bands so the prediction is a continuous latitude-by-time field.
Observed monitoring events are overlaid on the same colour scale.
"""
import io, contextlib
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import figstyle as fs
fs.setup()
with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
    import case_study as cs

R, F = "../results", "../figures"
TVAR = ["kd490", "meansummermIST", "maxsummermIST", "CVsummermIST"]
REG_ORDER = ["JBMP", "MMP", "CSMC", "SIMP", "NCMP"]
REG_LAB = {"JBMP": "Jurien Bay", "MMP": "Marmion", "CSMC": "Cockburn Sound",
           "SIMP": "Shoalwater", "NCMP": "Ngari Capes"}
REG_COL = {"JBMP": "#0072B2", "MMP": "#E69F00", "CSMC": "#009E73",
           "SIMP": "#CC79A7", "NCMP": "#5D3A9B"}
RF_KW = dict(n_estimators=600, min_samples_leaf=3, max_features=0.6, random_state=7, n_jobs=1)

agg = cs.agg.copy().reset_index(drop=True)
grid = pd.read_csv(f"{R}/case_map_grid.csv")
YEARS = np.arange(int(agg["year"].min()), int(agg["year"].max()) + 1)


def idw(site_lon, site_lat, vals, qlon, qlat, k=6, p=2):
    out = np.zeros(len(qlon))
    for i in range(len(qlon)):
        d = cs.haversine(qlon[i], qlat[i], site_lon, site_lat)
        j = np.argsort(d)[:k]; w = 1.0 / (d[j] ** p + 1e-6)
        out[i] = np.sum(w * vals[j]) / np.sum(w)
    return out


def predict_year_field():
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.preprocessing import StandardScaler
    from scipy.spatial.distance import cdist
    y = np.log1p(agg["dens"].to_numpy())
    Xtr = cs.build_design(agg, agg, "stRF")[0]
    rf = RandomForestRegressor(**RF_KW).fit(Xtr, y)
    raw = agg["dens"].to_numpy()
    Ltr = rf.apply(Xtr)
    leafmeans = [pd.Series(raw).groupby(Ltr[:, b]).mean() for b in range(Ltr.shape[1])]

    # ---- temporal AOA reference (Meyer & Pebesma 2021), matching the spatial version
    # but referenced to TRAINING SITE-YEARS with each year's covariates, so an unusual
    # environment-and-location combination in a given year is flagged as extrapolation.
    aoa_vars = cs.COVARS + ["longitude", "latitude"]
    ref = agg[aoa_vars].to_numpy()
    sc = StandardScaler().fit(ref)
    Zref = sc.transform(ref)
    imp = RandomForestRegressor(n_estimators=300, min_samples_leaf=3, random_state=1,
                                n_jobs=1).fit(ref, y).feature_importances_
    w = np.sqrt(imp + 1e-9)
    Zref_w = Zref * w
    # Threshold calibrated on BETWEEN-SITE nearest-neighbour distances. The manuscript's
    # spatial AOA referenced site means for the same reason: a site's repeated years share
    # identical coordinates and depth (the high-importance features), so within-site
    # neighbours are near-duplicates that would shrink the distance scale and over-flag
    # everything between sites. Excluding same-site pairs keeps the temporal covariate
    # variation while calibrating to the genuine between-location coverage.
    D = cdist(Zref_w, Zref_w)
    same = agg["site"].to_numpy()[:, None] == agg["site"].to_numpy()[None, :]
    D[same] = np.inf
    dref_nn = D.min(axis=1)
    dbar = dref_nn.mean()
    DI_thresh = float(np.percentile(dref_nn / dbar, 95))

    sites = agg[["site", "region", "latitude", "longitude"]].drop_duplicates("site")
    per = {s: agg[agg["site"] == s].set_index("year") for s in sites["site"]}
    stat = grid[["depth", "htyperank"]].to_numpy()
    glon, glat = grid["longitude"].to_numpy(), grid["latitude"].to_numpy()

    fields, di_fields = {}, {}
    for Y in YEARS:
        q = grid[["longitude", "latitude", "region"]].copy(); q["year"] = Y
        q["depth"] = stat[:, 0]; q["htyperank"] = stat[:, 1]
        for c in TVAR:
            sv = []
            for _, s in sites.iterrows():
                sy = per[s["site"]]
                sv.append(sy.loc[Y, c] if Y in sy.index
                          else sy[c].iloc[np.argmin(np.abs(sy.index - Y))])
            q[c] = idw(sites["longitude"].to_numpy(), sites["latitude"].to_numpy(),
                       np.array(sv), glon, glat)
        Xq = cs.build_design(agg, q, "stRF")[1]
        Lq = rf.apply(Xq)
        pred = np.zeros(len(q))
        for b in range(Lq.shape[1]):
            pred += leafmeans[b].reindex(Lq[:, b]).to_numpy()
        fields[Y] = pred / Lq.shape[1]
        # DI for this year's corridor cells
        Zq = sc.transform(q[aoa_vars].to_numpy()) * w
        di_fields[Y] = cdist(Zq, Zref_w).min(axis=1) / dbar
    return fields, di_fields, DI_thresh


def figT1():
    fields, di_fields, DI_thresh = predict_year_field()
    # collapse corridor to latitude bands (density and DI)
    nb = 70
    edges = np.linspace(grid["latitude"].min(), grid["latitude"].max(), nb + 1)
    bidx = np.clip(np.digitize(grid["latitude"].to_numpy(), edges) - 1, 0, nb - 1)
    M = np.full((nb, len(YEARS)), np.nan)
    Mdi = np.full((nb, len(YEARS)), np.nan)
    for yi, Y in enumerate(YEARS):
        for b in range(nb):
            sel = bidx == b
            if sel.any():
                M[b, yi] = np.nanmean(fields[Y][sel])
                Mdi[b, yi] = np.nanmean(di_fields[Y][sel])
    frac_out = float(np.mean(np.concatenate([di_fields[Y] for Y in YEARS]) > DI_thresh) * 100)
    traj = pd.read_csv(f"{R}/temporal_regional_traj.csv")
    site = agg[["site", "region", "latitude"]].drop_duplicates("site")
    vmax = float(np.nanpercentile(agg["dens"], 98))
    ext = [YEARS[0] - 0.5, YEARS[-1] + 0.5, grid["latitude"].min(), grid["latitude"].max()]
    xt = np.arange(YEARS[0], YEARS[-1] + 1, 3)

    def park_axis(ax, right=True):
        ticks = [site[site["region"] == rg]["latitude"].mean() for rg in REG_ORDER]
        ax.set_yticks(ticks); ax.set_yticklabels([REG_LAB[rg] for rg in REG_ORDER], fontsize=7.5)
        if right:
            sec = ax.secondary_yaxis("right")
            latt = np.arange(-33.0, -30.0, 1.0)
            sec.set_yticks(latt); sec.set_yticklabels([f"{abs(v):.0f}\u00b0S" for v in latt], fontsize=7)

    # per-latitude-band depth profile from the corridor cells (mean and cross-width range)
    depth_abs = -grid["depth"].to_numpy()
    yb = 0.5 * (edges[:-1] + edges[1:])
    dep_mean = np.array([np.nanmean(depth_abs[bidx == b]) if (bidx == b).any() else np.nan for b in range(nb)])
    dep_lo = np.array([np.nanmin(depth_abs[bidx == b]) if (bidx == b).any() else np.nan for b in range(nb)])
    dep_hi = np.array([np.nanmax(depth_abs[bidx == b]) if (bidx == b).any() else np.nan for b in range(nb)])
    dmax = float(np.ceil(np.nanmax(dep_hi)) + 1)

    fig = plt.figure(figsize=(9.4, 10.6))
    gs = GridSpec(3, 2, figure=fig, width_ratios=[0.20, 1.0],
                  height_ratios=[0.62, 1.0, 1.0], hspace=0.30, wspace=0.05,
                  left=0.10, right=0.99, top=0.965, bottom=0.06)

    def depth_panel(ax):
        ax.fill_betweenx(yb, dep_lo, dep_hi, color="#9a9a9a", alpha=0.30, lw=0)
        ax.plot(dep_mean, yb, color="#2b2b2b", lw=1.3)
        ax.set_xlim(dmax, 0); ax.set_ylim(ext[2], ext[3])       # 0 m at right (meets field)
        ax.set_xticks([0, 5, 10, 15]); ax.tick_params(labelsize=6.8)
        ax.set_yticks([site[site["region"] == rg]["latitude"].mean() for rg in REG_ORDER])
        ax.set_yticklabels([REG_LAB[rg] for rg in REG_ORDER], fontsize=7.5)
        ax.grid(axis="x", color="0.88", lw=0.5)
        ax.set_xlabel("Depth (m)", fontsize=8)

    def right_lat(ax):
        sec = ax.secondary_yaxis("right")
        latt = np.arange(-33.0, -30.0, 1.0)
        sec.set_yticks(latt); sec.set_yticklabels([f"{abs(v):.0f}\u00b0S" for v in latt], fontsize=7)

    # (a) regional trajectories (spans both columns)
    axa = fig.add_subplot(gs[0, :])
    for rg in REG_ORDER:
        d = traj[traj["region"] == rg].sort_values("year")
        axa.plot(d["year"], d["obs"], "o", color=REG_COL[rg], ms=3.5, alpha=0.8)
        axa.plot(d["year"], d["pred"], "-", color=REG_COL[rg], lw=1.8)
    axa.set_ylabel("Mean density 0.04 m$^{-2}$"); axa.set_xlim(YEARS[0] - 0.5, YEARS[-1] + 0.5)
    axa.set_xticks(xt)
    axa.set_title("(a)  Regional mean density: observed (points) vs stRF (lines)",
                  fontsize=9.5, loc="left")
    axa.legend(handles=[Line2D([], [], color=REG_COL[r], lw=2, marker="o", ms=3.5,
               label=REG_LAB[r]) for r in REG_ORDER], fontsize=7.2, ncol=5,
               loc="upper center", bbox_to_anchor=(0.5, -0.16), frameon=False)

    # (b) latitude x time predicted density, with depth profile at left
    axb = fig.add_subplot(gs[1, 1])
    im = axb.imshow(M, origin="lower", extent=ext, aspect="auto", cmap="YlGnBu",
                    vmin=0, vmax=vmax, interpolation="bilinear")
    axb.scatter(agg["year"], agg["latitude"], c=agg["dens"], cmap="YlGnBu", vmin=0, vmax=vmax,
                s=14, edgecolor="black", linewidth=0.35, zorder=3)
    axb.set_xticks(xt); axb.tick_params(labelleft=False)
    axb.set_title("(b)  stRF predicted density along the coast through time "
                  "(circles = observed, same scale)", fontsize=9.5, loc="left")
    right_lat(axb)
    cbb = fig.colorbar(im, ax=axb, fraction=0.045, pad=0.02)
    cbb.set_label("Shoots 0.04 m$^{-2}$", fontsize=8); cbb.ax.tick_params(labelsize=7)
    depth_panel(fig.add_subplot(gs[1, 0], sharey=axb))

    # (c) latitude x time area of applicability, with depth profile at left
    axc = fig.add_subplot(gs[2, 1])
    imc = axc.imshow(Mdi, origin="lower", extent=ext, aspect="auto", cmap="YlOrRd",
                     vmin=0, vmax=float(np.nanpercentile(Mdi, 98)), interpolation="bilinear")
    axc.contour(YEARS, yb, Mdi, levels=[DI_thresh], colors="black", linewidths=1.2)
    axc.scatter(agg["year"], agg["latitude"], s=5, color="black", alpha=0.5, zorder=3)
    axc.set_xlabel("Year"); axc.set_xticks(xt); axc.tick_params(labelleft=False)
    axc.set_title(f"(c)  Area of applicability: dissimilarity index "
                  f"(inside black contour = supported; {frac_out:.0f}% outside)",
                  fontsize=9.5, loc="left")
    right_lat(axc)
    cbc = fig.colorbar(imc, ax=axc, fraction=0.045, pad=0.02)
    cbc.set_label("Dissimilarity index", fontsize=8); cbc.ax.tick_params(labelsize=7)
    cbc.ax.axhline(DI_thresh, color="black", lw=1.2)
    depth_panel(fig.add_subplot(gs[2, 0], sharey=axc))

    fs.save(fig, f"{F}/FigureT1_spatiotemporal_reconstruction.png")
    print(f"AOA: DI threshold {DI_thresh:.2f}; {frac_out:.1f}% of corridor-years outside AOA")


if __name__ == "__main__":
    figT1()
    print("latitude-time Hovmoller done")
