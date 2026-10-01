import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
R = os.environ.get("MS_RESULTS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")) + "/"
FIG = os.environ.get("MS_FIGURES", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")) + "/"

"""Case-study figures from the R outputs (manuscript style: Okabe-Ito, no embedded titles, shared legends)."""
import json, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon
from matplotlib.tri import Triangulation

try:
    import figstyle as fs; fs.setup()
except Exception: fs = None
O = R + "case_study/"; F = FIG + "main/"; FS = FIG + "si/"; os.makedirs(F, exist_ok=True); os.makedirs(FS, exist_ok=True)
MC = {"GLM": "#999999", "GAM": "#009E73", "RF": "#E69F00", "SpatialRF": "#56B4E9", "stRF": "#D55E00"}
ML = {"GLM": "GLM", "GAM": "GAM", "RF": "Standard RF", "SpatialRF": "Spatial RF", "stRF": "stRF"}
REG = ["JBMP", "MMP", "CSMC", "SIMP", "NCMP"]; RL = {"JBMP": "Jurien Bay", "MMP": "Marmion", "CSMC": "Cockburn Sound", "SIMP": "Shoalwater", "NCMP": "Ngari Capes"}
RC = {"JBMP": "#0072B2", "MMP": "#E69F00", "CSMC": "#009E73", "SIMP": "#CC79A7", "NCMP": "#5D3A9B"}
def save(fig, name):
    F_ = FS if name.startswith("FigureS") else F
    fig.savefig(F_ + name, dpi=300); plt.close(fig); print(name)

agg = pd.read_csv(O + "agg.csv"); site = pd.read_csv(O + "sites.csv"); grid = pd.read_csv(O + "corridor.csv")
cv = pd.read_csv(O + "cv_regimes.csv"); imp = pd.read_csv(O + "importance.csv"); dc = pd.read_csv(O + "dropcol.csv")
coast = json.load(open(os.environ.get("CASE_COAST", R + "../data/coast_wa.json")))

# ------------------------------------------------------------------ Figure 6
fig, ax = plt.subplots(1, 3, figsize=(11, 3.9)); fig.subplots_adjust(left=0.06, right=0.99, top=0.95, bottom=0.30, wspace=0.32)
regimes = [("random 10-fold", "random\n10-fold"), ("site-grouped", "site\ngrouped"), ("LORO", "leave one\nregion out")]
xs = np.arange(3); w = 0.16; meths = ["GLM", "GAM", "RF", "SpatialRF", "stRF"]
for i, m in enumerate(meths):
    v = []
    for r, _ in regimes:
        row = cv[(cv.regime == r) & (cv.method == m)]
        v.append(row.Spearman_within_park.iloc[0] if r == "LORO" and pd.notna(row.Spearman_within_park.iloc[0]) else row.Spearman_pooled.iloc[0])
    ax[0].bar(xs + (i - 2) * w, v, w, color=MC[m])
ax[0].set_xticks(xs); ax[0].set_xticklabels([l for _, l in regimes], fontsize=8); ax[0].set_ylabel("Spearman rank correlation"); ax[0].axhline(0, color="0.5", lw=0.6)
ax[0].text(0.02, 0.97, "(a)", transform=ax[0].transAxes, va="top", fontweight="bold")
# (b) permutation importance, grouped
grp = {"longitude": ["COORD1"], "latitude": ["COORD2"], "spatial lags": ["LAG1", "LAG2", "LAG3"], "distance fields": [f"EDF{i}" for i in range(1, 9)],
       "temporal": ["TEMP1"], "depth": ["depth"], "light (kd490)": ["kd490"], "habitat rank": ["htyperank"], "thermal (3)": ["meansummermIST", "maxsummermIST", "CVsummermIST"]}
sh = {k: imp[imp.feature.isin(v)].share.sum() for k, v in grp.items()}; order = sorted(sh, key=lambda k: -sh[k])
SC = {k: ("#0072B2" if k in ("longitude", "latitude", "spatial lags", "distance fields") else "#D55E00" if k == "temporal" else "#999999") for k in grp}
ax[1].barh(range(len(order)), [100 * sh[k] for k in order], color=[SC[k] for k in order]); ax[1].set_yticks(range(len(order))); ax[1].set_yticklabels(order, fontsize=8); ax[1].invert_yaxis()
ax[1].set_xlabel("Permutation importance (% of total)"); ax[1].text(0.98, 0.03, "(b)", transform=ax[1].transAxes, ha="right", fontweight="bold")
# (c) drop-column unique contribution
lab = {"longitude": "longitude", "latitude": "latitude", "both_coords": "both coordinates", "spatial_block": "all spatial features", "depth": "depth", "kd490": "light (kd490)", "htyperank": "habitat rank", "meansummermIST": "mean summer SST", "maxsummermIST": "max summer SST", "CVsummermIST": "CV summer SST", "temporal": "temporal"}
d2 = dc[dc.removed != "full"].sort_values("loss", ascending=False)
ax[2].barh(range(len(d2)), d2.loss, color=["#0072B2" if r in ("longitude", "latitude", "both_coords", "spatial_block") else "#D55E00" if r == "temporal" else "#999999" for r in d2.removed])
ax[2].set_yticks(range(len(d2))); ax[2].set_yticklabels([lab[r] for r in d2.removed], fontsize=8); ax[2].invert_yaxis(); ax[2].axvline(0, color="0.5", lw=0.6)
ax[2].set_xlabel("Loss in CV rank correlation when removed"); ax[2].text(0.98, 0.03, "(c)", transform=ax[2].transAxes, ha="right", fontweight="bold")
fig.legend(handles=[Patch(fc=MC[m], label=ML[m]) for m in meths] + [Patch(fc="#0072B2", label="spatial feature (b, c)"), Patch(fc="#D55E00", label="temporal (b, c)"), Patch(fc="#999999", label="environmental (b, c)")],
           loc="lower center", ncol=8, frameon=False, fontsize=7.6, bbox_to_anchor=(0.5, 0.0))
save(fig, "Figure6.png")

# ------------------------------------------------------------------ Figure 7
def draw_coast(ax):
    for poly in coast:
        P = np.array(poly)
        if P.shape[0] < 3: continue
        ax.add_patch(Polygon(P, closed=True, fc="#e9e4d8", ec="0.55", lw=0.4, zorder=1))
lon0, lon1 = 114.75, 116.05; lat0, lat1 = agg.latitude.min() - 0.12, agg.latitude.max() + 0.12
fig, axes = plt.subplots(1, 3, figsize=(7.6, 7.4), sharey=True); fig.subplots_adjust(wspace=0.10, top=0.92, bottom=0.09, left=0.10, right=0.98)
for a in axes: draw_coast(a); a.set_xlim(lon0, lon1); a.set_ylim(lat0, lat1); a.set_facecolor("white"); a.set_xlabel("Longitude (°E)"); a.tick_params(labelsize=7.5)
axes[0].set_ylabel("Latitude (°S)"); axes[0].set_yticks(np.arange(-33.5, -30, 0.5)); axes[0].set_yticklabels([f"{abs(v):.1f}" for v in np.arange(-33.5, -30, 0.5)])
for rg in REG:
    s = site[site.region == rg]; axes[0].scatter(s.longitude, s.latitude, s=8 + 1.6 * s.dens, color=RC[rg], edgecolor="black", lw=0.4, zorder=3)
    axes[0].text(s.longitude.mean() - 0.42, s.latitude.mean(), rg, fontsize=7.5, fontweight="bold", va="center")
axes[0].text(0.02, 0.98, "(a)", transform=axes[0].transAxes, va="top", fontweight="bold")
vlo, vhi = float(np.floor(min(site.dens.min(), grid.pred.min()))), float(np.ceil(max(site.dens.quantile(0.98), grid.pred.max())))
sc = axes[1].scatter(grid.longitude, grid.latitude, c=grid.pred, cmap="YlGnBu", s=9, marker="s", edgecolors="none", vmin=vlo, vmax=vhi, zorder=2)
axes[1].scatter(site.longitude, site.latitude, c=site.dens, cmap="YlGnBu", vmin=vlo, vmax=vhi, s=26, edgecolor="black", lw=0.5, zorder=4)
tri = Triangulation(grid.longitude.values, grid.latitude.values); T = tri.triangles; gx, gy = grid.longitude.values, grid.latitude.values
sx = np.cos(np.radians(gy.mean()))
def el(a, b): return np.sqrt((sx * (gx[a] - gx[b])) ** 2 + (gy[a] - gy[b]) ** 2)
emax = np.maximum.reduce([el(T[:, 0], T[:, 1]), el(T[:, 1], T[:, 2]), el(T[:, 0], T[:, 2])]); tri.set_mask(emax > np.percentile(emax, 92))
iso = axes[1].tricontour(tri, grid.depth.values, levels=[-15, -10, -5], colors="0.15", linewidths=0.6, alpha=0.85, zorder=3)
axes[1].clabel(iso, fmt={-5: "5 m", -10: "10 m", -15: "15 m"}, fontsize=5.2, inline=True)
cax = axes[1].inset_axes([0.84, 0.58, 0.05, 0.25]); cb = fig.colorbar(sc, cax=cax); cb.ax.yaxis.set_ticks_position("left"); cb.ax.tick_params(labelsize=6, pad=1.5, length=2); cb.ax.set_title("Shoots\n0.04 m$^{-2}$", fontsize=6.2, pad=3)
axes[1].legend(handles=[Line2D([], [], marker="o", ls="", mfc="0.75", mec="black", ms=6, label="Observed site"), Line2D([], [], color="0.15", lw=0.8, label="Depth (m)")], loc="upper right", fontsize=6.6, framealpha=0.92)
axes[1].text(0.02, 0.98, "(b)", transform=axes[1].transAxes, va="top", fontweight="bold")
axes[2].scatter(grid.longitude[grid.inside], grid.latitude[grid.inside], color="#56B4E9", s=9, marker="s", edgecolors="none", zorder=2)
axes[2].scatter(grid.longitude[~grid.inside], grid.latitude[~grid.inside], color="#D55E00", s=9, marker="s", edgecolors="none", zorder=2)
axes[2].scatter(site.longitude, site.latitude, s=7, color="black", zorder=4)
axes[2].legend(handles=[Patch(fc="#56B4E9", label="Inside AOA"), Patch(fc="#D55E00", label="Outside AOA (extrapolation)")], loc="upper right", fontsize=6.6, framealpha=0.92)
axes[2].text(0.02, 0.98, "(c)", transform=axes[2].transAxes, va="top", fontweight="bold")
save(fig, "Figure7.png")

# ------------------------------------------------------------------ Figure 8
hv = pd.read_csv(O + "hovmoller.csv"); traj = pd.read_csv(O + "traj.csv"); dp = pd.read_csv(O + "depth_profile.csv"); meta = json.load(open(O + "meta.json"))
years = sorted(hv.year.unique()); lats = sorted(hv.lat.unique()); M = hv.pivot(index="lat", columns="year", values="pred").loc[lats, years].values; Mdi = hv.pivot(index="lat", columns="year", values="DI").loc[lats, years].values
thr = meta["temporal_DI_thresh"]; ext = [years[0] - 0.5, years[-1] + 0.5, min(lats), max(lats)]; xt = np.arange(years[0], years[-1] + 1, 3); vmax = float(np.nanpercentile(agg.dens, 98))
fig = plt.figure(figsize=(9.4, 10.6)); gs = GridSpec(3, 2, figure=fig, width_ratios=[0.2, 1], height_ratios=[0.62, 1, 1], hspace=0.30, wspace=0.05, left=0.10, right=0.99, top=0.965, bottom=0.06)
axa = fig.add_subplot(gs[0, :])
for rg in REG:
    d = traj[traj.region == rg].sort_values("year"); axa.plot(d.year, d.obs, "o", color=RC[rg], ms=3.5, alpha=0.8); axa.plot(d.year, d.pred, "-", color=RC[rg], lw=1.8)
axa.set_ylabel("Mean density 0.04 m$^{-2}$"); axa.set_xlim(years[0] - 0.5, years[-1] + 0.5); axa.set_xticks(xt); axa.text(0.01, 0.97, "(a)", transform=axa.transAxes, va="top", fontweight="bold")
axa.legend(handles=[Line2D([], [], color=RC[r], lw=2, marker="o", ms=3.5, label=RL[r]) for r in REG], fontsize=7.2, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.16), frameon=False)
def depth_panel(ax):
    ax.fill_betweenx(dp.lat, dp.dep_lo, dp.dep_hi, color="#9a9a9a", alpha=0.3, lw=0); ax.plot(dp.dep_mean, dp.lat, color="#2b2b2b", lw=1.3)
    ax.set_xlim(np.ceil(dp.dep_hi.max()) + 1, 0); ax.set_ylim(ext[2], ext[3]); ax.set_xticks([0, 5, 10, 15]); ax.tick_params(labelsize=6.8)
    ax.set_yticks([site[site.region == rg].latitude.mean() for rg in REG]); ax.set_yticklabels([RL[rg] for rg in REG], fontsize=7.5); ax.grid(axis="x", color="0.88", lw=0.5); ax.set_xlabel("Depth (m)", fontsize=8)
def right_lat(ax):
    sec = ax.secondary_yaxis("right"); latt = np.arange(-33.0, -30.0, 1.0); sec.set_yticks(latt); sec.set_yticklabels([f"{abs(v):.0f}°S" for v in latt], fontsize=7)
axb = fig.add_subplot(gs[1, 1]); im = axb.imshow(M, origin="lower", extent=ext, aspect="auto", cmap="YlGnBu", vmin=0, vmax=vmax, interpolation="bilinear")
axb.scatter(agg.year, agg.latitude, c=agg.dens, cmap="YlGnBu", vmin=0, vmax=vmax, s=14, edgecolor="black", linewidth=0.35, zorder=3); axb.set_xticks(xt); axb.tick_params(labelleft=False); axb.text(0.01, 0.97, "(b)", transform=axb.transAxes, va="top", fontweight="bold"); right_lat(axb)
cbb = fig.colorbar(im, ax=axb, fraction=0.045, pad=0.02); cbb.set_label("Shoots 0.04 m$^{-2}$", fontsize=8); cbb.ax.tick_params(labelsize=7); depth_panel(fig.add_subplot(gs[1, 0], sharey=axb))
axc = fig.add_subplot(gs[2, 1]); imc = axc.imshow(Mdi, origin="lower", extent=ext, aspect="auto", cmap="YlOrRd", vmin=0, vmax=float(np.nanpercentile(Mdi, 98)), interpolation="bilinear")
axc.contour(years, lats, Mdi, levels=[thr], colors="black", linewidths=1.2); axc.scatter(agg.year, agg.latitude, s=5, color="black", alpha=0.5, zorder=3)
axc.set_xlabel("Year"); axc.set_xticks(xt); axc.tick_params(labelleft=False); axc.text(0.01, 0.97, "(c)", transform=axc.transAxes, va="top", fontweight="bold"); right_lat(axc)
cbc = fig.colorbar(imc, ax=axc, fraction=0.045, pad=0.02); cbc.set_label("Dissimilarity index", fontsize=8); cbc.ax.tick_params(labelsize=7); cbc.ax.axhline(thr, color="black", lw=1.2); depth_panel(fig.add_subplot(gs[2, 0], sharey=axc))
save(fig, "Figure8.png")

# ------------------------------------------------------------------ Figure S5: latitudinal transect
g2 = grid.sort_values("latitude"); band = g2.groupby(pd.cut(g2.latitude, 80)).agg(lat=("latitude", "mean"), pred=("pred", "mean"), lo=("lo", "mean"), hi=("hi", "mean"), inside=("inside", "mean")).dropna()
fig, ax = plt.subplots(figsize=(8.6, 3.9)); fig.subplots_adjust(left=0.08, right=0.99, top=0.95, bottom=0.30)
ax.fill_between(band.lat, band.lo, band.hi, color="#D55E00", alpha=0.18, lw=0); ax.plot(band.lat, band.pred, color="#D55E00", lw=1.8)
out = band[band.inside < 0.5]
for _, r in out.iterrows(): ax.axvspan(r.lat - 0.022, r.lat + 0.022, color="0.85", lw=0, zorder=0)
for rg in REG:
    s = site[site.region == rg]; ax.scatter(s.latitude, s.dens, color=RC[rg], s=18, edgecolor="black", lw=0.4, zorder=3)
ax.set_xlabel("Latitude (°S)"); ax.set_ylabel("Shoot density 0.04 m$^{-2}$"); ax.set_xticks(np.arange(-33.5, -30, 0.5)); ax.set_xticklabels([f"{abs(v):.1f}" for v in np.arange(-33.5, -30, 0.5)])
fig.legend(handles=[Line2D([], [], color="#D55E00", lw=1.8, label="stRF prediction (corridor mean)"), Patch(fc="#D55E00", alpha=0.18, label="80% quantile forest band"), Patch(fc="0.85", label="outside area of applicability")] + [Line2D([], [], marker="o", ls="", mfc=RC[r], mec="black", ms=5, label=RL[r]) for r in REG], loc="lower center", ncol=4, frameon=False, fontsize=7.6, bbox_to_anchor=(0.5, 0.0))
save(fig, "FigureS5.png")

# ------------------------------------------------------------------ Figure S6: forecast skill
ro = pd.read_csv(O + "rolling_origin.csv"); lv = ro.pivot_table(index="horizon", columns="method", values="rho_level"); an = ro[ro.method == "stRF"].groupby("horizon").rho_anomaly.mean()
fig, ax = plt.subplots(1, 2, figsize=(9, 3.9)); fig.subplots_adjust(left=0.08, right=0.985, top=0.94, bottom=0.26, wspace=0.26)
FC = {"stRF": "#D55E00", "persistence": "#0072B2", "climatology": "#666666"}
for m in ["stRF", "persistence", "climatology"]: ax[0].plot(lv.index, lv[m], "-o", color=FC[m], ms=5, lw=1.8)
ax[0].set_xlabel("Forecast horizon (years ahead)"); ax[0].set_ylabel("Spearman (predicted vs observed)"); ax[0].set_xticks(lv.index); ax[0].text(0.02, 0.97, "(a)", transform=ax[0].transAxes, va="top", fontweight="bold")
ax[1].axhline(0, color="0.6", lw=1, ls="--"); ax[1].plot(an.index, an.values, "-o", color="#D55E00", ms=5, lw=1.8); ax[1].set_ylim(-0.35, 0.35); ax[1].set_xticks(an.index)
ax[1].set_xlabel("Forecast horizon (years ahead)"); ax[1].set_ylabel("Corr(predicted, observed) anomaly"); ax[1].text(0.02, 0.97, "(b)", transform=ax[1].transAxes, va="top", fontweight="bold")
fig.legend(handles=[Line2D([], [], color=FC[m], marker="o", ms=5, lw=1.8, label=m) for m in FC], loc="lower center", ncol=3, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, 0.0))
save(fig, "FigureS6.png")

# ------------------------------------------------------------------ Figure S7: interval coverage and width
iv = pd.read_csv(O + "intervals.csv"); s = iv.groupby(["regime", "interval"])[["coverage", "coverage_inside", "coverage_outside", "width"]].mean().reset_index()
IC = {"ensemble spread": "#999999", "QRF": "#0072B2", "split conformal": "#009E73"}
fig, ax = plt.subplots(1, 2, figsize=(9, 3.9)); fig.subplots_adjust(left=0.08, right=0.985, top=0.94, bottom=0.28, wspace=0.26)
regs = ["interpolation (random 5-fold)", "extrapolation (LORO)"]; xs = np.arange(2); w = 0.26
for i, itv in enumerate(IC):
    v = [s[(s.regime == r) & (s.interval == itv)].coverage.iloc[0] for r in regs]; ax[0].bar(xs + (i - 1) * w, v, w, color=IC[itv])
    wv = [s[(s.regime == r) & (s.interval == itv)].width.iloc[0] for r in regs]; ax[1].bar(xs + (i - 1) * w, wv, w, color=IC[itv])
ax[0].axhline(0.8, color="black", lw=1, ls="--"); ax[0].set_ylim(0, 1); ax[0].set_xticks(xs); ax[0].set_xticklabels(["interpolation\n(random 5-fold)", "extrapolation\n(leave one region out)"], fontsize=8); ax[0].set_ylabel("Empirical coverage of nominal 80% interval")
ax[0].text(0.02, 0.97, "(a)", transform=ax[0].transAxes, va="top", fontweight="bold")
ax[1].set_xticks(xs); ax[1].set_xticklabels(["interpolation\n(random 5-fold)", "extrapolation\n(leave one region out)"], fontsize=8); ax[1].set_ylabel("Mean interval width (shoots 0.04 m$^{-2}$)"); ax[1].text(0.02, 0.97, "(b)", transform=ax[1].transAxes, va="top", fontweight="bold")
fig.legend(handles=[Patch(fc=IC[k], label=k) for k in IC] + [Line2D([], [], color="black", ls="--", label="nominal 80%")], loc="lower center", ncol=4, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, 0.0))
save(fig, "FigureS7.png")
