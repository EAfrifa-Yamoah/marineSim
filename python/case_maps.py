"""
case_maps.py
============
(1) Drop-column ("unique contribution") importance for every stRF feature group,
    to set beside permutation importance and expose the spatial-feature redundancy
    behind longitude's high permutation ranking.
(2) A coastal-corridor prediction surface: a continuous ribbon following the WA
    coastline (a smooth longitude-vs-latitude spine, buffered cross-shelf), with
    covariates interpolated from the monitoring sites, stRF predictions, and the
    area-of-applicability dissimilarity index for every cell.
(3) Site-level means for an observed-density map.
"""
import io, contextlib, json
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, GroupKFold
from scipy.stats import spearmanr
from statsmodels.nonparametric.smoothers_lowess import lowess

with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
    import case_study as cs

RES = cs.RES
agg = cs.agg.copy()
COVARS = list(cs.COVARS)
ncov, nlag = len(COVARS), len(cs.RADII_KM)

# ============================================== (1) drop-column importance
Xall, _ = cs.build_design(agg, agg, "stRF")
y = np.log1p(agg["dens"].to_numpy())
nedf = Xall.shape[1] - ncov - nlag - 3
gi = {
    "Longitude": [Xall.shape[1] - 2], "Latitude": [Xall.shape[1] - 1],
    "Depth": [0], "Light (kd490)": [1], "Habitat rank": [2],
    "Summer SST": [3, 4, 5], "Spatial lags": list(range(ncov, ncov + nlag)),
    "Distance fields": list(range(ncov + nlag, ncov + nlag + nedf)),
    "Temporal lag": [ncov + nlag + nedf],
}

def cv_skill(drop):
    keep = [c for c in range(Xall.shape[1]) if c not in drop]
    X = Xall[:, keep]; kf = KFold(10, shuffle=True, random_state=11); yt, yp = [], []
    for tr, te in kf.split(X):
        rf = RandomForestRegressor(n_estimators=400, min_samples_leaf=3,
                                   max_features=0.6, random_state=7, n_jobs=1).fit(X[tr], y[tr])
        yp.append(rf.predict(X[te])); yt.append(y[te])
    return spearmanr(np.concatenate(yt), np.concatenate(yp)).correlation

full = cv_skill([])
perm = pd.read_csv(f"{RES}/case_importance.csv").set_index("feature")["importance_pct"]
rows = []
for g, cols in gi.items():
    loss = full - cv_skill(cols)
    rows.append(dict(feature=g, perm_pct=float(perm.get(g, np.nan)), unique_loss=float(loss)))
dc = pd.DataFrame(rows)
dc.to_csv(f"{RES}/case_dropcol.csv", index=False)
allspatial = full - cv_skill(gi["Longitude"] + gi["Latitude"] + gi["Spatial lags"] + gi["Distance fields"])

# ============================================== (3) site means
aggspec = {c: (c, "mean") for c in (["longitude", "latitude", "dens"] + COVARS)}
aggspec["region"] = ("region", "first")
site = agg.groupby("site").agg(**aggspec).reset_index()
site.to_csv(f"{RES}/case_sites.csv", index=False)

# ============================================== (2) coastal corridor grid
# coastal spine: smooth longitude as a function of latitude (captures the E-W bulge)
o = np.argsort(site["latitude"].to_numpy())
sp = lowess(site["longitude"].to_numpy()[o], site["latitude"].to_numpy()[o], frac=0.5, return_sorted=True)
def spine_lon(lat):
    return np.interp(lat, sp[:, 0], sp[:, 1])

lat_min, lat_max = site["latitude"].min() - 0.05, site["latitude"].max() + 0.05
lats = np.arange(lat_min, lat_max, 0.022)
lon_buf = 0.24   # ~22 km cross-shelf half-width
cells = []
for la in lats:
    c = spine_lon(la)
    for lo in np.arange(c - lon_buf, c + lon_buf + 1e-9, 0.022):
        cells.append((lo, la))
grid = pd.DataFrame(cells, columns=["longitude", "latitude"])

# ---- clip the corridor to ocean: drop cells that fall on land ----
# uses the same geo-maps coastline polygons rendered in Figure 7
from matplotlib.path import Path as MplPath
coast = json.load(open(f"{RES}/coast_wa.json"))
pts = grid[["longitude", "latitude"]].to_numpy()
on_land = np.zeros(len(grid), dtype=bool)
for ring in coast:
    ring = np.asarray(ring)
    if ring.shape[0] < 3:
        continue
    if (pts[:, 0].max() < ring[:, 0].min() or pts[:, 0].min() > ring[:, 0].max() or
            pts[:, 1].max() < ring[:, 1].min() or pts[:, 1].min() > ring[:, 1].max()):
        continue                      # ring bbox misses the grid entirely
    on_land |= MplPath(ring).contains_points(pts)
grid = grid[~on_land].reset_index(drop=True)
print(f"corridor cells: {len(on_land)} built, {int(on_land.sum())} on land removed, {len(grid)} in ocean")

# nearest site -> region; representative year
sx = site[["longitude", "latitude"]].to_numpy()
def nearest_site(lo, la):
    d = cs.haversine(lo, la, sx[:, 0], sx[:, 1]); return int(np.argmin(d))
ns = [nearest_site(lo, la) for lo, la in zip(grid["longitude"], grid["latitude"])]
grid["region"] = site["region"].to_numpy()[ns]
grid["year"] = int(agg["year"].max())

# IDW interpolation of covariates from site means (k nearest, inverse-square).
# k=4 rather than 8: averaging 8 sites retained only ~82% of the between-site SD in
# depth, flattening the covariate contrasts that drive the high-density predictions.
# k=4 retains ~94% while still giving a smooth surface between sites.
def idw(lo, la, col, k=4, p=2):
    d = cs.haversine(lo, la, sx[:, 0], sx[:, 1])
    idx = np.argsort(d)[:k]; w = 1.0 / (d[idx] ** p + 1e-6)
    return float(np.sum(w * site[col].to_numpy()[idx]) / np.sum(w))
for cov in COVARS:
    grid[cov] = [idw(lo, la, cov) for lo, la in zip(grid["longitude"], grid["latitude"])]

# stRF prediction surface.
#
# The forest is fitted on log1p(density) (splits are more stable on that scale), but the
# MAP must report density on the raw scale. Back-transforming with expm1(mean log) is a
# median-like quantity, not a mean, and is biased low by ~14% here; it was the dominant
# cause of the underprediction visible in the sampled areas. A Duan (1983) smearing
# factor removes the average bias but is heteroscedastic across parks (0.97-1.28), so a
# single factor over-corrects the low-density park and under-corrects the high-variance
# one.
#
# Instead, read the conditional mean directly off the leaves on the RAW scale: each tree
# predicts the mean observed density of the training rows sharing the query's leaf, and
# the forest averages those. This is the quantile-regression-forest view of a forest as
# an adaptive nearest-neighbour weighting (Meinshausen 2006). It removes the
# retransformation bias by construction, needs no correction factor, and leaves the
# fitted model (and every rank-based result in the paper) untouched.
Xtr, Xgr = cs.build_design(agg, grid, "stRF")
rf = RandomForestRegressor(n_estimators=600, min_samples_leaf=3, max_features=0.6,
                           random_state=7, n_jobs=1).fit(Xtr, y)

dens_raw = agg["dens"].to_numpy()
L_tr, L_gr = rf.apply(Xtr), rf.apply(Xgr)
pred = np.zeros(len(Xgr))
for b in range(L_tr.shape[1]):
    leaf_mean = pd.Series(dens_raw).groupby(L_tr[:, b]).mean()
    pred += leaf_mean.reindex(L_gr[:, b]).to_numpy()
grid["pred"] = pred / L_tr.shape[1]

# diagnostic only: the smearing factor the naive back-transform would have needed
oof = np.zeros(len(y))
for tr_i, te_i in GroupKFold(n_splits=5).split(Xtr, y, groups=agg["site"]):
    oof[te_i] = RandomForestRegressor(
        n_estimators=600, min_samples_leaf=3, max_features=0.6,
        random_state=7, n_jobs=1).fit(Xtr[tr_i], y[tr_i]).predict(Xtr[te_i])
SMEAR = float(np.mean(np.exp(y - oof)))

# AOA dissimilarity index (replicates case_study.py method) in {COVARS+lon+lat} space
aoa_vars = COVARS + ["longitude", "latitude"]
site_ref = agg.groupby("site")[aoa_vars].mean()
sc = StandardScaler().fit(site_ref.to_numpy())
Zref = sc.transform(site_ref.to_numpy()); Zgr = sc.transform(grid[aoa_vars].to_numpy())
imp = RandomForestRegressor(n_estimators=300, min_samples_leaf=3, random_state=1,
                            n_jobs=1).fit(agg[aoa_vars].to_numpy(), y).feature_importances_
w = np.sqrt(imp + 1e-9); Zref_w = Zref * w; Zgr_w = Zgr * w
def nn_min(A, B, excl=False):
    out = np.zeros(len(A))
    for i in range(len(A)):
        dd = np.sqrt(((B - A[i]) ** 2).sum(1))
        if excl: dd[i] = np.inf
        out[i] = dd.min()
    return out
d_ref = nn_min(Zref_w, Zref_w, excl=True); dbar = d_ref.mean()
DI_thresh = float(np.percentile(d_ref / dbar, 95))
grid["DI"] = nn_min(Zgr_w, Zref_w) / dbar
grid["inside"] = (grid["DI"] <= DI_thresh).astype(int)
grid.to_csv(f"{RES}/case_map_grid.csv", index=False)

with open(f"{RES}/case_maps_meta.json", "w") as f:
    json.dump(dict(full_cv=full, all_spatial_loss=float(allspatial),
                   DI_thresh=DI_thresh, grid_cells=len(grid),
                   pct_inside=float(grid["inside"].mean() * 100),
                   estimator="leaf raw-scale conditional mean",
                   naive_smearing_factor_would_have_been=SMEAR, idw_k=4,
                   pred_min=float(grid["pred"].min()),
                   pred_max=float(grid["pred"].max()),
                   pred_mean=float(grid["pred"].mean())), f, indent=2)

print(f"Full random-CV Spearman: {full:.3f}")
print("Drop-column unique contribution (Spearman loss) vs permutation reliance (%):")
for _, r in dc.sort_values("perm_pct", ascending=False).iterrows():
    print(f"  {r['feature']:16s} perm {r['perm_pct']:5.1f}%   unique loss {r['unique_loss']:+.3f}")
print(f"  {'ALL SPATIAL (group)':16s}            unique loss {allspatial:+.3f}")
print(f"\nCorridor grid: {len(grid)} cells; inside AOA {grid['inside'].mean()*100:.0f}%; DI thresh {DI_thresh:.2f}")
print("Wrote case_dropcol.csv, case_sites.csv, case_map_grid.csv, case_maps_meta.json")
