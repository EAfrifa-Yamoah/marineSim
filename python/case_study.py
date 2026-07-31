"""
case_study.py
=============
Posidonia sinuosa shoot density along the Western Australian latitudinal gradient.

Aim: use the pockets of (temporally uneven) monitoring data to extrapolate shoot
density across the full latitudinal extent, and quantify where that extrapolation
is and is not supported (area of applicability).

Two analyses:
  (1) Leave-One-Region-Out (LORO) CV: each marine park is held out in turn and its
      site-year mean shoot density predicted from the other four. This is a direct
      test of latitudinal transfer for GLM, GAM, RF, SpatialRF and stRF.
  (2) Full-extent prediction on a latitudinal transect with an area-of-applicability
      (AOA) overlay (Meyer & Pebesma 2021 dissimilarity index).
"""
import os, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

RES = os.path.join(os.path.dirname(__file__), "..", "results")
UP = "/mnt/user-data/uploads/Bayesiandataset_2025_final.csv"
rng = np.random.default_rng(424242)

FULLNAME = {"MMP": "Marmion Marine Park", "CSMC": "Cockburn Sound",
            "SIMP": "Shoalwater Islands Marine Park",
            "NCMP": "Ngari Capes Marine Park", "JBMP": "Jurien Bay Marine Park"}
COVARS = ["depth", "kd490", "htyperank", "meansummermIST", "maxsummermIST", "CVsummermIST"]
RADII_KM = [25, 75, 150]

# ----------------------------------------------------------------- load + aggregate
raw = pd.read_csv(UP)
agg = (raw.groupby(["Region_5", "Site", "YearSeagrassSampling"])
       .agg(dens=("rawcounts", "mean"),
            latitude=("latitude", "mean"), longitude=("longitude", "mean"),
            depth=("depth", "mean"), kd490=("kd490", "mean"),
            htyperank=("htyperank", "mean"),
            meansummermIST=("meansummermIST", "mean"),
            maxsummermIST=("maxsummermIST", "mean"),
            CVsummermIST=("CVsummermIST", "mean"))
       .reset_index().rename(columns={"YearSeagrassSampling": "year",
                                      "Region_5": "region", "Site": "site"}))
# impute residual covariate NAs with region-year then global means
for c in COVARS:
    agg[c] = agg.groupby("region")[c].transform(lambda s: s.fillna(s.mean()))
    agg[c] = agg[c].fillna(agg[c].mean())
agg = agg.dropna(subset=["dens"]).reset_index(drop=True)
print(f"site-year records: {len(agg)}  sites: {agg['site'].nunique()}  "
      f"years: {agg['year'].min()}-{agg['year'].max()}")

# great-circle distance (km) between two arrays of lon/lat
def haversine(lon1, lat1, lon2, lat2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1); dlmb = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2)**2
    return 2 * R * np.arcsin(np.sqrt(a))

# ----------------------------------------------------------------- stRF features
def spatial_lag_features(train, query):
    """IDW of training site-year densities at multiple radii, per query row,
    matched within the same year where possible (else all-year fallback)."""
    feats = np.zeros((len(query), len(RADII_KM)))
    tr_lon = train["longitude"].to_numpy(); tr_lat = train["latitude"].to_numpy()
    tr_y = train["year"].to_numpy(); tr_d = train["dens"].to_numpy()
    for i, (_, q) in enumerate(query.iterrows()):
        d = haversine(q["longitude"], q["latitude"], tr_lon, tr_lat)
        same = tr_y == q["year"]
        for j, r in enumerate(RADII_KM):
            for mask in (same, np.ones_like(same, dtype=bool)):
                sel = mask & (d <= r) & (d > 1e-6)
                if sel.sum() >= 2:
                    w = 1.0 / (d[sel] + 1.0)
                    feats[i, j] = np.sum(w * tr_d[sel]) / np.sum(w)
                    break
            else:
                # nearest-k fallback so the feature is always defined
                k = min(5, (d > 1e-6).sum())
                idx = np.argsort(d)[:k]
                feats[i, j] = tr_d[idx].mean()
    return feats

def edf_features(train_coords, query_coords, k=8, anchors=None):
    """Euclidean distance fields to k anchors placed by k-means on the training
    coordinates.

    The anchors MUST be shared between the training and the query design matrices.
    Drawing a fresh anchor set for the query (as an earlier version did) misaligns the
    columns: the forest learns splits on 'distance to anchor 3' for one set of anchors
    and is then asked to predict using distances to a different set, which turns the
    distance fields into noise at prediction time and shrinks predictions toward the
    mean. Pass the anchors returned by the training call back in for the query call.
    """
    if anchors is None:
        anchors = KMeans(n_clusters=min(k, len(train_coords)), n_init=5,
                         random_state=0).fit(train_coords).cluster_centers_
    def dmat(C):
        out = np.zeros((len(C), len(anchors)))
        for a in range(len(anchors)):
            out[:, a] = haversine(C[:, 0], C[:, 1], anchors[a, 0], anchors[a, 1])
        return out
    return dmat(query_coords), anchors

def temporal_feature(train, query):
    """previous-year regional mean density, broadcast to query rows."""
    reg_year = train.groupby(["region", "year"])["dens"].mean()
    vals = np.zeros(len(query))
    for i, (_, q) in enumerate(query.iterrows()):
        key = (q["region"], q["year"] - 1)
        if key in reg_year.index:
            vals[i] = reg_year[key]
        else:
            ry = reg_year[reg_year.index.get_level_values(0) == q["region"]]
            vals[i] = ry.mean() if len(ry) else train["dens"].mean()
    return vals

# ----------------------------------------------------------------- model fitting
def build_design(train, query, method):
    Xc_tr = train[COVARS].to_numpy(); Xc_q = query[COVARS].to_numpy()
    if method in ("GLM", "RF"):
        return Xc_tr, Xc_q
    if method == "SpatialRF":
        return (np.column_stack([Xc_tr, train[["longitude", "latitude"]].to_numpy()]),
                np.column_stack([Xc_q, query[["longitude", "latitude"]].to_numpy()]))
    if method == "stRF":
        lag_tr = spatial_lag_features(train, train)
        lag_q = spatial_lag_features(train, query)
        ed_tr, km = edf_features(train[["longitude", "latitude"]].to_numpy(),
                                 train[["longitude", "latitude"]].to_numpy())
        ed_q, _ = edf_features(train[["longitude", "latitude"]].to_numpy(),
                               query[["longitude", "latitude"]].to_numpy(), anchors=km)
        tmp_tr = temporal_feature(train, train)[:, None]
        tmp_q = temporal_feature(train, query)[:, None]
        Xtr = np.column_stack([Xc_tr, lag_tr, ed_tr, tmp_tr,
                               train[["longitude", "latitude"]].to_numpy()])
        Xq = np.column_stack([Xc_q, lag_q, ed_q, tmp_q,
                              query[["longitude", "latitude"]].to_numpy()])
        return Xtr, Xq

def fit_predict(train, query, method):
    ytr = np.log1p(train["dens"].to_numpy())
    if method == "GLM":
        Xtr, Xq = build_design(train, query, method)
        sc = StandardScaler().fit(Xtr)
        m = LinearRegression().fit(sc.transform(Xtr), ytr)
        pred = m.predict(sc.transform(Xq))
    elif method == "GAM":
        from pygam import LinearGAM, s, te
        import contextlib
        Xg_tr = np.column_stack([train[["longitude", "latitude"]].to_numpy(),
                                 train[COVARS].to_numpy()])
        Xg_q = np.column_stack([query[["longitude", "latitude"]].to_numpy(),
                                query[COVARS].to_numpy()])
        terms = te(0, 1) + s(2) + s(3) + s(5) + s(6)
        with open(os.devnull, "w") as dn, contextlib.redirect_stdout(dn), \
                contextlib.redirect_stderr(dn):
            gam = LinearGAM(terms, max_iter=60).fit(Xg_tr, ytr)
            pred = gam.predict(Xg_q)
        # spatial tensors extrapolate explosively beyond the sampled range;
        # constrain to a biologically sane band on the log scale
        pred = np.clip(pred, 0.0, np.log1p(150.0))
    else:
        Xtr, Xq = build_design(train, query, method)
        m = RandomForestRegressor(n_estimators=400, min_samples_leaf=3,
                                  max_features=0.6, random_state=7, n_jobs=1)
        m.fit(Xtr, ytr)
        pred = m.predict(Xq)
    return np.expm1(pred)

def metrics(y, yhat):
    y = np.asarray(y); yhat = np.asarray(yhat)
    ss_res = np.sum((y - yhat)**2); ss_tot = np.sum((y - y.mean())**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    rho = spearmanr(y, yhat).correlation
    rmse = np.sqrt(np.mean((y - yhat)**2))
    return r2, rho, rmse

# ----------------------------------------------------------------- (1) LORO CV
METHODS = ["GLM", "GAM", "RF", "SpatialRF", "stRF"]
regions = ["JBMP", "MMP", "CSMC", "SIMP", "NCMP"]   # north -> south
loro_rows = []
preds_store = {}
for held in regions:
    tr = agg[agg.region != held].copy()
    te = agg[agg.region == held].copy()
    for meth in METHODS:
        try:
            yhat = fit_predict(tr, te, meth)
            r2, rho, rmse = metrics(te["dens"].to_numpy(), yhat)
        except Exception as e:
            r2 = rho = rmse = np.nan
            yhat = np.full(len(te), np.nan)
        loro_rows.append(dict(held=held, held_name=FULLNAME[held], method=meth,
                              n_test=len(te), R2=r2, Spearman=rho, RMSE=rmse,
                              lat=te["latitude"].mean()))
        preds_store[(held, meth)] = (te["dens"].to_numpy(), yhat,
                                     te["latitude"].to_numpy())
loro = pd.DataFrame(loro_rows)
loro.to_csv(f"{RES}/case_loro.csv", index=False)

pivot_r2 = loro.pivot(index="held", columns="method", values="R2").reindex(regions)
pivot_rho = loro.pivot(index="held", columns="method", values="Spearman").reindex(regions)
print("\n=== LORO R^2 by held-out region ===")
print(pivot_r2[METHODS].round(3).to_string())
print("\n=== LORO Spearman by held-out region ===")
print(pivot_rho[METHODS].round(3).to_string())
print("\nPooled mean (across regions):")
print(loro.groupby("method")[["R2", "Spearman", "RMSE"]].mean().reindex(METHODS).round(3).to_string())

# pooled-prediction R2/Spearman (concatenate all held-out predictions)
pooled = {}
for meth in METHODS:
    ys = np.concatenate([preds_store[(h, meth)][0] for h in regions])
    yh = np.concatenate([preds_store[(h, meth)][1] for h in regions])
    ok = np.isfinite(yh)
    pooled[meth] = dict(zip(["R2", "Spearman", "RMSE"], metrics(ys[ok], yh[ok])))
print("\nPooled-across-all-held-out predictions:")
for m in METHODS:
    print(f"  {m:10s} R2={pooled[m]['R2']:+.3f}  rho={pooled[m]['Spearman']:+.3f}  RMSE={pooled[m]['RMSE']:.2f}")

# ----------------------------------------------------------------- (2) full-extent + AOA
# train stRF on all data; predict on a latitudinal transect
full = agg.copy()
# smooth lon(lat) trend along coast
lon_fit = np.polyfit(full["latitude"], full["longitude"], 2)
lat_grid = np.linspace(full["latitude"].min() - 0.15, full["latitude"].max() + 0.15, 160)
lon_grid = np.polyval(lon_fit, lat_grid)
# representative environment: median covariates, latest-year thermal climatology
rep = {c: float(full[c].median()) for c in COVARS}
grid = pd.DataFrame(dict(region="grid", site="g", year=full["year"].max(),
                         latitude=lat_grid, longitude=lon_grid,
                         **{c: rep[c] for c in COVARS}))
# let depth/kd490/thermal vary smoothly with latitude (local means)
for c in ["depth", "kd490", "meansummermIST", "maxsummermIST", "CVsummermIST", "htyperank"]:
    coef = np.polyfit(full["latitude"], full[c], 2)
    grid[c] = np.polyval(coef, lat_grid)

grid_pred = fit_predict(full, grid, "stRF")

# AOA: RF-importance-weighted dissimilarity index (Meyer & Pebesma 2021),
# computed in the {environmental covariates + coordinates} predictor space and
# referenced to UNIQUE SITES (site-mean predictors) so that repeated site-year
# sampling does not collapse the within-training distance scale.
aoa_vars = COVARS + ["longitude", "latitude"]
site_ref = full.groupby("site")[aoa_vars].mean()
Xref = site_ref.to_numpy()
Xgr = grid[aoa_vars].to_numpy()
sc = StandardScaler().fit(Xref)
Zref = sc.transform(Xref); Zgr = sc.transform(Xgr)
rf_imp = RandomForestRegressor(n_estimators=300, min_samples_leaf=3,
                               random_state=1, n_jobs=1).fit(
    full[aoa_vars].to_numpy(), np.log1p(full["dens"].to_numpy())).feature_importances_
w = np.sqrt(rf_imp + 1e-9)
Zref_w = Zref * w; Zgr_w = Zgr * w

def nn_min(A, B, exclude_self=False):
    out = np.zeros(len(A))
    for i in range(len(A)):
        dd = np.sqrt(np.sum((B - A[i])**2, axis=1))
        if exclude_self:
            dd[i] = np.inf
        out[i] = dd.min()
    return out

d_ref = nn_min(Zref_w, Zref_w, exclude_self=True)
dbar = d_ref.mean()
DI_ref = d_ref / dbar
DI_grid = nn_min(Zgr_w, Zref_w) / dbar
# outlier-removed 95th-percentile threshold (Meyer & Pebesma)
q75, q25 = np.percentile(DI_ref, [75, 25])
DI_thresh = float(np.percentile(DI_ref, 95))
inside_aoa = DI_grid <= DI_thresh

casegrid = pd.DataFrame(dict(latitude=lat_grid, longitude=lon_grid,
                             pred_density=grid_pred, DI=DI_grid,
                             inside_aoa=inside_aoa.astype(int)))
casegrid.to_csv(f"{RES}/case_grid.csv", index=False)

# region centroids / monitoring footprint for the figure
regsum = (agg.groupby("region").agg(lat=("latitude", "mean"),
          dens=("dens", "mean"), nyears=("year", "nunique"),
          nsite=("site", "nunique")).reindex(regions).reset_index())
regsum["name"] = regsum["region"].map(FULLNAME)
regsum.to_csv(f"{RES}/case_regions.csv", index=False)

frac_inside = float(inside_aoa.mean() * 100)
print(f"\nFull-extent transect: {len(lat_grid)} points; inside AOA = {frac_inside:.0f}%")
print(f"DI threshold = {DI_thresh:.3f}; gaps flagged where DI exceeds it.")

# ----------------------------------------------------------------- save derived
case_out = dict(
    n_site_years=int(len(agg)), n_sites=int(agg["site"].nunique()),
    year_min=int(agg["year"].min()), year_max=int(agg["year"].max()),
    loro_pooled_mean={m: {k: float(loro[loro.method == m][k].mean())
                          for k in ["R2", "Spearman", "RMSE"]} for m in METHODS},
    pooled_pred=pooled,
    loro_by_region_R2={h: {m: float(pivot_r2.loc[h, m]) for m in METHODS} for h in regions},
    aoa_frac_inside=frac_inside, DI_threshold=float(DI_thresh),
    region_summary=regsum.to_dict(orient="records"),
)
with open(f"{RES}/case_derived.json", "w") as f:
    json.dump(case_out, f, indent=2)
print("\nWrote case_loro.csv, case_grid.csv, case_regions.csv, case_derived.json")
