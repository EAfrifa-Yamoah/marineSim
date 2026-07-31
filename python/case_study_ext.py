"""
case_study_ext.py
=================
Additional empirical analyses for the seagrass case study, layered on top of the
leave-one-region-out extrapolation already computed in case_study.py:

  (A) Interpolation skill: random 10-fold cross validation over site-years, so the
      held-out observations sit inside the sampled domain. Contrasting this with
      the leave-one-region-out (extrapolation) skill quantifies how far predictive
      performance falls when a model is pushed beyond the latitudes it has seen.
  (B) Variable importance: permutation importance for the stRF feature set, grouped
      into interpretable families, to show which signals drive shoot density.
  (C) Prediction uncertainty: the spread of the random forest ensemble along the
      latitudinal transect, giving an 80% prediction band that can be read against
      the area of applicability.
  (D) Temporal trends: the per-region slope of mean shoot density over 2003-2024.
"""
import os, json, io, contextlib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold, train_test_split
from sklearn.inspection import permutation_importance

with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
    import case_study as cs   # reuse data + feature engineering, suppress its prints

RES = cs.RES
agg = cs.agg.copy()
COVARS = cs.COVARS
regions = cs.regions
rng = np.random.default_rng(99)

# ====================================================== (A) interpolation CV
METHODS = ["GLM", "RF", "SpatialRF", "stRF"]
kf = KFold(n_splits=10, shuffle=True, random_state=11)
store = {m: {"y": [], "yhat": []} for m in METHODS}
idx = np.arange(len(agg))
for tr_i, te_i in kf.split(idx):
    tr = agg.iloc[tr_i].reset_index(drop=True)
    te = agg.iloc[te_i].reset_index(drop=True)
    for m in METHODS:
        try:
            yhat = cs.fit_predict(tr, te, m)
        except Exception:
            yhat = np.full(len(te), np.nan)
        store[m]["y"].append(te["dens"].to_numpy())
        store[m]["yhat"].append(yhat)

interp = {}
for m in METHODS:
    y = np.concatenate(store[m]["y"]); yh = np.concatenate(store[m]["yhat"])
    ok = np.isfinite(yh)
    interp[m] = dict(zip(["R2", "Spearman", "RMSE"], cs.metrics(y[ok], yh[ok])))

extrap = {m: cs.pooled[m] for m in METHODS}   # from LORO in case_study.py

ie = pd.DataFrame([
    dict(method=m, regime="Interpolation (random CV)", **interp[m]) for m in METHODS
] + [
    dict(method=m, regime="Extrapolation (leave-region-out)", **extrap[m]) for m in METHODS
])
ie.to_csv(f"{RES}/case_interp_extrap.csv", index=False)

# ====================================================== (B) variable importance
Xall, _ = cs.build_design(agg, agg, "stRF")
n_edf = Xall.shape[1] - len(COVARS) - len(cs.RADII_KM) - 1 - 2
feat_names = (list(COVARS) + [f"lag{r}" for r in cs.RADII_KM] +
              [f"edf{i+1}" for i in range(n_edf)] + ["tlag", "lon", "lat"])
y_all = np.log1p(agg["dens"].to_numpy())
Xtr, Xte, ytr, yte = train_test_split(Xall, y_all, test_size=0.25, random_state=3)
rf_imp = RandomForestRegressor(n_estimators=500, min_samples_leaf=3,
                               max_features=0.6, random_state=7, n_jobs=1).fit(Xtr, ytr)
pi = permutation_importance(rf_imp, Xte, yte, n_repeats=20, random_state=5, n_jobs=1)
imp = np.clip(pi.importances_mean, 0, None)

groups = {
    "Spatial lags": [f"lag{r}" for r in cs.RADII_KM],
    "Distance fields": [f"edf{i+1}" for i in range(n_edf)],
    "Latitude": ["lat"], "Longitude": ["lon"],
    "Depth": ["depth"], "Light (kd490)": ["kd490"], "Habitat rank": ["htyperank"],
    "Summer SST": ["meansummermIST", "maxsummermIST", "CVsummermIST"],
    "Temporal lag": ["tlag"],
}
grouped = {g: float(sum(imp[feat_names.index(f)] for f in fs)) for g, fs in groups.items()}
gtot = sum(grouped.values()) or 1.0
grouped_pct = {g: 100 * v / gtot for g, v in grouped.items()}
imp_df = pd.DataFrame(sorted(grouped_pct.items(), key=lambda kv: -kv[1]),
                      columns=["feature", "importance_pct"])
imp_df.to_csv(f"{RES}/case_importance.csv", index=False)

# ====================================================== (C) prediction uncertainty
grid = cs.grid.copy()
Xtr_full, Xgrid = cs.build_design(agg, grid, "stRF")
rf_unc = RandomForestRegressor(n_estimators=600, min_samples_leaf=3,
                               max_features=0.6, random_state=7, n_jobs=1).fit(
    Xtr_full, np.log1p(agg["dens"].to_numpy()))
per_tree = np.stack([est.predict(Xgrid) for est in rf_unc.estimators_])  # log scale
pred_mean = np.expm1(per_tree.mean(0))
pred_lo = np.expm1(np.percentile(per_tree, 10, axis=0))
pred_hi = np.expm1(np.percentile(per_tree, 90, axis=0))

cg = pd.read_csv(f"{RES}/case_grid.csv")   # has latitude, DI, inside_aoa
cg = cg.sort_values("latitude").reset_index(drop=True)
order = np.argsort(grid["latitude"].to_numpy())
cg["pred_mean"] = pred_mean[order]
cg["pred_lo"] = pred_lo[order]
cg["pred_hi"] = pred_hi[order]
cg.to_csv(f"{RES}/case_grid_unc.csv", index=False)
band_in = float(np.mean((pred_hi - pred_lo)[cg.sort_values('latitude').index]))
# mean band width inside vs outside AOA
inside = cg["inside_aoa"].to_numpy().astype(bool)
width = (cg["pred_hi"] - cg["pred_lo"]).to_numpy()
band_inside = float(width[inside].mean())
band_outside = float(width[~inside].mean()) if (~inside).any() else float("nan")

# ====================================================== (D) temporal trends
trends = {}
for reg in regions:
    sub = agg[agg.region == reg]
    yr = sub.groupby("year")["dens"].mean()
    if len(yr) >= 3:
        slope = float(np.polyfit(yr.index.values, yr.values, 1)[0])
        trends[reg] = dict(slope_per_yr=slope, n_years=int(len(yr)),
                           mean=float(yr.mean()))
# pooled across all parks
yr_all = agg.groupby("year")["dens"].mean()
pooled_slope = float(np.polyfit(yr_all.index.values, yr_all.values, 1)[0])

# ====================================================== save + report
out = dict(
    interp=interp, extrap=extrap,
    interp_extrap_drop_spearman={m: interp[m]["Spearman"] - extrap[m]["Spearman"] for m in METHODS},
    importance_pct=grouped_pct,
    band_inside=band_inside, band_outside=band_outside,
    trends=trends, pooled_slope_per_yr=pooled_slope,
)
with open(f"{RES}/case_ext.json", "w") as f:
    json.dump(out, f, indent=2)

print("=== (A) Interpolation vs extrapolation (Spearman / RMSE) ===")
for m in METHODS:
    print(f"  {m:10s} interp rho={interp[m]['Spearman']:+.3f} RMSE={interp[m]['RMSE']:5.2f} "
          f"| extrap rho={extrap[m]['Spearman']:+.3f} RMSE={extrap[m]['RMSE']:5.2f}")
print("\n=== (B) Grouped stRF permutation importance (%) ===")
for _, r in imp_df.iterrows():
    print(f"  {r['feature']:16s} {r['importance_pct']:5.1f}%")
print(f"\n=== (C) 80% band width: inside AOA {band_inside:.1f} vs outside {band_outside:.1f} shoots ===")
print(f"=== (D) Pooled temporal slope: {pooled_slope:+.2f} shoots/yr; per-region: "
      + ", ".join(f"{k} {v['slope_per_yr']:+.2f}" for k, v in trends.items()))
print("\nWrote case_interp_extrap.csv, case_importance.csv, case_grid_unc.csv, case_ext.json")
