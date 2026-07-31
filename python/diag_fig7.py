"""Diagnose the apparent underprediction in Figure 7(b)."""
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold
import case_study as cs

R = "../results"
COVARS = cs.COVARS
agg = pd.read_csv(f"{R}/case_agg.csv")
site = pd.read_csv(f"{R}/case_sites.csv")
grid = pd.read_csv(f"{R}/case_map_grid.csv")
y = np.log1p(agg["dens"].to_numpy())

print("=" * 68)
print("D1. SCALES: what is actually being compared")
print("=" * 68)
print(f"observed site-mean dens : min {site['dens'].min():.1f}  med {site['dens'].median():.1f} "
      f" mean {site['dens'].mean():.1f}  max {site['dens'].max():.1f}")
print(f"observed site-YEAR dens : min {agg['dens'].min():.1f}  med {agg['dens'].median():.1f} "
      f" mean {agg['dens'].mean():.1f}  max {agg['dens'].max():.1f}")
print(f"grid predictions        : min {grid['pred'].min():.1f}  med {grid['pred'].median():.1f} "
      f" mean {grid['pred'].mean():.1f}  max {grid['pred'].max():.1f}")
print(f"panel (b) colour scale is HARD-CODED vmin=10 vmax=40")
print(f"  -> grid preds below vmin (saturate palest): {(grid['pred'] < 10).mean()*100:.1f}%")
print(f"  -> grid preds above vmax (saturate darkest): {(grid['pred'] > 40).mean()*100:.1f}%")
print(f"  -> observed site means above 40: {(site['dens'] > 40).mean()*100:.1f}%")

print()
print("=" * 68)
print("D2. RETRANSFORMATION BIAS: expm1(mean(log1p)) vs mean(raw)")
print("=" * 68)
naive = np.expm1(np.mean(np.log1p(agg["dens"])))
print(f"expm1(mean log1p(dens)) = {naive:.2f}   vs   mean(dens) = {agg['dens'].mean():.2f}")
print(f"  -> naive back-transform sits {100*(1-naive/agg['dens'].mean()):.1f}% BELOW the arithmetic mean")

print()
print("=" * 68)
print("D3. MODEL FIDELITY AT THE SITES (is the model itself low?)")
print("=" * 68)
Xtr, _ = cs.build_design(agg, agg, "stRF")
rf = RandomForestRegressor(n_estimators=600, min_samples_leaf=3, max_features=0.6,
                           random_state=7, n_jobs=1).fit(Xtr, y)
ins_log = rf.predict(Xtr)
ins = np.expm1(ins_log)
# Duan smearing factor on residuals (log scale)
resid = y - ins_log
smear = float(np.mean(np.exp(resid)))
print(f"in-sample: mean obs {agg['dens'].mean():.2f}  mean pred(expm1) {ins.mean():.2f} "
      f" bias {ins.mean()-agg['dens'].mean():+.2f}")
print(f"Duan smearing factor  = {smear:.3f}  (multiply expm1 preds by this to de-bias the mean)")

# out-of-fold by site group
gkf = GroupKFold(n_splits=5)
oof = np.zeros(len(y))
for tr, te in gkf.split(Xtr, y, groups=agg["site"]):
    m = RandomForestRegressor(n_estimators=400, min_samples_leaf=3, max_features=0.6,
                              random_state=7, n_jobs=1).fit(Xtr[tr], y[tr])
    oof[te] = m.predict(Xtr[te])
oof_raw = np.expm1(oof)
print(f"out-of-fold: mean obs {agg['dens'].mean():.2f}  mean pred {oof_raw.mean():.2f} "
      f" bias {oof_raw.mean()-agg['dens'].mean():+.2f}")
print()
print("  bias by region (out-of-fold, raw scale):")
d = agg.copy(); d["oof"] = oof_raw
for rg, sub in d.groupby("region"):
    print(f"    {rg:6s} n={len(sub):4d}  obs {sub['dens'].mean():6.2f}  pred {sub['oof'].mean():6.2f} "
          f" bias {sub['oof'].mean()-sub['dens'].mean():+6.2f}")

print()
print("=" * 68)
print("D4. GRID vs NEAREST OBSERVED SITE (is the GRID low, or the model?)")
print("=" * 68)
sx = site[["longitude", "latitude"]].to_numpy()
gl = grid[["longitude", "latitude"]].to_numpy()
rows = []
for i in range(len(grid)):
    dd = cs.haversine(gl[i, 0], gl[i, 1], sx[:, 0], sx[:, 1])
    j = int(np.argmin(dd))
    rows.append((dd[j], site["dens"].iloc[j], site["region"].iloc[j], grid["pred"].iloc[i]))
gd = pd.DataFrame(rows, columns=["d_km", "obs_site", "region", "pred"])
near = gd[gd["d_km"] <= 3.0]
print(f"grid cells within 3 km of a site: {len(near)}")
print(f"  mean observed at those sites {near['obs_site'].mean():.2f}   "
      f"mean grid pred {near['pred'].mean():.2f}   bias {near['pred'].mean()-near['obs_site'].mean():+.2f}")
print()
print("  by region (grid cells within 3 km of a site):")
for rg, sub in near.groupby("region"):
    print(f"    {rg:6s} n={len(sub):4d}  obs_site {sub['obs_site'].mean():6.2f}  "
          f"grid_pred {sub['pred'].mean():6.2f}  bias {sub['pred'].mean()-sub['obs_site'].mean():+6.2f}")

print()
print("=" * 68)
print("D5. TEMPORAL: is the grid's year (2024) simply a low year?")
print("=" * 68)
print(f"grid year assigned = {grid['year'].iloc[0]}  (temporal feature = region mean of year-1)")
yr = agg.groupby("year")["dens"].mean()
print("  regional mean density by year (last 8):")
for yy, v in yr.tail(8).items():
    print(f"    {int(yy)}: {v:6.2f}")
print(f"  long-term mean across all years: {agg['dens'].mean():.2f}")
print()
print("  region x recent-year means vs region long-term mean:")
for rg, sub in agg.groupby("region"):
    lt = sub["dens"].mean()
    recent = sub[sub["year"] >= 2022]["dens"]
    rv = recent.mean() if len(recent) else np.nan
    print(f"    {rg:6s} long-term {lt:6.2f}   2022+ {rv:6.2f}   diff {rv-lt:+6.2f}")
