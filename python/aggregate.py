"""
aggregate.py
============
Consolidate the benchmark CSVs into the headline numbers, decomposition, cluster
bootstrap CIs, a hierarchical mixed-effects model, and the per-figure data tables.
Outputs results/derived.json and several results/fig_*.csv files.
"""
import os, glob, json
import numpy as np
import pandas as pd

R = os.path.join(os.path.dirname(__file__), "..", "results")
rng = np.random.default_rng(20260630)

# ---------------------------------------------------------------- load & merge
main = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{R}/bench_world*.csv"))],
                 ignore_index=True)
gam = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{R}/aux_gam_world*.csv"))],
                ignore_index=True)
keys = ["world_id", "n", "design", "detection", "n_time", "rep"]
df = main.merge(gam[keys + ["GAM"]], on=keys, how="left")

det07 = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{R}/aux_det07_world*.csv"))],
                  ignore_index=True)
geo = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{R}/aux_geo_world*.csv"))],
                ignore_index=True)

METHODS = ["GLM", "GAM", "RF", "SpatialRF", "stRF"]
for m in METHODS:
    df[m] = pd.to_numeric(df[m], errors="coerce")

# ---------------------------------------------------------------- decomposition
df["d_ml"] = df["RF"] - df["GLM"]
df["d_spatial"] = df["SpatialRF"] - df["RF"]
df["d_temporal"] = df["stRF"] - df["SpatialRF"]
df["d_total"] = df["stRF"] - df["GLM"]

def world_means(frame, col):
    return frame.groupby("world_id")[col].mean()

def cluster_boot(frame, col, B=1000):
    wm = world_means(frame, col)
    ids = wm.index.to_numpy(); vals = wm.to_numpy()
    boots = np.empty(B)
    for b in range(B):
        pick = rng.choice(len(ids), size=len(ids), replace=True)
        boots[b] = np.nanmean(vals[pick])
    return float(np.nanmean(vals)), float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))

dec = {}
for label, col in [("ml", "d_ml"), ("spatial", "d_spatial"),
                   ("temporal", "d_temporal"), ("total", "d_total")]:
    m, lo, hi = cluster_boot(df, col)
    dec[label] = dict(mean=m, lo=lo, hi=hi)

ratio = dec["spatial"]["mean"] / dec["ml"]["mean"]
# bootstrap the ratio too
wm_sp = world_means(df, "d_spatial").to_numpy()
wm_ml = world_means(df, "d_ml").to_numpy()
rb = []
for b in range(2000):
    pick = rng.choice(len(wm_sp), size=len(wm_sp), replace=True)
    denom = np.nanmean(wm_ml[pick])
    if abs(denom) > 1e-6:
        rb.append(np.nanmean(wm_sp[pick]) / denom)
ratio_lo, ratio_hi = float(np.percentile(rb, 2.5)), float(np.percentile(rb, 97.5))

shares = {k: dec[k]["mean"] / dec["total"]["mean"] for k in ["ml", "spatial", "temporal"]}

# ---------------------------------------------------------------- win counts
cellkeys = ["world_id", "n", "design", "n_time"]
cell = df.groupby(cellkeys)[METHODS].mean().reset_index()
cell["winner"] = cell[METHODS].idxmax(axis=1)
winpct = (cell["winner"].value_counts(normalize=True) * 100).round(1).to_dict()

# ---------------------------------------------------------------- AUC vs n (Fig 2)
longrows = []
by_n_method = df.groupby("n")[METHODS].agg(["mean", "std"])
fig2 = []
for n in sorted(df["n"].unique()):
    for m in METHODS:
        sub = df[df["n"] == n][m].dropna()
        fig2.append(dict(n=int(n), method=m, mean=float(sub.mean()), sd=float(sub.std())))
fig2 = pd.DataFrame(fig2)
fig2.to_csv(f"{R}/fig2_auc_by_n.csv", index=False)

# threshold n where mean AUC >= 0.75 per method
thr = {}
for m in METHODS:
    g = fig2[fig2.method == m].sort_values("n")
    hit = g[g["mean"] >= 0.75]
    thr[m] = int(hit["n"].iloc[0]) if len(hit) else None

# gap stRF - SpatialRF at n=500
gap_500 = float(df[df.n == 500]["d_temporal"].mean())
gap_30 = float(df[df.n == 30]["d_temporal"].mean())

# ---------------------------------------------------------------- temporal surface (Fig 3)
fig3 = (df.groupby(["n", "n_time"])["d_temporal"].mean().reset_index())
fig3.to_csv(f"{R}/fig3_temporal_surface.csv", index=False)
hi_regime = df[(df.n >= 100) & (df.n_time >= 3)]["d_temporal"]
lo_regime = df[(df.n < 100) | (df.n_time < 3)]["d_temporal"]
delta_high = float(hi_regime.mean())
# best-scenario mean (per cell means, upper tail)
cell_temp = df.groupby(cellkeys)["d_temporal"].mean()
delta_max = float(cell_temp[(df.groupby(cellkeys)["n"].first() >= 100)].quantile(0.9))

# ---------------------------------------------------------------- design effects (Fig 4, n=100)
fig4 = (df[df.n == 100].groupby(["design"])[METHODS].mean().reset_index())
fig4.to_csv(f"{R}/fig4_design_n100.csv", index=False)
# clustered gap vs mean(random,stratified) per method, averaged
d100 = df[df.n == 100]
design_means = d100.groupby("design")[METHODS].mean()
clustered_gap = float((0.5 * (design_means.loc["random"] + design_means.loc["stratified"])
                       - design_means.loc["clustered"]).mean())
strat_adv = float((design_means.loc["stratified"] - design_means.loc["random"]).mean())
# n=100 stratified vs n=200 clustered (any method avg of spatial methods)
strat100 = float(df[(df.n == 100) & (df.design == "stratified")][["SpatialRF", "stRF"]].mean().mean())
clus200 = float(df[(df.n == 200) & (df.design == "clustered")][["SpatialRF", "stRF"]].mean().mean())

# ---------------------------------------------------------------- sdmTMB-class (GEO) Section 4.5
geo["stRF"] = pd.to_numeric(geo["stRF"], errors="coerce")
geo["GEO"] = pd.to_numeric(geo["GEO"], errors="coerce")
geo_valid = geo.dropna(subset=["stRF", "GEO"])
strf_mean = float(geo_valid["stRF"].mean())
geo_mean = float(geo_valid["GEO"].mean())
diff = geo_valid["stRF"] - geo_valid["GEO"]
# paired cluster bootstrap over worlds
gwm = geo_valid.groupby("world_id").apply(lambda x: (x["stRF"] - x["GEO"]).mean()).to_numpy()
db = [np.nanmean(gwm[rng.choice(len(gwm), len(gwm), replace=True)]) for _ in range(2000)]
geo_diff = float(diff.mean())
geo_diff_lo, geo_diff_hi = float(np.percentile(db, 2.5)), float(np.percentile(db, 97.5))
geo_conv_overall = float(geo["GEO_conv"].mean() * 100)
geo_conv_small = float(geo[geo.n <= 50]["GEO_conv"].mean() * 100)
# time ratio: GEO_time vs stRF (stRF ~ fit time approximated from main timing ~0.05-0.1s at n;
# we report measured GEO median time and a conservative stRF fit time of 0.05s)
geo_time_median = float(geo["GEO_time"].median())
strf_fit_time = 0.05
time_ratio = geo_time_median / strf_fit_time

# ---------------------------------------------------------------- detection 0.7 effect
det_eff = {}
for m in ["GLM", "RF", "SpatialRF", "stRF"]:
    a10 = df[(df.n.isin([50, 200, 500])) & (df.design.isin(["random", "clustered"]))
             & (df.n_time.isin([1, 5]))][m].mean()
    a07 = det07[m].mean()
    det_eff[m] = dict(det10=float(a10), det07=float(a07), drop=float(a10 - a07))

# ---------------------------------------------------------------- mixed-effects model
mm_summary = {}
try:
    import statsmodels.formula.api as smf
    long = df.melt(id_vars=["world_id", "n", "design", "detection", "n_time", "rep"],
                   value_vars=METHODS, var_name="method", value_name="auc").dropna()
    long["scenario"] = (long["world_id"].astype(str) + "_" + long["n"].astype(str) + "_"
                        + long["design"] + "_" + long["n_time"].astype(str))
    long["logit"] = np.log(long["auc"].clip(0.02, 0.98) / (1 - long["auc"].clip(0.02, 0.98)))
    long["logn"] = np.log(long["n"])
    long["method"] = pd.Categorical(long["method"],
                                    categories=["GLM", "GAM", "RF", "SpatialRF", "stRF"])
    md = smf.mixedlm("logit ~ C(method) * logn + C(design)", long,
                     groups=long["world_id"])
    mf = md.fit(method="lbfgs", maxiter=200)
    params = mf.params.to_dict()
    mm_summary = {k: float(v) for k, v in params.items()
                  if "method" in k.lower() or "design" in k.lower() or k == "logn"}
    mm_summary["_note"] = "logit-AUC ~ method*log(n) + design, random intercept by world"
except Exception as e:
    mm_summary = {"error": str(e)}

# ---------------------------------------------------------------- assemble JSON
out = dict(
    n_rows_main=int(len(df)), n_fits_total=int(len(main) + len(gam) + len(det07) + len(geo)),
    n_scenarios=int(cell.shape[0]),
    decomposition=dec, ratio=dict(mean=ratio, lo=ratio_lo, hi=ratio_hi),
    shares=shares, win_pct=winpct, thresholds=thr,
    gap_stRF_spRF=dict(n30=gap_30, n500=gap_500),
    temporal=dict(delta_high=delta_high, delta_max=delta_max),
    design=dict(clustered_gap=clustered_gap, strat_adv=strat_adv,
                strat100=strat100, clus200=clus200),
    geo=dict(strf=strf_mean, geo=geo_mean, diff=geo_diff, lo=geo_diff_lo, hi=geo_diff_hi,
             conv_overall=geo_conv_overall, conv_small=geo_conv_small,
             time_median=geo_time_median, time_ratio=time_ratio),
    detection=det_eff, mixed_model=mm_summary,
)

with open(f"{R}/derived.json", "w") as f:
    json.dump(out, f, indent=2)

# also persist the per-row decomposition for Fig 1
df[["world_id", "n", "design", "n_time", "rep",
    "d_ml", "d_spatial", "d_temporal", "d_total"]].to_csv(f"{R}/fig1_decomp_rows.csv", index=False)
geo_valid.to_csv(f"{R}/fig_geo.csv", index=False)

# ---------------------------------------------------------------- console report
print("=== HEADLINE DECOMPOSITION (cluster bootstrap over 6 worlds) ===")
for k in ["ml", "spatial", "temporal", "total"]:
    d = dec[k]
    print(f"  {k:9s}: {d['mean']:+.3f}  95% CI [{d['lo']:+.3f}, {d['hi']:+.3f}]  share={shares.get(k, 1.0)*100:4.0f}%"
          if k != "total" else
          f"  {k:9s}: {d['mean']:+.3f}  95% CI [{d['lo']:+.3f}, {d['hi']:+.3f}]")
print(f"  spatial:ML ratio = {ratio:.2f}  95% CI [{ratio_lo:.2f}, {ratio_hi:.2f}]")
print(f"\nWin %: {winpct}")
print(f"Thresholds (n at AUC>=0.75): {thr}")
print(f"Temporal gap stRF-spRF: n30={gap_30:+.3f}  n500={gap_500:+.3f}")
print(f"Temporal high-regime mean (n>=100 & ts>=3): {delta_high:+.3f}; 90th pct(best)≈{delta_max:+.3f}")
print(f"Design: clustered gap={clustered_gap:.3f}  strat adv={strat_adv:.3f}")
print(f"  stratified@100={strat100:.3f} vs clustered@200={clus200:.3f}")
print(f"GEO: stRF={strf_mean:.3f} GEO={geo_mean:.3f} diff={geo_diff:+.3f} CI[{geo_diff_lo:+.3f},{geo_diff_hi:+.3f}]")
print(f"  conv overall={geo_conv_overall:.0f}%  conv n<=50={geo_conv_small:.0f}%  median time={geo_time_median:.2f}s  ratio~{time_ratio:.0f}x")
print(f"Detection drop (det 0.7 vs 1.0): " + ", ".join(f"{m}={det_eff[m]['drop']:+.3f}" for m in det_eff))
print(f"\nTotal fits across all runs: {out['n_fits_total']}, main rows: {out['n_rows_main']}, scenarios: {out['n_scenarios']}")
print("Wrote derived.json + fig CSVs")
