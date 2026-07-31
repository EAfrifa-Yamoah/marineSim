"""
si_analyses.py
==============
Computes the three supplementary analyses the main text refers to but that were
not part of the headline results:

  S2  Hyperparameter sensitivity of the random forest family (RF, Spatial RF, stRF)
  S3  Full pairwise method contrasts by sample size and by sampling design (Table S1)
  S4  Probabilistic calibration (Brier score and reliability) for all methods

Outputs are written to ../results/si_*.{csv,json}. S3 reuses the stored benchmark
AUCs; S2 and S4 fit models with the engine's new rf_params / return_proba hooks.
"""
import glob, json
import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss
import sim_engine as se

RES = "../results"
METHODS = ["GLM", "GAM", "RF", "SpatialRF", "stRF"]
rng_master = np.random.default_rng(20260701)


# ============================================================ S3 contrasts
def s3_contrasts():
    bench = pd.concat([pd.read_csv(f) for f in glob.glob(f"{RES}/bench_world*.csv")],
                      ignore_index=True)
    gam = pd.concat([pd.read_csv(f) for f in glob.glob(f"{RES}/aux_gam_world*.csv")],
                    ignore_index=True)
    keys = ["world_id", "n", "design", "detection", "n_time", "rep"]
    df = bench.merge(gam[keys + ["GAM"]], on=keys, how="left")

    by_n = df.groupby("n")[METHODS].mean().round(4)
    by_n.to_csv(f"{RES}/si_contrasts_by_n.csv")
    by_design = df[df.n == 100].groupby("design")[METHODS].mean().round(4)
    by_design.to_csv(f"{RES}/si_contrasts_by_design.csv")

    # paired bootstrap: stRF - method within each sample size
    rows = []
    for n, sub in df.groupby("n"):
        for m in ["GLM", "GAM", "RF", "SpatialRF"]:
            d = (sub["stRF"] - sub[m]).dropna().to_numpy()
            bs = [np.mean(rng_master.choice(d, len(d), replace=True)) for _ in range(2000)]
            rows.append(dict(n=int(n), contrast=f"stRF - {m}", delta=round(float(d.mean()), 4),
                             lo=round(float(np.percentile(bs, 2.5)), 4),
                             hi=round(float(np.percentile(bs, 97.5)), 4)))
    pd.DataFrame(rows).to_csv(f"{RES}/si_contrasts_pairwise.csv", index=False)
    print("S3 done:", by_n.shape[0], "sample sizes,", by_design.shape[0], "designs")
    return by_n


# ============================================================ S2 tuning sensitivity
def s2_tuning():
    worlds = list(se.make_worlds().values())
    leaves = [1, 3, 5]; feats = ["sqrt", 0.5, 0.8]
    grid = [(l, mf) for l in leaves for mf in feats]
    rows = []
    for wi, w in enumerate(worlds):
        for n in [50, 100, 200]:
            for rep in range(2):
                seed = 5000 + wi * 100 + n + rep         # same data across configs
                for (leaf, mf) in grid:
                    r = se.fit_eval(w, n, "random", 1.0, 5,
                                    np.random.default_rng(seed),
                                    methods=["RF", "SpatialRF", "stRF"],
                                    rf_params={"min_samples_leaf": leaf, "max_features": mf})
                    if r is None:
                        continue
                    for m in ["RF", "SpatialRF", "stRF"]:
                        rows.append(dict(world=wi, n=n, rep=rep, leaf=leaf, mf=str(mf),
                                         method=m, auc=r.get(m, np.nan)))
    d = pd.DataFrame(rows)
    cell = d.groupby(["method", "leaf", "mf"])["auc"].mean().reset_index()
    cell.to_csv(f"{RES}/si_tuning.csv", index=False)

    # per-method: default (leaf=1, sqrt) vs best config; and ranking check per config
    summ = {}
    for m in ["RF", "SpatialRF", "stRF"]:
        cm = cell[cell.method == m]
        default = float(cm[(cm.leaf == 1) & (cm.mf == "sqrt")]["auc"].iloc[0])
        summ[m] = dict(default=round(default, 4), best=round(float(cm["auc"].max()), 4),
                       worst=round(float(cm["auc"].min()), 4),
                       spread=round(float(cm["auc"].max() - cm["auc"].min()), 4))
    # ranking preserved at every config?
    piv = cell.pivot_table(index=["leaf", "mf"], columns="method", values="auc")
    ranking_ok = bool((piv["stRF"] >= piv["SpatialRF"]).all() and
                      (piv["SpatialRF"] > piv["RF"]).all())
    summ["ranking_preserved_all_configs"] = ranking_ok
    summ["n_configs"] = len(grid)
    json.dump(summ, open(f"{RES}/si_tuning_summary.json", "w"), indent=2)
    print("S2 done: ranking preserved at all", len(grid), "configs:", ranking_ok)
    return summ


# ============================================================ S4 calibration
def s4_calibration():
    worlds = list(se.make_worlds().values())
    pooled = {m: {"p": [], "y": []} for m in METHODS}
    brier = {m: [] for m in METHODS}
    for wi, w in enumerate(worlds):
        for rep in range(3):
            r = se.fit_eval(w, 100, "random", 1.0, 5,
                            np.random.default_rng(8000 + wi * 10 + rep),
                            methods=METHODS, return_proba=True)
            if r is None or "_ytest" not in r:
                continue
            y = r["_ytest"]
            for m in METHODS:
                p = r.get(m + "_p")
                if p is None or np.any(~np.isfinite(p)):
                    continue
                pooled[m]["p"].append(np.asarray(p)); pooled[m]["y"].append(np.asarray(y))
                brier[m].append(brier_score_loss(y, p))

    cal_rows = [dict(method=m, brier_mean=round(float(np.mean(brier[m])), 4),
                     brier_sd=round(float(np.std(brier[m])), 4)) for m in METHODS if brier[m]]
    pd.DataFrame(cal_rows).to_csv(f"{RES}/si_calibration.csv", index=False)

    rel_rows = []
    edges = np.linspace(0, 1, 11)
    for m in METHODS:
        if not pooled[m]["p"]:
            continue
        p = np.concatenate(pooled[m]["p"]); y = np.concatenate(pooled[m]["y"])
        idx = np.digitize(p, edges) - 1
        for b in range(10):
            sel = idx == b
            if sel.sum() >= 20:
                rel_rows.append(dict(method=m, bin_mid=round((edges[b] + edges[b + 1]) / 2, 3),
                                     pred_mean=round(float(p[sel].mean()), 4),
                                     obs_freq=round(float(y[sel].mean()), 4),
                                     count=int(sel.sum())))
    pd.DataFrame(rel_rows).to_csv(f"{RES}/si_reliability.csv", index=False)
    print("S4 done:", cal_rows)
    return cal_rows


if __name__ == "__main__":
    s3_contrasts()
    s2_tuning()
    s4_calibration()
    print("SI ANALYSES DONE")
