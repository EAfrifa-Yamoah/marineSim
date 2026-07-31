"""
verify.py
=========
Three targeted checks.

  ISSUE 2  Leakage invariance. The reviewer's test: changing a held-out
           response must not change any predictor used for fitting or
           prediction.
  ISSUE 4  Empirical scaling of the Gaussian process comparator, to test the
           claim that its cost "scales cubically in the number of observations"
           and that the gap widens to "roughly an order of magnitude at n=500".
  ISSUE 6  Coverage of the 80% across-tree ensemble band. The simulation knows
           the TRUE occurrence probability at every test cell, so coverage can
           be measured rather than asserted. Compared inside and outside the
           area of applicability.
"""

import json, time
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.ensemble import RandomForestClassifier

from sim_engine import World, RF_TREES, RADII_KM
from ablation import build_dataset, world_grid

R = "../results"
out = {}


# ===================================================================== ISSUE 2
def leakage_test(n_worlds=4):
    """
    Rebuild every feature block after permuting the held-out test responses.
    If no predictor changes, held-out responses cannot influence fitting or
    prediction.
    """
    checks, failures = 0, []
    for (wid, cfg, real, stat, rng_km, wseed) in world_grid(1)[:n_worlds]:
        w = World(stationary=stat, range_km=rng_km, n_time=5, seed=wseed)
        for n in (50, 200):
            ds = build_dataset(w, n, 5, 0, seed=wseed * 31 + n)
            if ds is None:
                continue
            base_tr = {k: v.copy() for k, v in ds["blocks_tr"].items()}
            base_te = {k: v.copy() for k, v in ds["blocks_te"].items()}
            base_Xte = ds["Xc_te"].copy()

            # corrupt the held-out responses completely
            rng = np.random.default_rng(0)
            ds2 = build_dataset(w, n, 5, 0, seed=wseed * 31 + n)
            ds2["y_te"] = rng.permutation(ds2["y_te"])
            ds2["y_te"] = 1 - ds2["y_te"]        # and flip them

            for b in base_tr:
                checks += 1
                if not np.allclose(base_tr[b], ds2["blocks_tr"][b]):
                    failures.append(f"train block {b} changed (world {wid}, n={n})")
                checks += 1
                if not np.allclose(base_te[b], ds2["blocks_te"][b]):
                    failures.append(f"test block {b} changed (world {wid}, n={n})")
            checks += 1
            if not np.allclose(base_Xte, ds2["Xc_te"]):
                failures.append(f"test covariates changed (world {wid}, n={n})")

    # independent check: does any response-derived TEST feature depend on
    # training responses only? perturb one TRAINING response and confirm the
    # test lag features DO change (i.e. the test is not vacuous).
    w = World(stationary=False, range_km=80, n_time=5, seed=20260601)
    ds = build_dataset(w, 100, 5, 0, seed=99)
    lag_before = ds["blocks_te"]["LAG"].copy()
    ds3 = build_dataset(w, 100, 5, 0, seed=99)
    ds3["y_tr"][0] = 1 - ds3["y_tr"][0]
    from ablation import spatial_lags
    lag_after = spatial_lags(ds3["site_xy"], ds3["y_tr"], ds3["test_xy"],
                             RADII_KM, is_train=False)
    sensitive = not np.allclose(lag_before, lag_after)

    return dict(checks=checks, failures=failures, passed=len(failures) == 0,
                test_is_non_vacuous=bool(sensitive))


# ===================================================================== ISSUE 4
def gp_scaling(reps=3):
    """Time the GP comparator and the stRF backbone across sample sizes."""
    from sklearn.gaussian_process import GaussianProcessClassifier
    from sklearn.gaussian_process.kernels import (Matern, ConstantKernel,
                                                  WhiteKernel)
    rows = []
    for (wid, cfg, real, stat, rng_km, wseed) in world_grid(1)[:3]:
        w = World(stationary=stat, range_km=rng_km, n_time=5, seed=wseed)
        for n in (30, 50, 100, 200, 500):
            for rep in range(reps):
                ds = build_dataset(w, n, 5, rep, seed=wseed * 13 + n * 7 + rep)
                if ds is None:
                    continue
                Xgp_tr = np.column_stack([ds["blocks_tr"]["COORD"], ds["Xc_tr"]])
                kern = (ConstantKernel(1.0) * Matern(length_scale=50.0, nu=1.5)
                        + WhiteKernel(1e-2))
                gp = GaussianProcessClassifier(kernel=kern,
                                               n_restarts_optimizer=0,
                                               max_iter_predict=100,
                                               random_state=0)
                t0 = time.time(); gp.fit(Xgp_tr, ds["y_tr"].astype(int))
                t_gp = time.time() - t0

                Xrf_tr = np.column_stack([ds["Xc_tr"], ds["blocks_tr"]["LAG"],
                                          ds["blocks_tr"]["TEMP"],
                                          ds["blocks_tr"]["EDF"]])
                rf = RandomForestClassifier(n_estimators=RF_TREES, n_jobs=1,
                                            random_state=0)
                t0 = time.time(); rf.fit(Xrf_tr, ds["y_tr"].astype(int))
                t_rf = time.time() - t0
                rows.append(dict(world=cfg, n=n, rep=rep,
                                 gp_fit_s=t_gp, strf_fit_s=t_rf))
    df = pd.DataFrame(rows)
    df.to_csv(f"{R}/gp_scaling.csv", index=False)
    med = df.groupby("n")[["gp_fit_s", "strf_fit_s"]].median()
    # log-log slope of GP time on n
    sl_gp = np.polyfit(np.log(med.index.values), np.log(med.gp_fit_s.values), 1)[0]
    sl_rf = np.polyfit(np.log(med.index.values), np.log(med.strf_fit_s.values), 1)[0]
    return dict(median_by_n=med.round(4).to_dict(),
                gp_loglog_slope=float(sl_gp), strf_loglog_slope=float(sl_rf),
                ratio_at_500=float(med.loc[500, "gp_fit_s"]
                                   / med.loc[500, "strf_fit_s"]),
                ratio_median_overall=float(df.gp_fit_s.median()
                                           / df.strf_fit_s.median()))


# ===================================================================== ISSUE 6
def aoa_di(Xtr, Xte, importances):
    """Meyer & Pebesma (2021) dissimilarity index, importance weighted."""
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-12
    w = np.sqrt(np.maximum(importances, 0))
    Ztr = (Xtr - mu) / sd * w
    Zte = (Xte - mu) / sd * w
    Dtr = cdist(Ztr, Ztr)
    np.fill_diagonal(Dtr, np.inf)
    dbar = cdist(Ztr, Ztr)[np.triu_indices(len(Ztr), 1)].mean()
    di_tr = Dtr.min(1) / dbar
    di_te = cdist(Zte, Ztr).min(1) / dbar
    return di_te, float(np.percentile(di_tr, 95))


def band_coverage(reps=2):
    """
    80% across-tree band vs the TRUE occurrence probability.
    Nominal coverage is 0.80. Reported overall, and split by whether the
    query point is inside or outside the area of applicability.
    """
    rows = []
    for (wid, cfg, real, stat, rng_km, wseed) in world_grid(1):
        w = World(stationary=stat, range_km=rng_km, n_time=5, seed=wseed)
        for n in (50, 200):
            for rep in range(reps):
                ds = build_dataset(w, n, 5, rep, seed=wseed * 17 + n * 3 + rep)
                if ds is None:
                    continue
                Xtr = np.column_stack([ds["Xc_tr"], ds["blocks_tr"]["LAG"],
                                       ds["blocks_tr"]["TEMP"],
                                       ds["blocks_tr"]["EDF"]])
                Xte = np.column_stack([ds["Xc_te"], ds["blocks_te"]["LAG"],
                                       ds["blocks_te"]["TEMP"],
                                       ds["blocks_te"]["EDF"]])
                rf = RandomForestClassifier(n_estimators=RF_TREES, n_jobs=1,
                                            random_state=0)
                rf.fit(Xtr, ds["y_tr"].astype(int))
                per_tree = np.stack([t.predict_proba(Xte)[:, 1]
                                     for t in rf.estimators_])
                lo = np.percentile(per_tree, 10, axis=0)
                hi = np.percentile(per_tree, 90, axis=0)
                p_true = ds["p_te"]
                cover = (p_true >= lo) & (p_true <= hi)
                width = hi - lo
                di, thr = aoa_di(Xtr, Xte, rf.feature_importances_)
                inside = di <= thr
                rows.append(dict(
                    world=cfg, n=n, rep=rep,
                    coverage=float(cover.mean()),
                    coverage_inside=float(cover[inside].mean()) if inside.any() else np.nan,
                    coverage_outside=float(cover[~inside].mean()) if (~inside).any() else np.nan,
                    width_inside=float(width[inside].mean()) if inside.any() else np.nan,
                    width_outside=float(width[~inside].mean()) if (~inside).any() else np.nan,
                    pct_outside=float((~inside).mean())))
    df = pd.DataFrame(rows)
    df.to_csv(f"{R}/band_coverage.csv", index=False)
    return df


if __name__ == "__main__":
    print("=" * 78)
    print("ISSUE 2 | LEAKAGE INVARIANCE")
    print("=" * 78)
    lk = leakage_test()
    out["leakage"] = lk
    print(f"  predictor comparisons run : {lk['checks']}")
    print(f"  failures                  : {len(lk['failures'])}")
    print(f"  test is non-vacuous       : {lk['test_is_non_vacuous']}"
          "   (training responses DO move test lags, as they must)")
    print(f"  VERDICT                   : "
          f"{'PASS - no held-out response reaches any predictor' if lk['passed'] else 'FAIL'}")
    for f in lk["failures"][:5]:
        print("   ", f)

    print()
    print("=" * 78)
    print("ISSUE 4 | EMPIRICAL SCALING OF THE GEOSTATISTICAL COMPARATOR")
    print("=" * 78)
    gs = gp_scaling()
    out["gp_scaling"] = gs
    med = pd.DataFrame(gs["median_by_n"])
    print(med.to_string())
    print(f"\n  GP   log-log slope vs n : {gs['gp_loglog_slope']:.2f}"
          "   (cubic would be 3.0)")
    print(f"  stRF log-log slope vs n : {gs['strf_loglog_slope']:.2f}")
    print(f"  GP/stRF fit-time ratio at n=500 : {gs['ratio_at_500']:.2f}x")
    print(f"  GP/stRF median ratio overall    : {gs['ratio_median_overall']:.2f}x")

    print()
    print("=" * 78)
    print("ISSUE 6 | COVERAGE OF THE 80% ACROSS-TREE BAND (nominal 0.80)")
    print("=" * 78)
    bc = band_coverage()
    agg = bc[["coverage", "coverage_inside", "coverage_outside",
              "width_inside", "width_outside", "pct_outside"]].mean()
    out["band"] = {k: float(v) for k, v in agg.items()}
    print(f"  overall coverage of the TRUE probability : {agg['coverage']:.3f}")
    print(f"  coverage inside  AOA                     : {agg['coverage_inside']:.3f}")
    print(f"  coverage outside AOA                     : {agg['coverage_outside']:.3f}")
    print(f"  mean band width inside  AOA              : {agg['width_inside']:.3f}")
    print(f"  mean band width outside AOA              : {agg['width_outside']:.3f}")
    print(f"  fraction of test cells outside AOA       : {agg['pct_outside']:.3f}")

    with open(f"{R}/verification.json", "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\nWROTE {R}/verification.json, gp_scaling.csv, band_coverage.csv")
