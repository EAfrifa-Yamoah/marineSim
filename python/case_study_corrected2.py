"""
case_study_corrected2.py
========================
D. Interpolation (random CV) vs site-grouped CV vs LORO on one scale.
E. Rolling origin forecast skill: stRF against persistence and climatology,
   for levels and for anomalies.
F. Prediction intervals: across-tree ensemble spread (what Figure S5 plots)
   against split conformal and quantile regression forest, with empirical
   coverage under interpolation and under extrapolation.
G. Area of applicability, and interval behaviour inside vs outside it.
H. Leakage assertion for the case study feature builder.
"""
import json, warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.spatial.distance import cdist
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold, GroupKFold

from case_study_corrected import (load, build_blocks, design, metrics,
                                  COVARS, BLOCKS, REGIONS, RF_KW, RES)

warnings.filterwarnings("ignore")
out = {}
FULL = frozenset(BLOCKS)


def fit_full(tr, te):
    b_tr, b_q, diag = build_blocks(tr, te)
    Xtr, Xq = design(tr, te, b_tr, b_q, FULL)
    ytr = np.log1p(tr["dens"].to_numpy())
    m = RandomForestRegressor(**RF_KW).fit(Xtr, ytr)
    return m, Xtr, Xq, ytr, diag


# ============================================================ D. CV regimes
def cv_regimes(agg):
    rows = []
    ladder = {"RF": frozenset(), "SpatialRF": frozenset({"COORD"}),
              "stRF": FULL}

    # random 10-fold (interpolation)
    kf = KFold(n_splits=10, shuffle=True, random_state=1)
    store = {k: [] for k in ladder}; obs = []
    for tr_i, te_i in kf.split(agg):
        tr, te = agg.iloc[tr_i], agg.iloc[te_i]
        b_tr, b_q, _ = build_blocks(tr, te)
        ytr = np.log1p(tr["dens"].to_numpy())
        for name, S in ladder.items():
            Xtr, Xq = design(tr, te, b_tr, b_q, S)
            m = RandomForestRegressor(**RF_KW).fit(Xtr, ytr)
            store[name].append(np.expm1(m.predict(Xq)))
        obs.append(te["dens"].to_numpy())
    o = np.concatenate(obs)
    for name in ladder:
        r2, rho, rmse = metrics(o, np.concatenate(store[name]))
        rows.append(dict(regime="random 10-fold", method=name,
                         R2=r2, Spearman=rho, RMSE=rmse))

    # site-grouped (new site, same region)
    gk = GroupKFold(n_splits=10)
    store = {k: [] for k in ladder}; obs = []
    for tr_i, te_i in gk.split(agg, groups=agg["site_uid"]):
        tr, te = agg.iloc[tr_i], agg.iloc[te_i]
        b_tr, b_q, _ = build_blocks(tr, te)
        ytr = np.log1p(tr["dens"].to_numpy())
        for name, S in ladder.items():
            Xtr, Xq = design(tr, te, b_tr, b_q, S)
            m = RandomForestRegressor(**RF_KW).fit(Xtr, ytr)
            store[name].append(np.expm1(m.predict(Xq)))
        obs.append(te["dens"].to_numpy())
    o = np.concatenate(obs)
    for name in ladder:
        r2, rho, rmse = metrics(o, np.concatenate(store[name]))
        rows.append(dict(regime="site-grouped", method=name,
                         R2=r2, Spearman=rho, RMSE=rmse))
    return pd.DataFrame(rows)


# ==================================================== E. rolling origin
def rolling_origin(agg, horizons=(1, 2, 3, 5)):
    rows = []
    years = np.sort(agg.year.unique())
    origins = [y for y in years if (agg.year <= y).sum() >= 200 and y <= years[-2]]
    for origin in origins:
        tr = agg[agg.year <= origin].copy()
        if tr.empty:
            continue
        clim = tr.groupby("site_uid")["dens"].mean()
        last = (tr.sort_values("year").groupby("site_uid")["dens"].last())
        for h in horizons:
            te = agg[agg.year == origin + h].copy()
            te = te[te.site_uid.isin(clim.index)]
            if len(te) < 15:
                continue
            m, _, Xq, _, _ = fit_full(tr, te)
            pred_strf = np.expm1(m.predict(Xq))
            pred_clim = clim.reindex(te.site_uid).to_numpy()
            pred_pers = last.reindex(te.site_uid).to_numpy()
            y = te["dens"].to_numpy()
            base = clim.reindex(te.site_uid).to_numpy()   # site mean
            for nm, p in [("stRF", pred_strf), ("climatology", pred_clim),
                          ("persistence", pred_pers)]:
                ok = np.isfinite(p)
                rho = spearmanr(y[ok], p[ok]).correlation
                rho_a = (spearmanr(y[ok] - base[ok], p[ok] - base[ok]).correlation
                         if nm == "stRF" else np.nan)
                rows.append(dict(origin=int(origin), horizon=h, method=nm,
                                 n=int(ok.sum()), rho_level=rho,
                                 rho_anomaly=rho_a,
                                 rmse=float(np.sqrt(np.mean((y[ok]-p[ok])**2)))))
    return pd.DataFrame(rows)


# ============================================ F/G. intervals + AOA
def aoa_di(Xtr, Xq, imp):
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-12
    w = np.sqrt(np.maximum(imp, 0))
    Ztr = (Xtr - mu) / sd * w
    Zq = (Xq - mu) / sd * w
    D = cdist(Ztr, Ztr)
    dbar = D[np.triu_indices(len(Ztr), 1)].mean()
    np.fill_diagonal(D, np.inf)
    di_tr = D.min(1) / dbar
    return cdist(Zq, Ztr).min(1) / dbar, float(np.percentile(di_tr, 95))


def qrf_interval(m, Xtr, ytr, Xq, lo=0.10, hi=0.90):
    """Quantile regression forest via leaf membership (Meinshausen 2006)."""
    leaves_tr = m.apply(Xtr)          # (ntr, ntrees)
    leaves_q = m.apply(Xq)
    n_tree = leaves_tr.shape[1]
    L, U = np.zeros(len(Xq)), np.zeros(len(Xq))
    for i in range(len(Xq)):
        w = np.zeros(len(Xtr))
        for t in range(n_tree):
            same = leaves_tr[:, t] == leaves_q[i, t]
            s = same.sum()
            if s:
                w[same] += 1.0 / s
        w /= w.sum()
        o = np.argsort(ytr)
        cw = np.cumsum(w[o])
        L[i] = ytr[o][np.searchsorted(cw, lo)]
        U[i] = ytr[o][min(np.searchsorted(cw, hi), len(ytr) - 1)]
    return np.expm1(L), np.expm1(U)


def interval_study(agg, alpha=0.80):
    rows = []
    regimes = [("interpolation (random 10-fold)", "kf"),
               ("extrapolation (LORO)", "loro")]
    for label, kind in regimes:
        if kind == "kf":
            splits = list(KFold(n_splits=5, shuffle=True,
                                random_state=2).split(agg))
            splits = [(agg.iloc[a], agg.iloc[b]) for a, b in splits]
        else:
            splits = [(agg[agg.region != r], agg[agg.region == r])
                      for r in REGIONS]
        for tr, te in splits:
            if len(te) < 10:
                continue
            # calibration split for conformal
            idx = np.arange(len(tr))
            rs = np.random.default_rng(3).permutation(idx)
            cal_i, fit_i = rs[:max(30, len(idx)//4)], rs[max(30, len(idx)//4):]
            tr_fit, tr_cal = tr.iloc[fit_i], tr.iloc[cal_i]

            m, Xtr, Xq, ytr, _ = fit_full(tr, te)
            per_tree = np.stack([t.predict(Xq) for t in m.estimators_])
            lo_e = np.expm1(np.percentile(per_tree, 10, axis=0))
            hi_e = np.expm1(np.percentile(per_tree, 90, axis=0))

            lo_q, hi_q = qrf_interval(m, Xtr, ytr, Xq)

            m2, Xtr2, Xcal, ytr2, _ = fit_full(tr_fit, tr_cal)
            resid = np.abs(np.log1p(tr_cal["dens"].to_numpy()) - m2.predict(Xcal))
            q = np.quantile(resid, alpha)
            _, _, Xq2, _, _ = fit_full(tr_fit, te)
            ctr = m2.predict(Xq2)
            lo_c, hi_c = np.expm1(ctr - q), np.expm1(ctr + q)

            di, thr = aoa_di(Xtr, Xq, m.feature_importances_)
            inside = di <= thr
            y = te["dens"].to_numpy()
            for nm, lo_v, hi_v in [("ensemble spread", lo_e, hi_e),
                                   ("QRF", lo_q, hi_q),
                                   ("split conformal", lo_c, hi_c)]:
                cov = (y >= lo_v) & (y <= hi_v)
                rows.append(dict(
                    regime=label, interval=nm,
                    coverage=float(cov.mean()),
                    coverage_inside=float(cov[inside].mean()) if inside.any() else np.nan,
                    coverage_outside=float(cov[~inside].mean()) if (~inside).any() else np.nan,
                    width=float(np.mean(hi_v - lo_v)),
                    width_inside=float(np.mean((hi_v-lo_v)[inside])) if inside.any() else np.nan,
                    width_outside=float(np.mean((hi_v-lo_v)[~inside])) if (~inside).any() else np.nan,
                    pct_outside=float((~inside).mean())))
    return pd.DataFrame(rows)


# ======================================================== H. leakage
def leakage_check(agg):
    tr = agg[agg.region != "CSMC"].copy()
    te = agg[agg.region == "CSMC"].copy()
    _, b_q, _ = build_blocks(tr, te)
    te2 = te.copy()
    rng = np.random.default_rng(0)
    te2["dens"] = rng.permutation(te2["dens"].to_numpy()) * 3.0 + 11.0
    _, b_q2, _ = build_blocks(tr, te2)
    same = all(np.allclose(b_q[b], b_q2[b]) for b in BLOCKS)
    # non-vacuity: perturbing a TRAINING response must move the query lags
    tr2 = tr.copy(); tr2.iloc[0, tr2.columns.get_loc("dens")] += 50.0
    _, b_q3, _ = build_blocks(tr2, te)
    moved = not np.allclose(b_q["LAG"], b_q3["LAG"])
    return dict(no_leakage=bool(same), non_vacuous=bool(moved))


if __name__ == "__main__":
    raw, agg = load()

    print("=" * 78)
    print("D. VALIDATION REGIMES ON ONE SCALE")
    print("=" * 78)
    cv = cv_regimes(agg)
    cv.to_csv(f"{RES}/case_cv_regimes.csv", index=False)
    print(cv.round(3).to_string(index=False))
    out["cv_regimes"] = cv.round(4).to_dict("records")

    print()
    print("=" * 78)
    print("E. ROLLING ORIGIN FORECAST SKILL")
    print("=" * 78)
    ro = rolling_origin(agg)
    ro.to_csv(f"{RES}/case_rolling_origin.csv", index=False)
    lv = ro.pivot_table(index="horizon", columns="method", values="rho_level")
    print("  Rank correlation with observed density, by horizon:")
    print(lv.round(3).to_string())
    an = ro[ro.method == "stRF"].groupby("horizon")["rho_anomaly"].mean()
    print("\n  stRF skill at predicting departures from a site's mean:")
    print(an.round(3).to_string())
    out["rolling_origin"] = dict(levels=lv.round(4).to_dict(),
                                 anomaly=an.round(4).to_dict())

    print()
    print("=" * 78)
    print("F/G. PREDICTION INTERVALS (nominal 80%)")
    print("=" * 78)
    iv = interval_study(agg)
    iv.to_csv(f"{RES}/case_intervals.csv", index=False)
    agg_iv = iv.groupby(["regime", "interval"])[
        ["coverage", "width", "width_inside", "width_outside",
         "coverage_inside", "coverage_outside", "pct_outside"]].mean()
    print(agg_iv.round(3).to_string())
    out["intervals"] = agg_iv.round(4).reset_index().to_dict("records")

    print()
    print("=" * 78)
    print("H. LEAKAGE ASSERTION (case study feature builder)")
    print("=" * 78)
    lk = leakage_check(agg)
    print(f"  held-out responses cannot alter any predictor : {lk['no_leakage']}")
    print(f"  test is non-vacuous                           : {lk['non_vacuous']}")
    out["leakage"] = lk

    with open(f"{RES}/case_corrected2.json", "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\nWROTE case_cv_regimes.csv, case_rolling_origin.csv, "
          f"case_intervals.csv, case_corrected2.json")
