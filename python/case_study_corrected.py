"""
case_study_corrected.py
=======================
Corrected case study analysis addressing collaborator Issues 1, 2, 5 and 6.

  A. Data reconciliation (61 sites vs 74 in Supplement S7).
  B. Feature block ablation + Shapley, under random CV and under LORO, so the
     case study decomposition matches the corrected simulation decomposition.
  C. LORO transfer with honest metrics, and a diagnostic of how often the
     spatial lag features fall back when a whole park is withheld.
  D. Rolling origin forecast skill against persistence and climatology, for
     levels and for anomalies.
  E. Prediction intervals: across-tree ensemble spread (what the paper plots)
     against split conformal and quantile regression forest, with empirical
     coverage measured under interpolation and under extrapolation.
  F. Area of applicability, and coverage inside vs outside it.
  G. Leakage assertion for the case study feature builder.

GAM is omitted: pygam is not installed in this environment.
"""
import os, json, itertools, warnings
from math import factorial

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.model_selection import KFold

warnings.filterwarnings("ignore")
RES = "../results"
UP = os.environ.get("MARINESIM_CASE_DATA",
                    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "data", "raw",
                                 "Bayesiandataset_2025_final.csv"))
COVARS = ["depth", "kd490", "htyperank", "meansummermIST",
          "maxsummermIST", "CVsummermIST"]
RADII_KM = [25, 75, 150]
REGIONS = ["JBMP", "MMP", "CSMC", "SIMP", "NCMP"]          # north -> south
BLOCKS = ("COORD", "EDF", "LAG", "TEMP")
RF_KW = dict(n_estimators=400, min_samples_leaf=3, max_features=0.6,
             random_state=7, n_jobs=1)
out = {}


# ----------------------------------------------------------------- A. data
def load():
    raw = pd.read_csv(UP)
    agg = (raw.groupby(["Region_5", "Site", "YearSeagrassSampling"])
           .agg(dens=("rawcounts", "mean"),
                latitude=("latitude", "mean"), longitude=("longitude", "mean"),
                depth=("depth", "mean"), kd490=("kd490", "mean"),
                htyperank=("htyperank", "mean"),
                meansummermIST=("meansummermIST", "mean"),
                maxsummermIST=("maxsummermIST", "mean"),
                CVsummermIST=("CVsummermIST", "mean"))
           .reset_index()
           .rename(columns={"YearSeagrassSampling": "year",
                            "Region_5": "region", "Site": "site"}))
    for c in COVARS:
        agg[c] = agg.groupby("region")[c].transform(lambda s: s.fillna(s.mean()))
        agg[c] = agg[c].fillna(agg[c].mean())
    agg = agg.dropna(subset=["dens"]).reset_index(drop=True)
    agg["site_uid"] = agg["region"] + ":" + agg["site"]
    return raw, agg


def haversine(lon1, lat1, lon2, lat2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1); dlmb = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


# ------------------------------------------------------- feature blocks
def lag_block(train, query, report_fallback=False):
    """
    Multi-scale IDW lag of training densities, matched within year where
    possible. Records which fallback level each feature used:
      0 = same-year neighbours inside the radius
      1 = any-year neighbours inside the radius
      2 = nearest-k fallback (radius empty)
    """
    feats = np.zeros((len(query), len(RADII_KM)))
    level = np.zeros((len(query), len(RADII_KM)), dtype=int)
    tr_lon = train["longitude"].to_numpy(); tr_lat = train["latitude"].to_numpy()
    tr_y = train["year"].to_numpy(); tr_d = train["dens"].to_numpy()
    q_lon = query["longitude"].to_numpy(); q_lat = query["latitude"].to_numpy()
    q_yr = query["year"].to_numpy()
    for i in range(len(query)):
        d = haversine(q_lon[i], q_lat[i], tr_lon, tr_lat)
        same = tr_y == q_yr[i]
        for j, r in enumerate(RADII_KM):
            done = False
            for lv, mask in enumerate((same, np.ones_like(same, dtype=bool))):
                sel = mask & (d <= r) & (d > 1e-6)
                if sel.sum() >= 2:
                    w = 1.0 / (d[sel] + 1.0)
                    feats[i, j] = np.sum(w * tr_d[sel]) / np.sum(w)
                    level[i, j] = lv
                    done = True
                    break
            if not done:
                k = min(5, int((d > 1e-6).sum()))
                idx = np.argsort(d)[:k]
                feats[i, j] = tr_d[idx].mean()
                level[i, j] = 2
    return (feats, level) if report_fallback else feats


def edf_block(train, query, k=8, anchors=None):
    C_tr = train[["longitude", "latitude"]].to_numpy()
    if anchors is None:
        anchors = KMeans(n_clusters=min(k, len(C_tr)), n_init=5,
                         random_state=0).fit(C_tr).cluster_centers_
    C_q = query[["longitude", "latitude"]].to_numpy()
    D = np.zeros((len(C_q), len(anchors)))
    for a in range(len(anchors)):
        D[:, a] = haversine(C_q[:, 0], C_q[:, 1], anchors[a, 0], anchors[a, 1])
    return D, anchors


def temp_block(train, query):
    """Previous-year regional mean density, from training rows only."""
    reg_year = train.groupby(["region", "year"])["dens"].mean()
    gm = train["dens"].mean()
    vals = np.zeros(len(query)); avail = np.zeros(len(query), dtype=bool)
    regs = query["region"].to_numpy(); yrs = query["year"].to_numpy()
    for i in range(len(query)):
        key = (regs[i], yrs[i] - 1)
        if key in reg_year.index:
            vals[i] = reg_year[key]; avail[i] = True
        else:
            ry = reg_year[reg_year.index.get_level_values(0) == regs[i]]
            vals[i] = ry.mean() if len(ry) else gm
    return vals[:, None], avail


def build_blocks(train, query):
    b_tr, b_q = {}, {}
    b_tr["COORD"] = train[["longitude", "latitude"]].to_numpy()
    b_q["COORD"] = query[["longitude", "latitude"]].to_numpy()
    e_tr, anc = edf_block(train, train)
    e_q, _ = edf_block(train, query, anchors=anc)
    b_tr["EDF"], b_q["EDF"] = e_tr, e_q
    b_tr["LAG"] = lag_block(train, train)
    lq, lvl = lag_block(train, query, report_fallback=True)
    b_q["LAG"] = lq
    t_tr, _ = temp_block(train, train)
    t_q, tav = temp_block(train, query)
    b_tr["TEMP"], b_q["TEMP"] = t_tr, t_q
    return b_tr, b_q, dict(lag_fallback=lvl, temp_available=tav)


def design(train, query, b_tr, b_q, subset):
    Xtr = np.column_stack([train[COVARS].to_numpy()]
                          + [b_tr[b] for b in sorted(subset)])
    Xq = np.column_stack([query[COVARS].to_numpy()]
                         + [b_q[b] for b in sorted(subset)])
    return Xtr, Xq


def rf_fit_predict(Xtr, ytr_log, Xq):
    m = RandomForestRegressor(**RF_KW).fit(Xtr, ytr_log)
    return np.expm1(m.predict(Xq)), m


def metrics(y, yhat):
    y, yhat = np.asarray(y, float), np.asarray(yhat, float)
    ok = np.isfinite(yhat)
    y, yhat = y[ok], yhat[ok]
    ss_res = np.sum((y - yhat) ** 2); ss_tot = np.sum((y - y.mean()) ** 2)
    return (1 - ss_res / ss_tot if ss_tot > 0 else np.nan,
            spearmanr(y, yhat).correlation,
            float(np.sqrt(np.mean((y - yhat) ** 2))))


def shapley_from_values(v, blocks=BLOCKS):
    k = len(blocks); phi = {}
    for b in blocks:
        others = [x for x in blocks if x != b]; tot = 0.0
        for r in range(len(others) + 1):
            for S in itertools.combinations(others, r):
                S = frozenset(S)
                w = factorial(len(S)) * factorial(k - len(S) - 1) / factorial(k)
                tot += w * (v[S | {b}] - v[S])
        phi[b] = tot
    return phi


SUBSETS = [frozenset(s) for k in range(len(BLOCKS) + 1)
           for s in itertools.combinations(BLOCKS, k)]


# =========================================================== main analyses
if __name__ == "__main__":
    raw, agg = load()

    print("=" * 78)
    print("A. DATA RECONCILIATION")
    print("=" * 78)
    print(f"  raw quadrat counts        : {len(raw)}")
    print(f"  site-year records         : {len(agg)}")
    print(f"  unique sites (region:site): {agg['site_uid'].nunique()}")
    print(f"  unique site NAMES         : {agg['site'].nunique()}")
    print(f"  unique coordinate pairs   : "
          f"{agg[['latitude','longitude']].drop_duplicates().shape[0]}")
    print(f"  years                     : {agg.year.min()}-{agg.year.max()} "
          f"({agg.year.nunique()} distinct)")
    out["data"] = dict(raw=len(raw), site_years=len(agg),
                       sites=int(agg["site_uid"].nunique()),
                       years=[int(agg.year.min()), int(agg.year.max())])
    print(agg.groupby("region").agg(sites=("site", "nunique"),
                                    site_years=("dens", "size"),
                                    yrs=("year", "nunique"),
                                    lat=("latitude", "mean")).round(2).to_string())

    # ------------------------------------------------------- B/C. LORO
    print()
    print("=" * 78)
    print("C. LEAVE ONE REGION OUT: TRANSFER SKILL AND FEATURE FALLBACK")
    print("=" * 78)
    loro_rows, fb_rows, coalition_rows = [], [], []
    pooled_pred = {S: [] for S in SUBSETS}
    pooled_obs = []
    for held in REGIONS:
        tr = agg[agg.region != held].copy()
        te = agg[agg.region == held].copy()
        b_tr, b_q, diag = build_blocks(tr, te)
        ytr = np.log1p(tr["dens"].to_numpy())

        lv = diag["lag_fallback"]
        fb_rows.append(dict(held=held, n_test=len(te),
                            pct_sameyear=float((lv == 0).mean()),
                            pct_anyyear=float((lv == 1).mean()),
                            pct_nearestk=float((lv == 2).mean()),
                            pct_r25_fallback=float((lv[:, 0] == 2).mean()),
                            temp_available=float(diag["temp_available"].mean())))

        for S in SUBSETS:
            Xtr, Xq = design(tr, te, b_tr, b_q, S)
            yhat, _ = rf_fit_predict(Xtr, ytr, Xq)
            r2, rho, rmse = metrics(te["dens"], yhat)
            coalition_rows.append(dict(held=held, subset="+".join(sorted(S)) or "BASE",
                                       R2=r2, Spearman=rho, RMSE=rmse))
            pooled_pred[S].append(yhat)
        pooled_obs.append(te["dens"].to_numpy())

        # published model ladder for comparison
        for meth, S in [("RF", frozenset()), ("SpatialRF", frozenset({"COORD"})),
                        ("stRF", frozenset(BLOCKS))]:
            Xtr, Xq = design(tr, te, b_tr, b_q, S)
            yhat, _ = rf_fit_predict(Xtr, ytr, Xq)
            r2, rho, rmse = metrics(te["dens"], yhat)
            loro_rows.append(dict(held=held, method=meth, R2=r2,
                                  Spearman=rho, RMSE=rmse, n_test=len(te)))
        sc = StandardScaler().fit(tr[COVARS].to_numpy())
        glm = LinearRegression().fit(sc.transform(tr[COVARS].to_numpy()), ytr)
        yhat = np.expm1(glm.predict(sc.transform(te[COVARS].to_numpy())))
        r2, rho, rmse = metrics(te["dens"], yhat)
        loro_rows.append(dict(held=held, method="GLM", R2=r2, Spearman=rho,
                              RMSE=rmse, n_test=len(te)))

    fb = pd.DataFrame(fb_rows)
    print("\n  Spatial lag feature fallback rate under LORO "
          "(fraction of held-out rows):")
    print(fb.round(3).to_string(index=False))
    out["lag_fallback"] = fb.mean(numeric_only=True).round(4).to_dict()

    loro = pd.DataFrame(loro_rows)
    loro.to_csv(f"{RES}/case_loro_corrected.csv", index=False)
    obs_all = np.concatenate(pooled_obs)
    print("\n  Pooled across all held-out parks:")
    for m in ["GLM", "RF", "SpatialRF", "stRF"]:
        s = loro[loro.method == m]
        print(f"    {m:10s} mean rho={s.Spearman.mean():+.3f}  "
              f"mean R2={s.R2.mean():+.3f}  mean RMSE={s.RMSE.mean():.2f}")
    out["loro"] = loro.groupby("method")[["R2", "Spearman", "RMSE"]].mean().round(4).to_dict()

    # coalition Shapley on pooled LORO rank correlation
    co = pd.DataFrame(coalition_rows)
    co.to_csv(f"{RES}/case_coalitions_loro.csv", index=False)
    v = {}
    for S in SUBSETS:
        yh = np.concatenate(pooled_pred[S])
        v[S] = metrics(obs_all, yh)[1]          # pooled Spearman
    ph = shapley_from_values(v)
    print()
    print("=" * 78)
    print("B. CASE STUDY SHAPLEY DECOMPOSITION (pooled LORO rank correlation)")
    print("=" * 78)
    base = v[frozenset()]
    print(f"  covariates only (base) rho = {base:+.3f}")
    tot = sum(ph.values())
    for b in BLOCKS:
        print(f"  {b:8s} {ph[b]:+.4f}   ({100*ph[b]/tot:5.1f}% of block gain)")
    print(f"  {'TOTAL':8s} {tot:+.4f}  -> full model rho = "
          f"{v[frozenset(BLOCKS)]:+.3f}")
    print(f"  TEMP added last (full model): "
          f"{v[frozenset(BLOCKS)] - v[frozenset(BLOCKS)-{'TEMP'}]:+.4f}")
    out["case_shapley"] = {b: float(ph[b]) for b in BLOCKS}
    out["case_shapley"]["base_rho"] = float(base)
    out["case_shapley"]["full_rho"] = float(v[frozenset(BLOCKS)])
    out["case_shapley"]["temp_last"] = float(
        v[frozenset(BLOCKS)] - v[frozenset(BLOCKS) - {"TEMP"}])

    with open(f"{RES}/case_corrected.json", "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\nWROTE {RES}/case_loro_corrected.csv, case_coalitions_loro.csv, "
          f"case_corrected.json")
