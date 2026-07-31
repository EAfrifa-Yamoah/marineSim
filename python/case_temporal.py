"""
case_temporal.py
================
Spatiotemporal evolution of the stRF case study.

Three experiments, all using the corrected leaf raw-scale predictor:

  E1  Hindcast reconstruction. Site-grouped 5-fold CV: every site's whole trajectory is
      predicted while that site is held out, then aggregated to regional means. Shows the
      model reconstructs the divergent regional dynamics for unseen sites.

  E2  Rolling-origin forecast. For each origin year t0, train on years <= t0 and forecast
      each later year (horizon h = 1..H). The year's environmental covariates are treated
      as known forcing; the density memory is restricted to <= t0. Compared against a
      persistence baseline (a site's last observed value) and a climatology baseline (a
      site's mean up to t0).

  E3  Blocked per-year CV. Hold out each year in turn (train on all other years). This is
      temporal interpolation - it may use future context - and forms the horizon-0
      reference (an upper bound) on the skill curve.

  E4  Full reconstruction field for the space-time (site x year) diagram and the annual
      corridor surfaces.

Outputs: ../results/temporal_*.csv
"""
import io, contextlib, json
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold
from scipy.stats import spearmanr

with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
    import case_study as cs

RES = "../results"
agg = cs.agg.copy().reset_index(drop=True)
YEARS = np.arange(int(agg["year"].min()), int(agg["year"].max()) + 1)
RF_KW = dict(n_estimators=600, min_samples_leaf=3, max_features=0.6, random_state=7, n_jobs=1)


def strf_leafraw(train, query):
    """Fit stRF on log1p(density); predict the raw-scale conditional mean from the leaves
    (removes retransformation bias; matches the corrected Figure 7 surface)."""
    y = np.log1p(train["dens"].to_numpy())
    Xtr, Xq = cs.build_design(train, query, "stRF")
    rf = RandomForestRegressor(**RF_KW).fit(Xtr, y)
    raw = train["dens"].to_numpy()
    Ltr, Lq = rf.apply(Xtr), rf.apply(Xq)
    pred = np.zeros(len(query))
    for b in range(Ltr.shape[1]):
        lm = pd.Series(raw).groupby(Ltr[:, b]).mean()
        pred += lm.reindex(Lq[:, b]).to_numpy()
    return pred / Ltr.shape[1]


# ============================================================ E1 reconstruction
def e1_reconstruction():
    oof = np.full(len(agg), np.nan)
    for tr_i, te_i in GroupKFold(n_splits=5).split(agg, groups=agg["site"]):
        oof[te_i] = strf_leafraw(agg.iloc[tr_i], agg.iloc[te_i])
    d = agg[["region", "year", "site", "dens", "latitude"]].copy()
    d["pred"] = oof
    d.to_csv(f"{RES}/temporal_reconstruction.csv", index=False)
    # regional trajectories
    traj = (d.groupby(["region", "year"])
              .agg(obs=("dens", "mean"), pred=("pred", "mean"), n=("dens", "size"))
              .reset_index())
    traj.to_csv(f"{RES}/temporal_regional_traj.csv", index=False)
    rho = spearmanr(d["dens"], d["pred"]).correlation
    rmse = float(np.sqrt(np.mean((d["dens"] - d["pred"]) ** 2)))
    print(f"E1 reconstruction (site-grouped OOF): Spearman {rho:.3f}  RMSE {rmse:.2f}")
    return traj


# ============================================================ E2 rolling forecast
def e2_forecast(H=5):
    rows = []
    origins = [t for t in YEARS if (agg["year"] <= t).sum() >= 120 and t <= YEARS.max() - 1]
    for t0 in origins:
        train = agg[agg["year"] <= t0]
        # per-site memory available at the origin
        last = train.sort_values("year").groupby("site")["dens"].last()
        clim = train.groupby("site")["dens"].mean()
        for h in range(1, H + 1):
            Y = t0 + h
            if Y > YEARS.max():
                continue
            tgt = agg[(agg["year"] == Y) & (agg["site"].isin(train["site"].unique()))]
            if len(tgt) < 20:
                continue
            pred = strf_leafraw(train, tgt)
            for i, (_, q) in enumerate(tgt.reset_index(drop=True).iterrows()):
                rows.append(dict(origin=int(t0), horizon=h, year=int(Y),
                                 region=q["region"], site=q["site"], obs=q["dens"],
                                 stRF=pred[i],
                                 persistence=float(last.get(q["site"], np.nan)),
                                 climatology=float(clim.get(q["site"], np.nan))))
    f = pd.DataFrame(rows)
    f.to_csv(f"{RES}/temporal_forecast_raw.csv", index=False)
    # skill by horizon (pooled), all three predictors
    sk = []
    for h, sub in f.groupby("horizon"):
        for m in ["stRF", "persistence", "climatology"]:
            s = sub.dropna(subset=[m])
            sk.append(dict(horizon=int(h), method=m, n=len(s),
                           spearman=spearmanr(s["obs"], s[m]).correlation,
                           rmse=float(np.sqrt(np.mean((s["obs"] - s[m]) ** 2))),
                           mae=float(np.mean(np.abs(s["obs"] - s[m])))))
    print(f"E2 forecast: {len(origins)} origins {origins[0]}-{origins[-1]}, {len(f)} forecasts")
    return pd.DataFrame(sk)


# ============================================================ E3 blocked per-year
def e3_blocked():
    rows = []
    for Y in YEARS:
        tgt = agg[agg["year"] == Y]
        if len(tgt) < 20:
            continue
        train = agg[agg["year"] != Y]
        pred = strf_leafraw(train, tgt)
        last = agg[agg["year"] < Y].sort_values("year").groupby("site")["dens"].last()
        clim = agg[agg["year"] != Y].groupby("site")["dens"].mean()
        t = tgt.reset_index(drop=True)
        for i in range(len(t)):
            rows.append(dict(year=int(Y), region=t["region"].iloc[i], obs=t["dens"].iloc[i],
                             stRF=pred[i],
                             persistence=float(last.get(t["site"].iloc[i], np.nan)),
                             climatology=float(clim.get(t["site"].iloc[i], np.nan))))
    b = pd.DataFrame(rows)
    sk = []
    for m in ["stRF", "persistence", "climatology"]:
        s = b.dropna(subset=[m])
        sk.append(dict(horizon=0, method=m, n=len(s),
                       spearman=spearmanr(s["obs"], s[m]).correlation,
                       rmse=float(np.sqrt(np.mean((s["obs"] - s[m]) ** 2))),
                       mae=float(np.mean(np.abs(s["obs"] - s[m])))))
    print(f"E3 blocked per-year: {b['year'].nunique()} year-folds, {len(b)} predictions")
    return pd.DataFrame(sk)


# ============================================================ E4 reconstruction field
def e4_field():
    # fit on all data; predict every site at every year (regular annual grid)
    sites = agg[["site", "region", "latitude", "longitude"]].drop_duplicates("site")
    # static covariates per site (mean); time-varying handled below
    stat = agg.groupby("site")[["depth", "htyperank"]].mean()
    tvar = ["kd490", "meansummermIST", "maxsummermIST", "CVsummermIST"]
    # per (site, year) time-varying covariate: observed if present else nearest-year value
    obs_lookup = agg.set_index(["site", "year"])
    grid = []
    for _, s in sites.iterrows():
        sy = agg[agg["site"] == s["site"]].set_index("year")
        for Y in YEARS:
            row = dict(site=s["site"], region=s["region"], year=Y,
                       latitude=s["latitude"], longitude=s["longitude"],
                       depth=stat.loc[s["site"], "depth"], htyperank=stat.loc[s["site"], "htyperank"])
            for c in tvar:
                if Y in sy.index:
                    row[c] = sy.loc[Y, c]
                else:  # nearest observed year for this site
                    j = sy.index[np.argmin(np.abs(sy.index - Y))]
                    row[c] = sy.loc[j, c]
            grid.append(row)
    G = pd.DataFrame(grid)
    G["pred"] = strf_leafraw(agg, G)
    # attach observed where available
    G = G.merge(agg[["site", "year", "dens"]], on=["site", "year"], how="left")
    G.to_csv(f"{RES}/temporal_field.csv", index=False)
    print(f"E4 field: {sites.shape[0]} sites x {len(YEARS)} years = {len(G)} cells "
          f"({G['dens'].notna().sum()} observed)")
    return G


if __name__ == "__main__":
    traj = e1_reconstruction()
    sk2 = e2_forecast()
    sk3 = e3_blocked()
    skill = pd.concat([sk3, sk2], ignore_index=True).sort_values(["method", "horizon"])
    skill.to_csv(f"{RES}/temporal_skill.csv", index=False)
    e4_field()
    print("\nSKILL CURVE (Spearman by horizon):")
    print(skill.pivot(index="horizon", columns="method", values="spearman").round(3).to_string())
    print("\nSKILL CURVE (RMSE by horizon):")
    print(skill.pivot(index="horizon", columns="method", values="rmse").round(2).to_string())
    print("\nTEMPORAL ANALYSIS DONE")
