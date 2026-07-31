"""
ablation.py
===========
Corrected feature-block decomposition for the stRF benchmark.

Addresses Major Issue 1. The published contrast

    d_temporal := AUC(stRF) - AUC(SpatialRF)

changes four feature blocks at once, because in sim_engine.fit_eval:

    SpatialRF = covariates + COORD
    stRF      = covariates + LAG + TEMP + EDF          (no COORD)

so the contrast ADDS {LAG, TEMP, EDF} and REMOVES {COORD}. It cannot be
read as the temporal contribution.

This script fits a common RF backbone over all 2^5 subsets of

    COORD  x, y
    EDF    distances to 15 k-means anchors
    LAG    multi-scale spatial lags at 25 / 75 / 150 km
    TEMP   y(nearest training site, t-1) and the site history mean
    NN0    CONTROL: same-timestep nearest-neighbour spatial lag

NN0 exists because TEMP is evaluated at query points by a nearest-training-
site lookup (sim_engine.temporal_lag). Under a temporally persistent site
effect, y(nearest site, t-1) is a close proxy for y(nearest site, t), so
TEMP can act as a fine-scale SPATIAL neighbour feature rather than as
temporal borrowing. Comparing TEMP's contribution with and without NN0 in
the coalition separates the two.

Covariates are always present. Detection = 1.0. Random design.
"""

import os, sys, time, itertools
from math import factorial

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim_engine import (World, sample_sites, spatial_lags, temporal_lag, edf,
                        GRID, RADII_KM, N_ANCHORS, RF_TREES, N_TEST)

BLOCKS = ("COORD", "EDF", "LAG", "TEMP", "NN0")
CORE_BLOCKS = ("COORD", "EDF", "LAG", "TEMP")   # the four in the published model
BASE_SEED = 20260601


# --------------------------------------------------------------------------- #
def nn_lag_same_timestep(train_xy, train_y, query_xy, is_train):
    """
    CONTROL block. Response at the nearest OTHER training site at the same
    timestep. Self-excluded at training points. This is the same-time
    analogue of TEMP's nearest-site lookup.
    """
    D = cdist(query_xy, train_xy)
    if is_train:
        np.fill_diagonal(D, np.inf)
    j = np.argmin(D, axis=1)
    return train_y[j]


def build_dataset(world, n, n_time, rep, design="random", seed=0):
    """Draw one replicate and construct every feature block for train and test."""
    rng = np.random.default_rng(seed)
    t_eval = n_time - 1

    site_idx = sample_sites(n, design, rng)
    site_xy = world.coords(site_idx)

    obs_y = {t: world.y[t].ravel()[site_idx].astype(float) for t in range(n_time)}
    y_tr = obs_y[t_eval]
    if y_tr.sum() == 0 or y_tr.sum() == len(y_tr):
        return None

    pool = np.setdiff1d(np.arange(GRID * GRID), site_idx)
    test_idx = rng.choice(pool, size=min(N_TEST, len(pool)), replace=False)
    test_xy = world.coords(test_idx)
    y_te = world.y[t_eval].ravel()[test_idx]
    p_te = world.p[t_eval].ravel()[test_idx]          # TRUE probability
    if y_te.sum() == 0 or y_te.sum() == len(y_te):
        return None

    Xc_tr = world.covariates(site_idx, t_eval)
    Xc_te = world.covariates(test_idx, t_eval)
    gmean = float(y_tr.mean())

    blocks_tr, blocks_te = {}, {}

    blocks_tr["COORD"] = site_xy
    blocks_te["COORD"] = test_xy

    k = min(N_ANCHORS, len(site_xy))
    anchors = KMeans(n_clusters=k, n_init=3, random_state=0).fit(site_xy).cluster_centers_
    blocks_tr["EDF"] = edf(site_xy, site_xy, anchors)
    blocks_te["EDF"] = edf(site_xy, test_xy, anchors)

    blocks_tr["LAG"] = spatial_lags(site_xy, y_tr, site_xy, RADII_KM, is_train=True)
    blocks_te["LAG"] = spatial_lags(site_xy, y_tr, test_xy, RADII_KM, is_train=False)

    if t_eval >= 1:
        prev_y = obs_y[t_eval - 1]
        tl_tr = temporal_lag(site_xy, prev_y, site_xy, gmean)
        tl_te = temporal_lag(site_xy, prev_y, test_xy, gmean)
        hist = np.stack([obs_y[u] for u in range(t_eval)], axis=0).mean(0)
        tree = cKDTree(site_xy)
        hm_tr = hist[tree.query(site_xy, k=1)[1]]
        hm_te = hist[tree.query(test_xy, k=1)[1]]
    else:
        tl_tr = np.full(len(site_xy), gmean); tl_te = np.full(len(test_xy), gmean)
        hm_tr = np.full(len(site_xy), gmean); hm_te = np.full(len(test_xy), gmean)
    blocks_tr["TEMP"] = np.column_stack([tl_tr, hm_tr])
    blocks_te["TEMP"] = np.column_stack([tl_te, hm_te])

    blocks_tr["NN0"] = nn_lag_same_timestep(site_xy, y_tr, site_xy, True)[:, None]
    blocks_te["NN0"] = nn_lag_same_timestep(site_xy, y_tr, test_xy, False)[:, None]

    return dict(Xc_tr=Xc_tr, Xc_te=Xc_te, y_tr=y_tr, y_te=y_te, p_te=p_te,
                blocks_tr=blocks_tr, blocks_te=blocks_te,
                site_xy=site_xy, test_xy=test_xy, rng_seed=seed)


def rf_auc(ds, subset, seed=0, return_spread=False):
    """Fit the common RF backbone on covariates + subset; return test AUC."""
    Xtr = np.column_stack([ds["Xc_tr"]] + [ds["blocks_tr"][b] for b in subset])
    Xte = np.column_stack([ds["Xc_te"]] + [ds["blocks_te"][b] for b in subset])
    rf = RandomForestClassifier(n_estimators=RF_TREES, n_jobs=1, random_state=seed)
    rf.fit(Xtr, ds["y_tr"].astype(int))
    pr = rf.predict_proba(Xte)[:, 1]
    auc = roc_auc_score(ds["y_te"], pr)
    if return_spread:
        per_tree = np.stack([t.predict_proba(Xte)[:, 1] for t in rf.estimators_])
        return auc, pr, per_tree
    return auc


def glm_auc(ds):
    m = LogisticRegression(max_iter=400, C=1.0)
    m.fit(ds["Xc_tr"], ds["y_tr"].astype(int))
    return roc_auc_score(ds["y_te"], m.predict_proba(ds["Xc_te"])[:, 1])


# --------------------------------------------------------------------------- #
def shapley(values, blocks):
    """Shapley value per block given v(S) for every subset S of `blocks`."""
    k = len(blocks)
    phi = {}
    for b in blocks:
        others = [x for x in blocks if x != b]
        tot = 0.0
        for r in range(len(others) + 1):
            for S in itertools.combinations(others, r):
                S = frozenset(S)
                w = factorial(len(S)) * factorial(k - len(S) - 1) / factorial(k)
                tot += w * (values[S | {b}] - values[S])
        phi[b] = tot
    return phi


# --------------------------------------------------------------------------- #
SAMPLE_SIZES = [30, 50, 100, 200, 500]
TIMESTEPS = [1, 3, 5]


def world_grid(n_realisations):
    """6 world configurations x n_realisations independent field realisations."""
    out, wid = [], 0
    for stat in (True, False):
        for rng_km in (30, 80, 200):
            cfg = f"{'stat' if stat else 'nonstat'}_r{rng_km}"
            for r in range(n_realisations):
                out.append((wid, cfg, r, stat, rng_km,
                            BASE_SEED + 1000 * r + wid))
            wid += 1
    return out


def required_subsets():
    """
    16 subsets of the four published blocks (full Shapley), plus four
    coalitions containing the NN0 control needed to test whether TEMP is
    acting as a spatial nearest-neighbour feature.
    """
    subs = [frozenset(s) for k in range(len(CORE_BLOCKS) + 1)
            for s in itertools.combinations(CORE_BLOCKS, k)]
    subs += [frozenset({"NN0"}),
             frozenset({"NN0", "TEMP"}),
             frozenset({"COORD", "EDF", "LAG", "NN0"}),
             frozenset({"COORD", "EDF", "LAG", "NN0", "TEMP"})]
    return subs


def main(n_realisations=3, reps=1, w0=0, w1=None,
         out="../results/ablation_rows.csv"):
    subsets = required_subsets()
    rows, t0, nfit = [], time.time(), 0

    grid = world_grid(n_realisations)
    grid = grid[w0:(w1 if w1 is not None else len(grid))]
    for (wid, cfg, real, stat, rng_km, wseed) in grid:
        world = World(stationary=stat, range_km=rng_km, n_time=5, seed=wseed)
        for n in SAMPLE_SIZES:
            for nt in TIMESTEPS:
                for rep in range(reps):
                    dseed = wseed * 31 + n * 131 + nt * 7 + rep
                    ds = build_dataset(world, n, nt, rep, seed=dseed)
                    if ds is None:
                        continue
                    rec = dict(world_cfg=cfg, world_id=wid, realisation=real,
                               stationary=stat, range_km=rng_km,
                               n=n, n_time=nt, rep=rep,
                               GLM=glm_auc(ds))
                    for S in subsets:
                        key = "v_" + ("BASE" if not S else
                                      "+".join(sorted(S)))
                        rec[key] = rf_auc(ds, sorted(S), seed=int(dseed) % 10000)
                        nfit += 1
                    rows.append(rec)
        el = time.time() - t0
        print(f"[world {wid} {cfg} real{real}] rows={len(rows)} fits={nfit} "
              f"{el:.0f}s", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print(f"WROTE {out}: {len(df)} datasets, {nfit} RF fits, "
          f"{time.time()-t0:.0f}s")
    return df


if __name__ == "__main__":
    nr = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    w0 = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    w1 = int(sys.argv[3]) if len(sys.argv) > 3 else None
    tag = sys.argv[4] if len(sys.argv) > 4 else "all"
    main(nr, 1, w0, w1, out=f"../results/ablation_rows_{tag}.csv")
