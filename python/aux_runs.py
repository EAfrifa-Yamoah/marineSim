"""
aux_runs.py
===========
Auxiliary benchmark runs that complement the main 4-method grid.

  gam   : GAM (pygam) AUC on the same cells as the main grid (det=1.0).
  geo   : Matern-GP geostatistical comparator (sdmTMB-class) paired with stRF
          on a reduced subset (random design, det=1.0, ts=5) + fit time + conv.
  det07 : detection = 0.7 robustness subset for the 4 fast methods.

Usage: python3 aux_runs.py <gam|geo|det07> <world_start> <world_end>
"""
import sys, os, time
import numpy as np
import pandas as pd
from sim_engine import World, fit_eval

OUT = os.path.join(os.path.dirname(__file__), "..", "results")
BASE_SEED = 20260601
SAMPLE_SIZES = [30, 50, 100, 200, 500]
DESIGNS = ["random", "clustered", "stratified"]
TIMESTEPS = [1, 3, 5]


def specs():
    out = []; wid = 0
    for stat in (True, False):
        for rng_km in (30, 80, 200):
            out.append((wid, f"{'stat' if stat else 'nonstat'}_r{rng_km}",
                        stat, rng_km, BASE_SEED + wid)); wid += 1
    return out


def cellseed(seed, n, design, det, nt, rep):
    return (seed * 7919 + n * 131 + DESIGNS.index(design) * 17
            + int(det * 10) * 53 + nt * 7 + rep)


def run_gam(world, name, wid, seed, reps=6):
    rows = []
    for n in SAMPLE_SIZES:
        for design in DESIGNS:
            for nt in TIMESTEPS:
                for rep in range(reps):
                    rng = np.random.default_rng(cellseed(seed, n, design, 1.0, nt, rep))
                    r = fit_eval(world, n=n, design=design, detection=1.0,
                                 n_time=nt, rng=rng, methods=["GAM"])
                    if r is None:
                        continue
                    rows.append(dict(world=name, world_id=wid, n=n, design=design,
                                     detection=1.0, n_time=nt, rep=rep,
                                     GAM=r.get("GAM", np.nan)))
    return pd.DataFrame(rows)


def run_geo(world, name, wid, seed, reps=4):
    rows = []
    for n in SAMPLE_SIZES:
        for rep in range(reps):
            rng = np.random.default_rng(cellseed(seed, n, "random", 1.0, 5, rep) + 999)
            r = fit_eval(world, n=n, design="random", detection=1.0, n_time=5,
                         rng=rng, methods=["stRF"], include_geo=True)
            if r is None:
                continue
            rows.append(dict(world=name, world_id=wid, n=n, rep=rep,
                             stRF=r.get("stRF", np.nan), GEO=r.get("GEO", np.nan),
                             GEO_time=r.get("GEO_time", np.nan),
                             GEO_conv=r.get("GEO_conv", 0)))
    return pd.DataFrame(rows)


def run_det07(world, name, wid, seed, reps=6):
    rows = []
    methods = ["GLM", "RF", "SpatialRF", "stRF"]
    for n in [50, 200, 500]:
        for design in ["random", "clustered"]:
            for nt in [1, 5]:
                for rep in range(reps):
                    rng = np.random.default_rng(cellseed(seed, n, design, 0.7, nt, rep))
                    r = fit_eval(world, n=n, design=design, detection=0.7,
                                 n_time=nt, rng=rng, methods=methods)
                    if r is None:
                        continue
                    rows.append(dict(world=name, world_id=wid, n=n, design=design,
                                     detection=0.7, n_time=nt, rep=rep,
                                     **{m: r.get(m, np.nan) for m in methods}))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    mode = sys.argv[1]
    w0 = int(sys.argv[2]); w1 = int(sys.argv[3])
    fn = {"gam": run_gam, "geo": run_geo, "det07": run_det07}[mode]
    for (wid, name, stat, rng_km, seed) in specs()[w0:w1]:
        t0 = time.time()
        world = World(stationary=stat, range_km=rng_km, n_time=5, seed=seed)
        df = fn(world, name, wid, seed)
        fp = os.path.join(OUT, f"aux_{mode}_world{wid}.csv")
        df.to_csv(fp, index=False)
        print(f"[{mode} world {wid} {name}] {len(df)} rows {time.time()-t0:.1f}s",
              flush=True)
