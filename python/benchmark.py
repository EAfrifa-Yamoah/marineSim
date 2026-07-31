"""
benchmark.py
============
Full factorial simulation benchmark for the stRF paper.

Grid: 6 worlds x 5 sample sizes x 3 designs x 2 detection x 3 timesteps x R reps,
fitting GLM, GAM, RF, SpatialRF, stRF on every draw. Results are written one CSV
per world to results/ so the run is restartable and chunkable.

Usage:
    python3 benchmark.py <world_start> <world_end> [reps]
e.g. python3 benchmark.py 0 2      # worlds 0,1 (end exclusive), default reps
     python3 benchmark.py 0 6 6    # all worlds, 6 reps
"""
import sys, time, os
import numpy as np
import pandas as pd
from sim_engine import World, fit_eval, make_worlds

OUT = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(OUT, exist_ok=True)

SAMPLE_SIZES = [30, 50, 100, 200, 500]
DESIGNS = ["random", "clustered", "stratified"]
DETECTION = [1.0]            # main grid; detection=0.7 robustness run separately
TIMESTEPS = [1, 3, 5]
METHODS = ["GLM", "RF", "SpatialRF", "stRF"]   # GAM/GEO run on a reduced grid
BASE_SEED = 20260601


def world_specs():
    """Mirror make_worlds() ordering so world ids/seeds are reproducible."""
    specs = []
    wid = 0
    for stat in (True, False):
        for rng_km in (30, 80, 200):
            name = f"{'stat' if stat else 'nonstat'}_r{rng_km}"
            specs.append((wid, name, stat, rng_km, BASE_SEED + wid))
            wid += 1
    return specs


def run_world(wid, name, stat, rng_km, seed, reps):
    # build the ground-truth world once (n_time=5 is the max we evaluate)
    world = World(stationary=stat, range_km=rng_km, n_time=5, seed=seed)
    rows = []
    t0 = time.time()
    nfit = 0
    for n in SAMPLE_SIZES:
        for design in DESIGNS:
            for det in DETECTION:
                for nt in TIMESTEPS:
                    for rep in range(reps):
                        # deterministic per-cell seed
                        s = (seed * 7919 + n * 131 + DESIGNS.index(design) * 17
                             + int(det * 10) * 53 + nt * 7 + rep)
                        rng = np.random.default_rng(s)
                        r = fit_eval(world, n=n, design=design, detection=det,
                                     n_time=nt, rng=rng, methods=METHODS)
                        nfit += 1
                        if r is None:
                            continue
                        rows.append(dict(
                            world=name, world_id=wid, stationary=stat,
                            range_km=rng_km, n=n, design=design,
                            detection=det, n_time=nt, rep=rep,
                            **{m: r.get(m, np.nan) for m in METHODS}))
    df = pd.DataFrame(rows)
    fp = os.path.join(OUT, f"bench_world{wid}_{name}.csv")
    df.to_csv(fp, index=False)
    print(f"[world {wid} {name}] {len(df)} rows, {nfit} fits, "
          f"{time.time()-t0:.1f}s -> {os.path.basename(fp)}", flush=True)
    return df


if __name__ == "__main__":
    w0 = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    w1 = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    reps = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    specs = world_specs()
    T0 = time.time()
    for (wid, name, stat, rng_km, seed) in specs[w0:w1]:
        run_world(wid, name, stat, rng_km, seed, reps)
    print(f"TOTAL worlds {w0}:{w1} in {time.time()-T0:.1f}s", flush=True)
