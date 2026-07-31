"""
decompose.py
============
Turn ablation_rows.csv into the corrected decomposition.

Produces
  1. The published contrast  d_temporal = stRF - SpatialRF, recomputed on the
     same backbone, to show what it actually measures.
  2. Sequential (order dependent) increments along a nested ladder.
  3. Shapley values over the four published blocks (order free).
  4. The NN0 confound test: TEMP's marginal contribution with and without a
     same timestep nearest neighbour lag already in the model.
  5. Bootstrap intervals resampling INDEPENDENT FIELD REALISATIONS rather
     than the six designed world types.
"""

import itertools, json
from math import factorial

import numpy as np
import pandas as pd

CORE = ("COORD", "EDF", "LAG", "TEMP")
R = "../results"


def key(S):
    return "v_" + ("BASE" if not S else "+".join(sorted(S)))


def load():
    df = pd.read_csv(f"{R}/ablation_rows.csv")
    # unique identifier for an independent field realisation
    df["field"] = df["world_cfg"] + "_r" + df["realisation"].astype(str)
    return df


# ------------------------------------------------------------------ Shapley
def shapley_row(row, blocks=CORE):
    k = len(blocks)
    phi = {}
    for b in blocks:
        others = [x for x in blocks if x != b]
        tot = 0.0
        for r in range(len(others) + 1):
            for S in itertools.combinations(others, r):
                S = frozenset(S)
                w = factorial(len(S)) * factorial(k - len(S) - 1) / factorial(k)
                tot += w * (row[key(S | {b})] - row[key(S)])
        phi[b] = tot
    return phi


def add_decompositions(df):
    # --- 1. the published contrast, recomputed on a common backbone ---------
    df["published_d_temporal"] = (df[key(frozenset({"LAG", "TEMP", "EDF"}))]
                                  - df[key(frozenset({"COORD"}))])
    df["published_d_spatial"] = (df[key(frozenset({"COORD"}))] - df["v_BASE"])
    df["published_d_ml"] = df["v_BASE"] - df["GLM"]
    df["published_d_total"] = (df[key(frozenset({"LAG", "TEMP", "EDF"}))]
                               - df["GLM"])

    # --- 2. corrected: TEMP added LAST to an otherwise complete model -------
    df["temporal_last"] = (df[key(frozenset({"COORD", "EDF", "LAG", "TEMP"}))]
                           - df[key(frozenset({"COORD", "EDF", "LAG"}))])
    # TEMP added FIRST to covariates only
    df["temporal_first"] = df[key(frozenset({"TEMP"}))] - df["v_BASE"]

    # --- 3. sequential ladder: covariates -> COORD -> EDF -> LAG -> TEMP ----
    df["seq_COORD"] = df[key(frozenset({"COORD"}))] - df["v_BASE"]
    df["seq_EDF"] = (df[key(frozenset({"COORD", "EDF"}))]
                     - df[key(frozenset({"COORD"}))])
    df["seq_LAG"] = (df[key(frozenset({"COORD", "EDF", "LAG"}))]
                     - df[key(frozenset({"COORD", "EDF"}))])
    df["seq_TEMP"] = df["temporal_last"]

    # --- 4. Shapley over the four published blocks --------------------------
    sh = df.apply(lambda r: shapley_row(r), axis=1, result_type="expand")
    sh.columns = [f"shap_{c}" for c in sh.columns]
    df = pd.concat([df, sh], axis=1)
    df["shap_ML"] = df["v_BASE"] - df["GLM"]

    # --- 5. NN0 confound test ----------------------------------------------
    # TEMP's marginal value when the model already has COORD+EDF+LAG
    df["temp_marg_noNN0"] = (df[key(frozenset({"COORD", "EDF", "LAG", "TEMP"}))]
                             - df[key(frozenset({"COORD", "EDF", "LAG"}))])
    # TEMP's marginal value when a same timestep 1-NN lag is ALSO present
    df["temp_marg_withNN0"] = (
        df[key(frozenset({"COORD", "EDF", "LAG", "NN0", "TEMP"}))]
        - df[key(frozenset({"COORD", "EDF", "LAG", "NN0"}))])
    # what NN0 alone buys over the same baseline (the spatial resolution effect)
    df["nn0_marg"] = (df[key(frozenset({"COORD", "EDF", "LAG", "NN0"}))]
                      - df[key(frozenset({"COORD", "EDF", "LAG"}))])
    return df


# --------------------------------------------------------------- bootstrap
def boot_ci(df, col, group="field", B=2000, seed=1):
    """Cluster bootstrap resampling independent field realisations."""
    rng = np.random.default_rng(seed)
    groups = df[group].unique()
    means = df.groupby(group)[col].mean()
    idx = rng.integers(0, len(groups), size=(B, len(groups)))
    draws = means.values[idx].mean(axis=1)
    return (float(df[col].mean()), float(np.percentile(draws, 2.5)),
            float(np.percentile(draws, 97.5)))


def boot_ci_worldtype(df, col, B=2000, seed=1):
    """The published approach: resample the six designed world types."""
    rng = np.random.default_rng(seed)
    means = df.groupby("world_cfg")[col].mean()
    idx = rng.integers(0, len(means), size=(B, len(means)))
    draws = means.values[idx].mean(axis=1)
    return (float(df[col].mean()), float(np.percentile(draws, 2.5)),
            float(np.percentile(draws, 97.5)))


def fmt(t):
    return f"{t[0]:+.4f}  [{t[1]:+.4f}, {t[2]:+.4f}]"


# -------------------------------------------------------------------- main
if __name__ == "__main__":
    df = add_decompositions(load())
    df.to_csv(f"{R}/ablation_decomposed.csv", index=False)
    out = {}

    print("=" * 78)
    print("A. WHAT THE PUBLISHED CONTRAST ACTUALLY MEASURES")
    print("=" * 78)
    print(f"  n datasets              : {len(df)}")
    print(f"  independent realisations: {df['field'].nunique()}")
    for lab, c in [("ML uplift (RF - GLM)", "published_d_ml"),
                   ("Spatial (RF+COORD - RF)", "published_d_spatial"),
                   ("'Temporal' (stRF - SpatialRF)", "published_d_temporal"),
                   ("Total (stRF - GLM)", "published_d_total")]:
        t = boot_ci(df, c)
        out[c] = t
        print(f"  {lab:32s} {fmt(t)}")

    print()
    print("=" * 78)
    print("B. CORRECTED TEMPORAL INCREMENT (only TEMP changes)")
    print("=" * 78)
    for lab, c in [("TEMP added last (full model)", "temporal_last"),
                   ("TEMP added first (covariates only)", "temporal_first"),
                   ("TEMP Shapley (order free)", "shap_TEMP")]:
        t = boot_ci(df, c)
        out[c] = t
        print(f"  {lab:36s} {fmt(t)}")

    print()
    print("=" * 78)
    print("C. SHAPLEY DECOMPOSITION OVER THE FOUR PUBLISHED BLOCKS")
    print("=" * 78)
    tot = 0.0
    shap_res = {}
    for b in ("COORD", "EDF", "LAG", "TEMP"):
        t = boot_ci(df, f"shap_{b}")
        shap_res[b] = t
        tot += t[0]
        print(f"  {b:8s} {fmt(t)}")
    tml = boot_ci(df, "shap_ML")
    shap_res["ML"] = tml
    print(f"  {'ML':8s} {fmt(tml)}")
    grand = tot + tml[0]
    print(f"  {'-'*60}")
    print(f"  total gain over GLM: {grand:+.4f}")
    for b in ("COORD", "EDF", "LAG", "TEMP"):
        print(f"    share {b:6s}: {100*shap_res[b][0]/grand:5.1f}%")
    print(f"    share {'ML':6s}: {100*tml[0]/grand:5.1f}%")
    print(f"    share spatial (COORD+EDF+LAG): "
          f"{100*sum(shap_res[b][0] for b in ('COORD','EDF','LAG'))/grand:5.1f}%")
    out["shapley"] = shap_res
    out["shapley_total"] = grand

    print()
    print("=" * 78)
    print("D. NN0 CONTROL: IS 'TEMPORAL' REALLY SPATIAL RESOLUTION?")
    print("=" * 78)
    for lab, c in [("TEMP marginal, no NN0 present", "temp_marg_noNN0"),
                   ("TEMP marginal, NN0 present", "temp_marg_withNN0"),
                   ("NN0 marginal (same timestep 1-NN)", "nn0_marg")]:
        t = boot_ci(df, c)
        out[c] = t
        print(f"  {lab:36s} {fmt(t)}")
    red = df["temp_marg_noNN0"].mean() - df["temp_marg_withNN0"].mean()
    print(f"  --> TEMP contribution absorbed by a same timestep 1-NN lag: "
          f"{red:+.4f} AUC "
          f"({100*red/max(df['temp_marg_noNN0'].mean(),1e-9):.0f}% of it)")

    print()
    print("=" * 78)
    print("E. TEMPORAL INCREMENT BY TEMPORAL DEPTH  (T = 1 is the placebo)")
    print("=" * 78)
    print(f"  {'T':>3} {'published (stRF-spRF)':>22} {'corrected (TEMP last)':>22}"
          f" {'TEMP | NN0':>12}")
    for T in sorted(df.n_time.unique()):
        s = df[df.n_time == T]
        print(f"  {T:>3} {s['published_d_temporal'].mean():>22.4f}"
              f" {s['temporal_last'].mean():>22.4f}"
              f" {s['temp_marg_withNN0'].mean():>12.4f}")

    print()
    print("=" * 78)
    print("F. TEMPORAL INCREMENT BY SAMPLE SIZE (tests the n >= 100 threshold)")
    print("=" * 78)
    print(f"  {'n':>4} {'published':>12} {'corrected':>12} {'corrected, T>=3':>16}")
    for n in sorted(df.n.unique()):
        s = df[df.n == n]
        s3 = s[s.n_time >= 3]
        print(f"  {n:>4} {s['published_d_temporal'].mean():>12.4f}"
              f" {s['temporal_last'].mean():>12.4f}"
              f" {s3['temporal_last'].mean():>16.4f}")

    print()
    print("=" * 78)
    print("G. BOOTSTRAP BASIS: 6 WORLD TYPES vs 18 FIELD REALISATIONS")
    print("=" * 78)
    for c in ("published_d_temporal", "temporal_last", "shap_TEMP"):
        a = boot_ci_worldtype(df, c)
        b = boot_ci(df, c)
        print(f"  {c:24s} worlds  {fmt(a)}  width {a[2]-a[1]:.4f}")
        print(f"  {'':24s} fields  {fmt(b)}  width {b[2]-b[1]:.4f}")

    with open(f"{R}/decomposition_corrected.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWROTE {R}/ablation_decomposed.csv and "
          f"{R}/decomposition_corrected.json")
