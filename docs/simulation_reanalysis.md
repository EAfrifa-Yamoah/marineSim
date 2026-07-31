# stRF revision: re-analysis results

**Status.** Run against your own `python/sim_engine.py` — the ground truth generator, feature functions and RF backbone are yours, unmodified. What is new is the experimental design wrapped around them: 32 feature coalitions instead of 5 fixed models, 3 independent field realisations per world configuration instead of 1, and three targeted verification experiments. Nothing here is a reimplementation.

**Scope.** 18 field realisations (6 configurations × 3 realisations) × 5 sample sizes × 3 temporal spans, random design, detection 1.0 = 270 datasets, 5,400 forest fits, plus 45 GP timing fits and 72 coverage fits. Absolute values differ from the published grid because this run uses the random design only and a different realisation set; the published contrast reproduces at +0.052 here against +0.039 in `derived.json`. Signs, orderings and the qualitative conclusions are what transfer.

---

## Issue 1 — the temporal contribution

**What the published contrast changes.** From `sim_engine.fit_eval`:

- Spatial RF (L328–334) = covariates + `COORD`
- stRF (L336–375) = covariates + `LAG` + `TEMP` + `EDF`, **no coordinates**

`aggregate.py` L35 sets `d_temporal = stRF − SpatialRF`, which adds three blocks and removes one.

**The placebo test settles it.** At T = 1 the temporal features are set to the training mean by construction (`sim_engine.py` L360–362), so they carry no information and any honest estimator must return zero.

| Timesteps | Published contrast | Corrected increment (temporal block only) |
|---:|---:|---:|
| 1 | **+0.0473** | +0.0011 |
| 3 | +0.0540 | +0.0103 |
| 5 | +0.0531 | +0.0104 |

The published contrast attributes +0.047 AUC to temporal features in the one condition where they are provably uninformative. The corrected increment passes the placebo and shows the real temporal effect switching on at T ≥ 3, which is the mechanism Section 2.2 claims.

**Sequential attribution is not identified.** The increment assigned to the temporal block depends entirely on where it enters:

| Attribution | Temporal increment |
|---|---:|
| Added first, to covariates only | +0.0912 [+0.0811, +0.1010] |
| Published contrast | +0.0515 [+0.0441, +0.0586] |
| **Shapley (order free)** | **+0.0317 [+0.0281, +0.0352]** |
| Added last, to a complete model | +0.0072 [+0.0051, +0.0094] |

A thirteenfold spread across orderings. Any nested ladder — including the one the collaborator proposed — inherits this. The Shapley value is the defensible headline; the added-last figure is the one a practitioner needs, since it answers "I already have a spatial model, is it worth building temporal features?"

**Corrected Shapley decomposition** (total gain over the aspatial GLM = +0.189 AUC):

| Block | Shapley ΔAUC | 95% CI | Share |
|---|---:|---|---:|
| Distance fields | +0.0621 | [+0.0564, +0.0677] | 32.9% |
| Spatial lags | +0.0358 | [+0.0314, +0.0404] | 19.0% |
| Coordinates | +0.0338 | [+0.0296, +0.0384] | 17.9% |
| Temporal features | +0.0317 | [+0.0281, +0.0352] | 16.8% |
| Algorithmic flexibility | +0.0253 | [+0.0129, +0.0384] | 13.4% |
| **Spatial total** | **+0.1317** | | **69.8%** |

The paper's central claim survives and is if anything cleaner: spatial representation is 69.8% of the gain against 13.4% for algorithm choice, a ratio of 5.2 to 1. Temporal drops from 24% to 16.8%. The most valuable single block is the distance fields, which the manuscript currently never isolates.

Replacement Figure 2 is in `Figure2_corrected_decomposition.png`.

## Issue 2 — leakage: PASS

88 predictor comparisons across 4 worlds and 2 sample sizes. Held-out responses were permuted and inverted, then every feature block rebuilt: **zero changes**. A non-vacuity control confirms the test can detect dependence (perturbing a *training* response does move the test lag features, as it must).

Verified by reading: `spatial_lags(..., is_train=False)` uses training responses only; `temporal_lag` uses training site prior responses only; `y_test` reaches nothing but `roc_auc_score`. Self exclusion is `np.fill_diagonal(mask, False)` (L219).

State this as a positive result with the assertion test cited, not as a defensive paragraph. One genuine omission to fix: the forest trains on **n rows at the evaluation timestep only** (`y_tr = obs_y[t_eval]`); earlier timesteps feed the temporal features and nothing else. The manuscript never says this.

## Issue 3 — accounting

Exact, from `benchmark.py` (`DETECTION = [1.0]`, 4 methods) and `aux_runs.py`:

| Run | Datasets | Fits |
|---|---:|---:|
| Core grid, 4 methods (6 × 5 × 3 × 3 × 6) | 1,620 | 6,480 |
| GAM auxiliary, same cells | 1,620 | 1,620 |
| GP comparator (random, T = 5, 4 reps) | 120 | 240 |
| Detection 0.7 robustness (3 n × 2 designs × 2 T × 6 reps × 6 worlds) | 432 | 1,728 |
| **Total** | **3,792** | **10,068** |

"~10,000 fits" is correct. Core scenarios = 6 × 5 × 3 × 3 = **270 exactly**. Fix the abstract's design string (detection is not in the crossed product) and Table 1's "fully crossed". `derived.json`'s `n_fits_total: 3791` counts datasets, not fits (`aggregate.py` L173) — relabel.

`N_TEST = 1500` (`sim_engine.py` L38): the Supplement is right, the main text's 600 is wrong.

`win_pct` is computed on the 4 method main grid, which excludes GAM — Table 2's 0.0% for GAM comes from a different grid than the other entries.

**Bootstrap basis.** The collaborator is right in principle but the practical impact is small. Resampling 18 independent field realisations instead of 6 designed world types barely moves the intervals (width 0.0145 against 0.0154 for the published contrast). Report the realisation based version because it is defensible, not because it changes anything.

## Issue 4 — the geostatistical comparator

The comparator is `sklearn.gaussian_process.GaussianProcessClassifier` with `ConstantKernel × Matern(ν=1.5) + WhiteKernel` on [coords, covariates] (L377–397). Nothing resembling sdmTMB is fitted. Rename or fit the real thing.

Timing experiment, 45 paired fits over 3 worlds:

| n | GP fit (s) | stRF fit (s) |
|---:|---:|---:|
| 30 | 0.039 | 0.150 |
| 50 | 0.051 | 0.154 |
| 100 | 0.081 | 0.158 |
| 200 | 0.267 | 0.183 |
| 500 | 2.743 | 0.289 |

- GP log-log slope on n = **1.49**, not 3. The cubic claim is not supported over the tested range.
- stRF slope = 0.22, dominated by fixed tree building overhead — so "linear scaling" is also not what is observed here.
- GP/stRF ratio at n = 500 = **9.5×**, so "roughly an order of magnitude at n = 500" *is* supported.
- The GP is **faster** than stRF below about n = 200. The claim that its median fit time was approximately twice that of stRF is not supportable in either direction; report the curve.

Delete the unmeasured "~30 s to several minutes" for sdmTMB in Section 2.7.

## Issue 6 — the prediction band and the thresholds

**Band coverage**, measured against the known true occurrence probability (nominal 80%):

| Quantity | Value |
|---|---:|
| Overall coverage | **0.725** |
| Coverage inside AOA | 0.721 |
| Coverage outside AOA | 0.938 |
| Mean width inside AOA | 0.526 |
| Mean width outside AOA | 0.895 |

The band is miscalibrated in both regions and in opposite directions: it under-covers where the data are and over-covers where they are not. Width does respond to support in the simulation, which differs from the case study finding — worth stating, since the two are measuring different targets (true probability against observed density). Either relabel as ensemble spread or build and validate a real interval.

**The n ≥ 100 threshold is refuted by the corrected increment.** The temporal gain declines monotonically with sample size:

| n | Published | Corrected | Corrected, T ≥ 3 |
|---:|---:|---:|---:|
| 30 | +0.0625 | +0.0119 | +0.0168 |
| 50 | +0.0645 | +0.0136 | +0.0188 |
| 100 | +0.0576 | +0.0060 | +0.0089 |
| 200 | +0.0441 | +0.0042 | +0.0060 |
| 500 | +0.0285 | +0.0006 | +0.0010 |

Temporal features are most valuable at the **smallest** samples, which inverts the Section 6.1 and Conclusions recommendation. Keep the ≥ 3 timesteps condition; drop the n ≥ 100 condition and replace it with the opposite guidance, stated as a pattern within the tested grid.

## One hypothesis I raised and then refuted

I predicted that `TEMP` was partly acting as a fine scale spatial nearest neighbour feature, because `temporal_lag` looks up the previous response at the nearest *training site*. Adding a same timestep 1-NN control block (`NN0`) tests this. It does not hold:

| Quantity | ΔAUC |
|---|---:|
| TEMP marginal, no NN0 present | +0.0072 |
| TEMP marginal, NN0 present | +0.0074 |
| NN0 marginal alone | −0.0000 |

A same timestep nearest neighbour lag absorbs none of the temporal contribution and buys nothing itself — the 25 km ring already captures it. The residual temporal effect, though small, is **genuinely temporal**. This strengthens the mechanism in Section 2.2 and is worth reporting as a control in the Supplement.

---

## Files

- `Figure2_corrected_decomposition.png` — replacement Figure 2
- `ablation.py`, `decompose.py`, `verify.py`, `fig_corrected_decomposition.py` — drop into `python/`
- `ablation_decomposed.csv` — per dataset coalition AUCs and all derived increments
- `decomposition_corrected.json`, `verification.json`, `gp_scaling.csv`, `band_coverage.csv`

## Not addressed

Issue 5 in full, and the case study half of Issue 6: the raw DBCA seagrass data are not in the archive. The derived files that *are* present (`case_loro.csv`, `case_dropcol.csv`, `case_grid_unc.csv`, `case_interp_extrap.csv`, `temporal_skill.csv`) are enough to re-check the transfer skill and forecast claims without them. The R package under `marineSim/R/` has not been brought into line with these changes — no R runtime here — so `features.R` and `fit.R` will need the block coalition interface added by hand.
