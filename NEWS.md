# marineSim 3.0.1 (2026-10-09)

* Figure 4 now resolves the temporal block by span, sample size and spatial
  coverage (tables T6 to T8, `02b_temporal_coverage.R`); the former pooled
  figure is Supplement Figure S9.
* Documented that the clustered design saturates its three cluster footprints at
  large n and is topped up with randomly placed cells (14% of sites at n = 100,
  25% at n = 200, 45% at n = 500); `results/clustered_random_fill.csv`.

# marineSim 3.0.0 (2026-10-01)

Release accompanying the submitted article. R is the sole computational engine.
The version follows the existing tags v2.0.0 and v2.0.1 of the July revision (whose
DESCRIPTION still read 0.3.0); the major bump marks the replacement of the engine.

* The package now contains the harmonised simulation engine of Section 3.1 to 3.4
  (ground truth generator with occurrence and Poisson abundance, three sampling
  designs, the four feature blocks, forest and GLM scoring, exact Shapley
  decomposition with cluster bootstrap, area of applicability). Function bodies
  are those of the scripts that produced the reported results; a rerun through
  the installed package reproduces the archived rows exactly.
* `analysis/` holds the numbered scripts that generate every number in the
  article and Supplement, including the expanded grid (three designs, six
  realisations, T = 10, detection 0.7, Poisson abundance), the GAM and sdmTMB
  comparators, and the R implementation of the seagrass case study.
* `results/` and `figures/` hold the reported outputs; `figures/scripts/` draws
  the figures from `results/` with ggplot2, so the whole repository is R.
* Removed: the earlier Python pipeline (`python/`), the `{targets}` pipeline and
  the previous `R/` sources, all of which implemented an earlier generator and
  are superseded. They remain in the git history before this release.
