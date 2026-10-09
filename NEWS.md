# marineSim 2.0.0 (2026-10-09)

Results change in this release: every clustered and stratified design number in
the article is replaced. Version 1.0.0 (Zenodo 10.5281/zenodo.23103920) remains
the record of the earlier generator.

* `sample_sites()` now implements the clustered design as preferential sampling
  (cells drawn without replacement with probability proportional to a three
  component Gaussian intensity), so exactly n distinct cells are returned and the
  design stays clustered at every n. The 1.0.0 generator rounded continuous
  draws to cells and topped up collisions with uniformly random cells, which made
  14% (n = 100), 25% (n = 200) and 45% (n = 500) of "clustered" sites random; it
  is retained as `"clustered_legacy"`. The stratified design now draws one cell
  per chosen block (previously about 3% of sites collided at n = 500 and were
  topped up at random; retained as `"stratified_legacy"`).
* The clustered and stratified arms of the benchmark (2,520 datasets) and their
  GAM comparator fits were rerun; random design rows are unchanged and reproduce
  the 1.0.0 archive exactly. `results/` and all affected figures and tables are
  updated. Six clustered draws with a degenerate training response were
  discarded (3,774 datasets, 79,254 fits).
* Figure 4 now resolves the temporal block by span, sample size and spatial
  coverage (tables T6 to T8, `02b_temporal_coverage.R`); the former pooled figure
  is Supplement Figure S9.
* `01_ablation_expanded.R` and `03_gam_core.R` accept `ABL_DESIGNS` to rerun a
  subset of designs; `05_si_figures_data.R`, `06_s6_sdmtmb.R` and `07_s6_pooled.R`
  attach ranger explicitly. The vignette is updated for the corrected designs.

# marineSim 1.0.0 (2026-10-02)

Release accompanying the submitted article. R is the sole computational engine.
The earlier tags v2.0.0 and v2.0.1 belong to the July revision (whose DESCRIPTION
still read 0.3.0); this release of the replaced engine is numbered 1.0.0.

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
