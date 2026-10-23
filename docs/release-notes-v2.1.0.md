No simulation result changes. This release adds the Shapley decomposition along the sample size gradient and reorders the main text figures so that the data limited result leads. Every results file from 2.0.0 is unchanged; two tables are added.

## New

* `analysis/02c_decompose_by_n.R` computes the exact Shapley decomposition over the 16 coalitions of {COORD, EDF, LAG, TEMP} within each per timestep sample size n (binary response, detection 1, spans 1, 3 and 5), with 95% intervals from a cluster bootstrap over field realisations (B = 2000), as in `02_decompose.R`. It writes:
  * `results/expanded_T9_decomposition_by_n.csv`, pooled over designs (main text Table 3, Figure 1a, Appendix S9 Table S16);
  * `results/expanded_T10_decomposition_by_n_design.csv`, within each design (Appendix S9 Table S17).
* `make decompose` now also runs `02c_decompose_by_n.R`.

## Figure numbering

Main text:

| Content | 2.0.0 | 2.1.0 | Script |
|---|---|---|---|
| Shapley decomposition by design; NN0 control | Figure 2 | Figure 1 (b, c), with a new panel (a) against sample size | `make_figures_1_3.R` |
| Temporal block by span, sample size and spatial coverage | Figure 4 | Figure 2 | `make_figure_2.R` |
| Robustness: detection; occurrence and abundance shares | Figure 5 | Figure 3 | `make_figures_1_3.R` |
| Five method comparison | Figure 3 | Figure 4 | `make_figure_4.R` |
| Case study | Figures 6 to 8 | Figures 5 to 7 | `make_figures_case.R` |

The schematic that was Figure 1 in 2.0.0 is no longer a main text figure; Supplement Figure S1 is the schematic.

Supplement:

| Content | 2.0.0 | 2.1.0 | Script |
|---|---|---|---|
| stRF versus sdmTMB (Appendix S6) | Figure S8 | Figure S5 | `make_figure_S5.R` |
| Case study (Appendix S7) | Figures S5 to S7 | Figures S6 to S8 | `make_figures_case.R` |

Figures S1 to S4 and S9 keep their numbers.

## Figure changes

* Figure 4: panel (a) shows mean AUC pooled over designs as bars drawn from the chance level (0.5), with 2.5 to 97.5 percentile whiskers of field realisation means. The clustered and stratified design effects sit below it in separate panels, (b) and (c), each on its own scale.
* Figure 2: points and intervals that share an x position are offset horizontally, and the dataset level points in panel (c) are jittered with a fixed seed so that the figure is reproducible.

## Scripts renamed

* `make_figures_2_5.R` is now `make_figures_1_3.R`.
* The temporal figure script `make_figure_4.R` is now `make_figure_2.R`.
* The method comparison script `make_figure_3.R` is now `make_figure_4.R`. A call to `make_figure_4.R` therefore draws the method comparison, not the temporal figure as in 2.0.0.
* `make_figure_S8.R` is now `make_figure_S5.R`.

`make figures` runs the renamed scripts.

## Metadata and documentation

* Version 2.1.0 in `DESCRIPTION`, `CITATION.cff` (released 2026-10-22) and `.zenodo.json`.
* README, NEWS, Makefile, vignette and script headers follow the new numbering. The 2.0.0 entry in NEWS is kept as released.

## Reproducing

```bash
make decompose   # Shapley tables T1 to T10 and the coverage analysis
make figures     # all figures from results/
```

## Citation

The concept DOI [10.5281/zenodo.23103919](https://doi.org/10.5281/zenodo.23103919) resolves to the latest archived version.

**Full changelog:** https://github.com/EAfrifa-Yamoah/marineSim/compare/v2.0.0...v2.1.0
