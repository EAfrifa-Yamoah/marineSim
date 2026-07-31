# stRF revision — completed changes

`Working_manuscript_REVISED.docx` and `Supplementary_Material_REVISED.docx` are fully updated: new results, corrected tables, re-rendered figures, new sections, renumbered headings and added references. Both validate against the originals (manuscript +5 paragraphs, supplement +2 — the new content). All equations, existing figures and formatting are preserved.

## Figures re-rendered and swapped into the documents

| Figure | Change |
|---|---|
| **Figure 2** | Replaced. (a) Shapley decomposition over all sixteen subsets of the four feature blocks; (b) placebo test by temporal depth, showing the nested contrast attributing +0.047 AUC where the temporal features are provably uninformative against +0.001 for the isolated increment. Frame resized to 6.5 × 3.55 in. Caption rewritten. |
| **Figure 6** | Replaced. Panel (a) now shows three validation regimes of increasing difficulty (random ten-fold, site-grouped, leave one region out) rather than a two-way interpolation/extrapolation split. Panels (b) and (c) rebuilt from the permutation and drop-column results. Frame resized. Caption rewritten. |
| **Figure S6** | Replaced with the refit rolling-origin results. stRF now shown below both persistence and climatology at every horizon. Caption corrected (the old horizon-0 panel no longer exists). |
| **Figure S7** | New. Empirical coverage and width of 80% intervals for ensemble spread, quantile regression forest and split conformal, under interpolation and under leave one region out. |

All figures use the Okabe-Ito colour-blind-safe palette, carry no embedded titles, and use a single shared legend per figure.

## Manuscript changes

**Abstract** — design string corrected (detection is not in the crossed product); 270 scenarios, 1,620 datasets, 10,068 fits; Shapley decomposition replaces the nested increments; case study claims moderated; rolling-origin null result added.

**Section 2.5** — additive expression relabelled a conceptual decomposition; order dependence stated; Shapley named as primary.

**Section 2.7** — unmeasured sdmTMB timing claim deleted.

**Section 3.2** — exact scenario and fit accounting.

**Section 3.4** — 1,500 test cells (not 600), redrawn per replicate, realised binary labels stated; bootstrap re-based on eighteen independent field realisations.

**Section 4.1** — retitled "Validation design and spatial prediction". Webster et al. (2024) now cited at first mention of the data, with the 61 sites / 74 coordinate pairs explanation.

**Section 5.1** — new paragraph on order dependence (+0.091 first, +0.052 nested, +0.007 last) and the placebo test.

**Section 5.5** — retitled "Comparison with a Matérn Gaussian process comparator"; cubic scaling replaced with the measured 1.5 power and the actual paired timings.

**Section 5.6** — interpolation numbers updated (0.75 / 0.80 / 0.81) with the site-grouped intermediate added; LORO reported both within-park (0.31, 0.35 against 0.24 and 0.21) and pooled, with the pooling artefact explained; R² corrected from "close to zero" to negative throughout.

**Section 5.7** — renumbered from the duplicate 5.6.

**Sections 5.8 and 5.9** — new. Forecast skill, and predictive uncertainty.

**Section 6.1** — misstated 0.04 AUC uplift corrected; threshold recommendation inverted.

**Conclusions** — five to one ratio, distance fields named as the largest block, case study limits stated plainly.

**Tables** — Table 1 detection and methods rows plus caption corrected; Table 2 feature descriptions corrected (stRF carries no raw coordinates) and the comparator relabelled a dense Matérn GP; Table 3 rebuilt as Shapley values with shares 13.4 / 69.8 / 16.8%.

**References** — Fotheringham 2002, Anderson 2025 (both), Afrifa-Yamoah 2026a; Valavi et al. (2019), MacKenzie et al. (2002) and tinyVAST added; data statement placeholders removed.

## Supplement changes

Opening section list corrected to S1–S8; S2 test cells clarified; S3 and S5 figure cross references fixed; S7 site count reconciled; S7 prediction band rewritten around measured coverage; forecast text updated with the refit numbers; Figure S7 added.

## Headline numbers

**Simulation** — Shapley shares: distance fields 32.9%, spatial lags 19.0%, coordinates 17.9%, temporal 16.8%, algorithmic flexibility 13.4%. Spatial total 69.8%, a ratio of 5.2 to 1 over algorithm choice. Temporal +0.032 Shapley, +0.007 added last, +0.001 at one timestep. Temporal gain declines with sample size: +0.017 at n = 30 to +0.001 at n = 500.

**Case study** — 61 named sites, 74 coordinate pairs, 748 site years. LORO within-park ρ: stRF 0.353, Spatial RF 0.313, GLM 0.244, RF 0.206; R² negative for all. Temporal feature unavailable for 100% of held-out rows under LORO; spatial lags fall back for 94%, 33% and 67% of rows in Jurien Bay, Marmion and Ngari Capes. Rolling origin: stRF below both baselines at every horizon. Interval coverage: 0.743 / 0.832 / 0.826 under interpolation, 0.628 / 0.680 / 0.472 under extrapolation, with 65.3% of held-out rows outside the area of applicability.

## Still outstanding

The R package under `marineSim/R/` has not been updated — there is no R runtime in this environment, so `features.R` and `fit.R` still need the block coalition interface added to match what the paper now reports. The case study `stRF` includes coordinates while the simulation `stRF` does not; the manuscript now states the simulation specification correctly in Table 2, but you may want to harmonise the two implementations.
