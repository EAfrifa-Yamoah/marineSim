#!/usr/bin/env python3
"""
revise.py — apply the collaborator-driven corrections to document.xml.

Each edit is (label, old, new). The script reports match counts so any edit
that fails to land is visible rather than silently skipped.
"""
import re, sys, html

MS = "/home/claude/work/ms/word/document.xml"
SUPP = "/home/claude/work/supp/word/document.xml"


def apply(path, edits, label):
    s = open(path, encoding="utf-8").read()
    ok, miss = [], []
    for name, old, new in edits:
        c = s.count(old)
        if c == 1:
            s = s.replace(old, new)
            ok.append(name)
        else:
            miss.append((name, c))
    open(path, "w", encoding="utf-8").write(s)
    print(f"\n### {label}: {len(ok)} applied, {len(miss)} unmatched")
    for n, c in miss:
        print(f"    UNMATCHED ({c}x): {n}")
    return miss


# ============================================================ MANUSCRIPT
MS_EDITS = [

# ---- Abstract: design string and fit accounting (Issue 3) ----------------
("abstract design string",
 "Using a fully factorial simulation benchmark (6 ecological worlds × 5 sample sizes × 3 sampling designs × 2 detection probabilities × 3 temporal spans × 6 replicates × 5 methods; ~10,000 model fits), we quantified the independent contributions of spatial structure, temporal information, and algorithmic flexibility through hierarchical mixed-effects modelling and cluster-bootstrap confidence intervals.",
 "Using a factorial simulation benchmark (6 ecological worlds × 5 sample sizes × 3 sampling designs × 3 temporal spans × 6 replicates = 270 scenarios and 1,620 datasets, with imperfect detection, the generalised additive model and the geostatistical comparator each run as separate auxiliary grids; 10,068 model fits in total), we quantified the contributions of spatial structure, temporal information, and algorithmic flexibility. Because the feature blocks are correlated, contributions were estimated as Shapley values over all subsets of the four engineered blocks rather than as increments along a single nested ladder, which we show is not identified."),

# ---- Abstract: corrected decomposition (Issue 1) -------------------------
("abstract decomposition",
 "Spatial structure was the dominant driver of predictive performance, increasing AUC by +0.110 (95% CI [0.086, 0.135]), a gain approximately seven times larger than that achieved through increased machine-learning flexibility over a baseline GLM (+0.015 AUC). Temporal information provided a further average improvement of +0.039 AUC, but only when sufficient temporal depth was available, with gains emerging from three or more time points and being most pronounced at small to moderate sample sizes.",
 "Spatial structure was the dominant driver of predictive performance. The three spatial feature blocks together contributed +0.132 AUC, 69.8% of the total gain over an aspatial GLM, against 13.4% (+0.025 AUC) attributable to algorithmic flexibility, a ratio of about five to one. Temporal features contributed +0.032 AUC (95% CI [0.028, 0.035]), 16.8% of the gain, and only +0.007 AUC when added to an otherwise complete spatial model. The temporal contribution was conditional on temporal depth, appearing only from three or more time points, and was largest at the smallest sample sizes rather than at moderate ones."),

# ---- Abstract: case study claims (Issue 5) -------------------------------
("abstract case study",
 "Under random cross-validation, all methods performed well; however, when entire latitudinal bands were withheld, only spatial models retained substantial predictive skill (rank correlation = 0.35 for stRF versus 0.01 for aspatial models).",
 "Under random cross-validation all methods performed well (rank correlation 0.81 for stRF). When entire latitudinal bands were withheld, spatial models retained limited rank-transfer skill (mean within-park rank correlation 0.35 for stRF and 0.31 for a spatial random forest, against 0.21 for an aspatial forest), but absolute calibration failed in every case, with negative coefficients of determination throughout. Under rolling-origin evaluation stRF did not exceed a site-level climatology or a persistence baseline at any horizon."),

# ---- Section 2.5: additive decomposition is conceptual (minor fix) -------
("2.5 decomposition framing",
 "This motivates the skill decomposition, in which the AUC gain over the aspatial baseline is split into machine learning, spatial and temporal increments through a sequence of nested models,",
 "This motivates a skill decomposition of the AUC gain over the aspatial baseline. We stress that the expression above is a conceptual decomposition of the target, not a proof that the fitted forest is additive. Because the engineered blocks are correlated, an increment computed along a single nested ladder depends on the order in which blocks enter; Section 5.1 shows this dependence spans an order of magnitude. We therefore report Shapley values over all subsets of the four blocks as the primary decomposition, retaining the nested form below only to define the quantities,"),

# ---- Section 3.2: scenario accounting (Issue 3) --------------------------
("3.2 scenario totals",
 "At each combination of parameters, we drew 6 independent replicate samples from the ground-truth surface, yielding approximately 270 unique scenarios and ~10,000 model fits across all methods.",
 "Sample size, sampling design and the number of timesteps were fully crossed with the six worlds, giving exactly 270 core scenarios; each was replicated six times, yielding 1,620 datasets. Imperfect detection was examined as a separate reduced robustness grid rather than as a fully crossed factor, and the generalised additive model and the geostatistical comparator were run on their own auxiliary grids. Table 1 gives the complete accounting: 3,792 datasets and 10,068 model fits in total."),

# ---- Section 3.4: test sample size (Issue 3) -----------------------------
("3.4 test cells",
 "All methods were evaluated by prediction onto a large held-out test sample from the ground-truth surface (600 cells).",
 "All methods were evaluated by prediction onto a held-out sample of 1,500 cells drawn from the ground-truth surface, disjoint from the monitoring network and redrawn for every replicate. Test cells carry realised Bernoulli outcomes with detection fixed at 1.0, so AUC is computed against binary labels; the Brier score and the reliability diagrams are assessed against the same realised outcomes."),

# ---- Section 3.4: bootstrap basis (Issue 3) -----------------------------
("3.4 bootstrap",
 "Headline effect sizes (spatial-structure contribution, temporal-feature contribution, and pure-ML uplift) were estimated via cluster bootstrap over ground-truth worlds with $B\\  = \\ 1,000$ replicates. This treats the set of simulated worlds as a sample from a hypothetical population of marine ecosystems, producing confidence intervals that reflect the across-world variability in the decomposition.",
 "Headline effect sizes were estimated by cluster bootstrap with 2,000 replicates. Rather than resampling the six designed world types, which would treat a fixed factorial as a sample from a population, we generated three independent field realisations of each world configuration and resampled those eighteen realisations. Intervals from the two bases are close in width, so this is a change of justification more than of magnitude, but only the realisation-based version supports a population-level reading."),

# ---- Section 4.1 heading (minor fix) ------------------------------------
("4.1 heading",
 "Transfer learning strategies",
 "Validation design and spatial prediction"),

# ---- Section 5.5 heading and scaling claim (Issue 4) --------------------
("5.5 heading",
 "Comparison with sdmTMB",
 "Comparison with a Matérn Gaussian process comparator"),

("5.5 scaling claim",
 "Its median fit time was approximately twice that of stRF at the sample sizes studied; because its cost scales cubically in the number of observations, the gap widened to roughly an order of magnitude at $n\\  = \\ 500$ (See S6 of Supplementary material).",
 "Relative fit times depended strongly on sample size and no single ratio describes them. In a paired timing experiment the comparator was faster than stRF below about 200 observations (0.04 s against 0.15 s at n = 30) and slower above it, reaching 2.74 s against 0.29 s at n = 500, a gap of roughly an order of magnitude. Its empirical cost grew as approximately the 1.5 power of sample size over the tested range rather than cubically, because fixed overhead still dominates at these sizes; the asymptotic scaling of a dense Gaussian process was not exercised here (See S6 of Supplementary material)."),

# ---- Section 2.7: unmeasured sdmTMB timing (Issue 4) --------------------
("2.7 sdmTMB timing",
 "For comparison, fitting sdmTMB on the same data with a spatial random field takes ~30 s to several minutes depending on mesh complexity.",
 "We did not fit sdmTMB and therefore make no timing claim for it."),

# ---- Section 6.1: misstated ML uplift -----------------------------------
("6.1 misstated uplift",
 "The average gain in predictive performance ($\\sim 0.10$ AUC) substantially exceeded the improvement obtained by selecting a more sophisticated modelling algorithm ($\\sim 0.04$ AUC), indicating that representing spatial dependence is often more important than the choice of algorithm itself.",
 "The spatial feature blocks together contributed about 0.13 AUC, roughly five times the 0.025 AUC attributable to algorithmic flexibility, indicating that representing spatial dependence matters more than the choice of algorithm."),

# ---- Section 6.1: threshold recommendation (Issue 6) -------------------
("6.1 thresholds",
 "Our results suggest that temporal predictors become valuable only when datasets contain at least three timesteps and approximately 100 observations per timestep. Below these thresholds, temporal features are more likely to contribute noise than signal. Although conservative, this guideline provides a practical starting point for marine monitoring programmes and ecological forecasting applications, where temporal replication is frequently limited.",
 "Within the tested grid, temporal features contributed nothing at a single timestep and a small positive increment from three timesteps onwards. Their value was greatest at the smallest sample sizes and declined monotonically as sampling density increased, from +0.017 AUC at 30 observations per timestep to +0.001 AUC at 500. Temporal replication therefore substitutes for spatial replication rather than requiring it, and the practical condition is temporal depth, not a minimum sample size. We caution that these are patterns within the tested grid rather than universal thresholds."),

# ---- Conclusions (Issues 1, 5, 6) ---------------------------------------
("conclusions",
 "Through a fully factorial simulation benchmark of ~10,000 model fits, we have shown that in marine SDM under data limitation: (i) spatial structure contributes roughly seven times more to predictive performance than algorithm choice (about 67% versus 9% of the total gain); (ii) temporal features within a tree-based framework deliver conditional value, helping reliably only above a joint threshold of n ≥ 100 per timestep and ≥ 3 timesteps; (iii) sampling design effects are large and comparable in magnitude to method-family differences; and (iv) stRF is competitive on accuracy with a correctly specified geostatistical model while offering operational advantages (no mesh, linear scaling and reliable convergence), whereas a geostatistical model such as sdmTMB remains preferable when formal parametric uncertainty is required.",
 "Through a simulation benchmark of 10,068 model fits, we have shown that in marine SDM under data limitation: (i) spatial structure contributes about five times more to predictive performance than algorithm choice (69.8% versus 13.4% of the total gain), with the distance fields the single most valuable block; (ii) temporal features deliver a small conditional gain that requires at least three timesteps and is largest at the smallest sample sizes, and that falls to +0.007 AUC once a complete spatial representation is already present; (iii) sampling design effects are large and comparable in magnitude to method-family differences; and (iv) stRF is competitive on accuracy with a Matérn Gaussian process comparator while offering operational advantages, whereas a likelihood-based geostatistical model remains preferable when formal parametric uncertainty is required."),

("conclusions case study",
 "The seagrass case study demonstrates the framework in action, showing that the spatiotemporal features recover a latitudinal gradient that aspatial models miss, while area of applicability diagnostics honestly flag the unmonitored latitudinal gaps as unsupported extrapolation.",
 "The seagrass case study delimits what the framework can and cannot do on real data. Spatial features recover a latitudinal ranking that aspatial models miss, but absolute calibration fails outside the sampled bands, no interval construction we tested achieves nominal coverage there, and the method does not beat a site-level climatology at forecasting. Area of applicability diagnostics flag two thirds of the held-out predictions as unsupported extrapolation, which is the operative result: the diagnostic, not the model, is what makes the predictions safe to use."),

# ---- Data statement placeholders (minor fix) ----------------------------
("data statement",
 "All code is released as the R package marineSim with a reproducible  pipeline (Landau 2021) at [https://github.com/github.com/eafrifayamoah/marineSim]. The full simulation results (~10k rows) are archived at Github and Zenado under CC-BY-4.0 with DOI 10.xxxx/figshare.placeholder.",
 "All code is released as the R package marineSim with a reproducible pipeline (Landau 2021) at https://github.com/eafrifayamoah/marineSim. The full simulation results are archived on Zenodo under CC-BY-4.0."),

# ---- Citation fixes -----------------------------------------------------
("Fotheringham year",
 "(Bell et al. 2014; Fotheringham et al. 2009)",
 "(Bell et al. 2014; Fotheringham et al. 2002)"),

("Anderson 2022 a",
 "sdmTMB (Anderson et al. 2022), GAM tensor smooths",
 "sdmTMB (Anderson et al. 2025), GAM tensor smooths"),

("Anderson 2022 b",
 "and sdmTMB (Anderson et al. 2022), stRF offers",
 "and sdmTMB (Anderson et al. 2025), stRF offers"),

("Afrifa-Yamoah lettering",
 "(Afrifa-Yamoah et al., 2026; Costanza et al. 2014; Halpern et al. 2015)",
 "(Afrifa-Yamoah et al. 2026a; Costanza et al. 2014; Halpern et al. 2015)"),
]


# ============================================================ SUPPLEMENT
SUPP_EDITS = [
("S opening list",
 "This document provides supporting detail for the simulation design and ground-truth generator (S1), a hyperparameter sensitivity analysis (S2), full pairwise method contrasts (S3), probabilistic calibration diagnostics (S4), an extended comparison with the geostatistical model (S5), the case study data sources and harmonisation protocol (S6), and a reference summary of the marineSim R package (S7).",
 "This document provides the study design and stRF framework overview (S1), the simulation design and ground-truth generator (S2), a hyperparameter sensitivity analysis (S3), full pairwise method contrasts (S4), probabilistic calibration diagnostics (S5), an extended comparison with the geostatistical comparator (S6), the case study data sources and harmonisation protocol (S7), and a reference summary of the marineSim R package (S8)."),

("S2 test cells",
 "Model skill is evaluated on a fixed set of 1,500 cells held out from every fitted model, so training and evaluation never share a location.",
 "Model skill is evaluated on 1,500 cells held out from every fitted model and redrawn for each replicate, so training and evaluation never share a location."),

("S3 figure ref",
 "the qualitative ranking stRF > Spatial RF > standard RF held at all nine settings (Figure S1).",
 "the qualitative ranking stRF > Spatial RF > standard RF held at all nine settings (Figure S3)."),

("S5 figure ref",
 "(Table S5, Figure S2)",
 "(Table S5, Figure S4)"),

("S7 site count",
 "yielding 748 site-years across 74 unique monitoring sites.",
 "yielding 748 site-years across 61 named monitoring sites. The 61 sites carry 74 distinct coordinate pairs because some sites were re-positioned between survey rounds; spatial features use the coordinate recorded for each site-year, and the count of 74 in earlier versions referred to coordinate pairs rather than sites."),

("S7 band caveat",
 "Critically, the prediction band did not widen across this gap (mean width $13.8$ shoots inside versus $22.8$ shoots in the well sampled bands): the random forest reports false confidence where it extrapolates, so its own spread cannot be used to detect unsupported predictions.",
 "The across-tree spread is not a calibrated predictive interval and should be read as ensemble spread. Under random cross-validation its empirical coverage was 0.74 against a nominal 0.80, while a quantile regression forest achieved 0.83 and split conformal 0.83. Under leave-one-region-out none of the three reached nominal coverage: 0.63 for ensemble spread, 0.68 for the quantile regression forest and 0.47 for split conformal. The forest's own spread therefore cannot be used to detect unsupported predictions, and neither can a correctly calibrated interval once prediction moves outside the sampled bands."),
]

if __name__ == "__main__":
    m1 = apply(MS, MS_EDITS, "MANUSCRIPT")
    m2 = apply(SUPP, SUPP_EDITS, "SUPPLEMENT")
    sys.exit(0)
