#!/usr/bin/env python3
"""revise2.py — edits that must work around inline OMML equation fields."""
from revise import apply, MS, SUPP

MS2 = [
("3.4 bootstrap basis",
 "estimated via cluster bootstrap over ground-truth worlds with ",
 "estimated by cluster bootstrap over independent field realisations with "),
("3.4 bootstrap rationale",
 " replicates. This treats the set of simulated worlds as a sample from a hypothetical population of marine ecosystems, producing confidence intervals that reflect the across-world variability in the decomposition.",
 " replicates. Resampling the six designed world types would treat a fixed factorial as a sample from a population. We therefore generated three independent field realisations of each world configuration and resampled those eighteen realisations instead. Intervals from the two bases are close in width, so this is a change of justification more than of magnitude, but only the realisation-based version supports a population-level reading."),

("5.5 scaling lead",
 "Its median fit time was approximately twice that of stRF at the sample sizes studied; because its cost scales cubically in the number of observations, the gap widened to roughly an order of magnitude at ",
 "Relative fit times depended strongly on sample size, and no single ratio describes them. In a paired timing experiment the comparator was faster than stRF below about 200 observations (0.04 s against 0.15 s at 30 observations) and slower above it, reaching 2.74 s against 0.29 s, a gap of roughly an order of magnitude, at "),
("5.5 scaling tail",
 " (See S6 of Supplementary material).",
 ". Its empirical cost grew as approximately the 1.5 power of sample size across the tested range rather than cubically, because fixed overhead still dominates at these sizes; the asymptotic scaling of a dense Gaussian process was not exercised here (See S6 of Supplementary material)."),

("6.1 uplift lead",
 "The average gain in predictive performance (",
 "The Shapley decomposition in Section 5.1 attributes 69.8% of the total gain to the three spatial feature blocks and 13.4% to algorithmic flexibility, a ratio of about five to one. The nested increments reported in earlier drafts ("),
("6.1 uplift tail",
 " AUC), indicating that representing spatial dependence is often more important than the choice of algorithm itself.",
 " AUC) are order dependent and are retained only for comparability. Representing spatial dependence matters considerably more than the choice of algorithm."),

("github url",
 "https://github.com/github.com/eafrifayamoah/marineSim",
 "https://github.com/eafrifayamoah/marineSim"),
("zenodo placeholder",
 "The full simulation results (~10k rows) are archived at Github and Zenado under CC-BY-4.0 with DOI 10.xxxx/figshare.placeholder.",
 "The full simulation results are archived on Zenodo under CC-BY-4.0; the archive DOI will be inserted on acceptance."),

# ---- Table 2: correct the feature descriptions (Issue 1) ----------------
("Table 2 stRF row",
 "ranger + spatial lags + temporal lags + EDF",
 "ranger, covariates + spatial lags + temporal lags + EDF (no raw coordinates)"),
("Table 2 SpatialRF row",
 "ranger, x, y + covariates",
 "ranger, covariates + x, y (coordinates only)"),

# ---- Table 3: replace the decomposition with Shapley values -------------
("T3 row1 label", "ML uplift (RF − GLM)", "Algorithmic flexibility"),
("T3 row2 label", "Spatial structure (sp. RF − RF)", "Spatial blocks (coordinates + distance fields + lags)"),
("T3 row3 label", "Temporal features (stRF − sp. RF)", "Temporal features"),
("T3 row4 label", "Total (stRF − GLM)", "Total gain over the aspatial GLM"),
("T3 ml est", "+0.015", "+0.025"),
("T3 ml ci", "[+0.004, +0.028]", "[+0.013, +0.038]"),
("T3 sp est", "+0.110", "+0.132"),
("T3 sp ci", "[+0.086, +0.135]", "[+0.118, +0.147]"),
("T3 tp est", "+0.039", "+0.032"),
("T3 tp ci", "[+0.028, +0.051]", "[+0.028, +0.035]"),
("T3 tot est", "+0.165", "+0.189"),
("T3 tot ci", "[+0.131, +0.198]", "[+0.167, +0.211]"),
("T3 share ml", "9%", "13.4%"),
("T3 share sp", "67%", "69.8%"),
("T3 share tp", "24%", "16.8%"),
("T3 mech tp", "Information borrowed across timesteps",
 "Information borrowed across timesteps, isolated by adding the temporal block alone"),
]

SUPP2 = [
("S3 figure ref",
 "held at all nine settings (Figure S1).",
 "held at all nine settings (Figure S3)."),
("S7 band lead",
 "Critically, the prediction band did not widen across this gap (mean width ",
 "The shaded band is across-tree ensemble spread, not a calibrated predictive interval. Its width did not increase across this gap (mean width "),
("S7 band tail",
 " shoots in the well sampled bands): the random forest reports false confidence where it extrapolates, so its own spread cannot be used to detect unsupported predictions.",
 " shoots in the well sampled bands). Empirical coverage confirms the problem. Against a nominal 80%, the band covered 74% of held-out observations under random cross-validation, while a quantile regression forest reached 83% and split conformal 83%. Under leave-one-region-out none of the three reached nominal coverage, at 63%, 68% and 47% respectively. The forest's own spread therefore cannot be used to detect unsupported predictions, and neither can a properly calibrated interval once prediction moves beyond the sampled bands: the area of applicability, not the interval, is what identifies where predictions are unsupported."),
("S7 sites",
 "yielding 748 site-years across 74 unique monitoring sites.",
 "yielding 748 site-years across 61 named monitoring sites, which carry 74 distinct coordinate pairs because some sites were re-positioned between survey rounds."),
]

if __name__ == "__main__":
    apply(MS, MS2, "MANUSCRIPT pass 2")
    apply(SUPP, SUPP2, "SUPPLEMENT pass 2")
