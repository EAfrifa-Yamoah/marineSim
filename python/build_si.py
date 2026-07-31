import re, shutil, os
import json
from PIL import Image

U = "si_build"; DOCPATH = f"{U}/word/document.xml"
doc = open(DOCPATH).read()
preamble = doc[:doc.index("<w:body>") + len("<w:body>")]
sectpr = doc[doc.rfind("<w:sectPr"):]           # sectPr + </w:body></w:document>
img_tpl = open("/tmp/img_template.xml").read()

def esc(s): return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

def run(t, b=False, sz=None):
    rpr = ""
    if b or sz:
        rpr = "<w:rPr>" + ("<w:b/><w:bCs/>" if b else "") + \
              (f'<w:sz w:val="{sz}"/><w:szCs w:val="{sz}"/>' if sz else "") + "</w:rPr>"
    return f'<w:r>{rpr}<w:t xml:space="preserve">{esc(t)}</w:t></w:r>'

def para(runs, jc="both", before=0, after=140):
    sp = f'<w:spacing w:after="{after}" w:line="276"' + (f' w:before="{before}"' if before else "") + "/>"
    return f'<w:p><w:pPr>{sp}<w:jc w:val="{jc}"/></w:pPr>{"".join(runs)}</w:p>'

def P(text): return para([run(text)])
def PB(parts): return para([run(t, b=b) for (t, b) in parts])
def H(text): return para([run(text, b=True, sz=26)], jc="left", before=240, after=120)
def TITLE(text): return para([run(text, b=True, sz=32)], jc="center", after=120)
def CAP(label, text): return para([run(label, b=True), run(text)])

def cell(t, w, header=False):
    if header:
        rpr = '<w:b/><w:bCs/><w:color w:val="FFFFFF"/><w:sz w:val="19"/><w:szCs w:val="19"/>'
        shd = '<w:shd w:fill="2C3E50" w:val="clear"/>'
    else:
        rpr = '<w:b w:val="false"/><w:bCs w:val="false"/><w:color w:val="000000"/><w:sz w:val="19"/><w:szCs w:val="19"/>'
        shd = ''
    bd = ('<w:tcBorders><w:top w:val="single" w:color="666666" w:sz="4"/>'
          '<w:left w:val="single" w:color="666666" w:sz="4"/>'
          '<w:bottom w:val="single" w:color="666666" w:sz="4"/>'
          '<w:right w:val="single" w:color="666666" w:sz="4"/></w:tcBorders>')
    mar = ('<w:tcMar><w:top w:type="dxa" w:w="80"/><w:left w:type="dxa" w:w="120"/>'
           '<w:bottom w:type="dxa" w:w="80"/><w:right w:type="dxa" w:w="120"/></w:tcMar>')
    return (f'<w:tc><w:tcPr><w:tcW w:type="dxa" w:w="{w}"/>{bd}{shd}{mar}</w:tcPr>'
            f'<w:p><w:pPr><w:spacing w:after="0" w:line="264"/><w:jc w:val="left"/></w:pPr>'
            f'<w:r><w:rPr>{rpr}</w:rPr><w:t xml:space="preserve">{esc(t)}</w:t></w:r></w:p></w:tc>')

def table(widths, header, rows):
    grid = "".join(f'<w:gridCol w:w="{w}"/>' for w in widths)
    tblpr = ('<w:tblPr><w:tblW w:type="dxa" w:w="9360"/><w:tblBorders>'
             '<w:top w:val="single" w:color="auto" w:sz="4"/><w:left w:val="single" w:color="auto" w:sz="4"/>'
             '<w:bottom w:val="single" w:color="auto" w:sz="4"/><w:right w:val="single" w:color="auto" w:sz="4"/>'
             '<w:insideH w:val="single" w:color="auto" w:sz="4"/><w:insideV w:val="single" w:color="auto" w:sz="4"/>'
             '</w:tblBorders></w:tblPr>')
    hdr = f'<w:tr><w:trPr><w:tblHeader/></w:trPr>{"".join(cell(t, w, True) for t, w in zip(header, widths))}</w:tr>'
    body = "".join(f'<w:tr><w:trPr><w:tblHeader w:val="false"/></w:trPr>{"".join(cell(t, w) for t, w in zip(r, widths))}</w:tr>' for r in rows)
    return f'<w:tbl>{tblpr}<w:tblGrid>{grid}</w:tblGrid>{hdr}{body}</w:tbl>'

IMAGES = []
def add_image(src, media, cx_fixed=None, dpi=None, cap=5029200):
    w, h = Image.open(src).size
    cx = cx_fixed if cx_fixed is not None else min(round(w / dpi * 914400), cap)
    cy = round(cx * h / w)
    rid = f"rId{200 + len(IMAGES)}"; docpr = 200 + len(IMAGES)
    IMAGES.append((rid, media, src))
    x = (img_tpl.replace('r:embed="rId17"', f'r:embed="{rid}"')
                .replace('wp:docPr id="8"', f'wp:docPr id="{docpr}"')
                .replace('name="fig8_design_overview.png"', f'name="{media}"')
                .replace('title="fig8_design_overview.png"', f'title="{media}"'))
    x = re.sub(r'descr="[^"]*"', f'descr="{media}"', x)
    x = re.sub(r'cx="\d+" cy="\d+"', f'cx="{cx}" cy="{cy}"', x)
    return x

EQ_OMML = json.load(open("eq_omml.json"))
def EQ(name):
    return (f'<w:p><w:pPr><w:spacing w:before="60" w:after="140" w:line="276"/></w:pPr>'
            f'{EQ_OMML[name]}</w:p>')

MINUS = "\u2212"; TIMES = "\u00d7"; GE = "\u2265"; NU = "\u03bd"; RHO = "\u03c1"
body = []
A = body.append

A(TITLE("Supplementary Material"))
A(para([run("Spatio-temporal random forests with adaptive local features for marine species "
            "distribution modelling under data limitation", b=True)], jc="center", after=60))
A(para([run("Ebenezer Afrifa-Yamoah, Francky Fouedjio, Ute Mueller, Kathryn McMahon, Glenn A. Hyndes",
            sz=20)], jc="center", after=200))
A(P("This document provides the mathematical formulation of the stRF model (S1), supporting detail for "
    "the simulation design and ground-truth generator (S2), a hyperparameter sensitivity analysis (S3), "
    "full pairwise method contrasts (S4), probabilistic calibration diagnostics (S5), an extended "
    "comparison with the geostatistical model (S6), the case study data sources and harmonisation "
    "protocol (S7), and a reference summary of the marineSim R package (S8). All analyses are "
    "reproducible from the code repository."))

# ---------------- S1 (mathematical formulation) ----------------
A(H("S1 | Mathematical formulation of the spatio-temporal random forest (stRF)"))
A(P("Let the study region be a subset of the plane and let Y(s, t) denote the response at location s and "
    "timestep t; in the simulation Y is a binary occurrence indicator and in the case study Y is the "
    "natural logarithm of one plus shoot density. Training data are observations y_i = Y(s_i, t) recorded "
    "at n network sites that are revisited over T timesteps. stRF is an ordinary random forest applied to "
    "an augmented feature map that couples the measured covariates with explicit encodings of spatial and "
    "temporal structure. The feature map at a query point (s, t) is"))
A(EQ("eq1_features"))
A(P("where z(s, t) collects the p environmental covariates, the terms L_1 to L_M are multi scale spatial "
    "lags, the terms tau_1 and tau_H are temporal features, and e_1 to e_K are Euclidean distance fields. "
    "Each block is defined in turn."))
A(PB([("Spatial lags. ", True), ("The lag at radius r_m is the mean observed response over the "
      "neighbourhood of s within that radius,", False)]))
A(EQ("eq2_lag"))
A(P("evaluated at radii 25, 75 and 150 km. These features generalise the nearest observation predictors of "
    "random forest spatial interpolation (Sekulic et al. 2020): short radii act as a local autocovariate "
    "that captures fine scale clustering and long radii encode the broad spatial trend. Self inclusion is "
    "excluded at training points to prevent leakage."))
A(PB([("Temporal features. ", True), ("Two features summarise a site's own history, the response at the "
      "previous timestep and the mean over all earlier timesteps at the same site,", False)]))
A(EQ("eq3_temporal"))
A(P("For a single timestep both default to the training mean, so their informativeness grows with the "
    "number of revisits; this is the mechanism behind the timestep gradient in the main results."))
A(PB([("Distance fields. ", True), ("The distance fields are Euclidean distances to K anchor points,", False)]))
A(EQ("eq4_edf"))
A(P("with the anchors placed at the centroids of a k-means partition of the training coordinates (K = 15). "
    "Distance to a fixed set of reference points is the device introduced by RFsp (Hengl et al. 2018) and "
    "by the bounding box distances of Behrens et al. (2018); it gives the trees a smooth, low dimensional "
    "coordinate embedding to split on, approximating a spatial trend surface without exposing raw "
    "coordinates that extrapolate poorly."))
A(PB([("Random forest predictor. ", True), ("The augmented features are supplied to a forest of B trees. "
      "The ensemble prediction and the piecewise constant form of each tree are", False)]))
A(EQ("eq5_forest"))
A(P("where the regions R are the leaves of tree b and their fitted values are leaf class proportions for "
    "occurrence or leaf means for log density. Trees are grown on bootstrap samples; at each node a random "
    "subset of the p + M + 2 + K features is considered and the split minimises the node loss (Gini "
    "impurity for classification, residual variance for regression),"))
A(EQ("eq6_split"))
A(P("Because spatial and temporal signal enter as ordinary features, standard forest machinery and default "
    "hyperparameters apply; S3 shows the ranking is insensitive to these settings."))
A(PB([("Additive interpretation and skill decomposition. ", True), ("Conceptually the augmented forest "
      "approximates an additive decomposition of the expected response into a covariate effect, a residual "
      "spatial component and a temporal component,", False)]))
A(EQ("eq7_additive"))
A(P("where g is the nonlinear covariate effect captured by the tree partitions, the spatial term is the "
    "residual autocorrelation absorbed by the spatial lags and distance fields, and the temporal term is "
    "the cross timestep signal absorbed by the temporal features. This motivates the skill decomposition "
    "in the main text, in which the AUC gain over the aspatial baseline is split into machine learning, "
    "spatial and temporal increments through a sequence of nested models,"))
A(EQ("eq8_decomp"))
A(P("with A denoting the area under the ROC curve of each model (GLM, standard RF, spatial RF and stRF); "
    "the three bracketed terms are the machine learning, spatial and temporal contributions respectively."))
A(PB([("Area of applicability. ", True), ("Predictions are accompanied by an area of applicability "
      "diagnostic (Meyer and Pebesma 2021) that flags query points outside the region of feature space "
      "spanned by the training data. Features are standardised and weighted by their permutation "
      "importance, and a dissimilarity index is the distance to the nearest training point relative to the "
      "mean pairwise training distance,", False)]))
A(EQ("eq9_di"))
A(P("with the importance weighted distance"))
A(EQ("eq10_diw"))
A(P("where the mean and standard deviation are those of feature q in the training data, I_q is its "
    "permutation importance, and the tilde denotes the standardised weighted feature vector. A query point "
    "is flagged as extrapolation when its dissimilarity index exceeds the 95th percentile of the training "
    "dissimilarities, the threshold used to delineate the unmonitored gaps in the case study."))

# ---------------- S2 ----------------
A(H("S2 | Simulation design and ground-truth generator"))
A(P(f"Each ground-truth world is generated on a 50 {TIMES} 50 lattice over a 400 {TIMES} 400 km domain "
    f"(8 km cells). Spatial fields are drawn as Matern Gaussian random fields by fast Fourier spectral "
    f"synthesis with smoothness {NU} = 3/2. Four measured covariates (depth, substrate, wave exposure and "
    f"sea surface temperature) are simulated as short to moderate range fields (35 to 60 km), deliberately "
    f"rougher than the broad unmeasured spatial process so that tree ensembles cannot use the covariate "
    f"vector as a smooth proxy for location. Sea surface temperature evolves over time as an AR(1) process "
    f"({RHO} = 0.6)."))
A(P("The occurrence probability at each cell and timestep is the logistic transform of a linear predictor "
    "that combines the covariates (with a mild depth optimum and a depth by wave interaction, the controlled "
    "source of pure machine learning uplift), an unmeasured spatial random effect at the world's "
    "autocorrelation range, and a spatially rough, temporally persistent site effect that only a site's own "
    "temporal history can resolve. A spatially clustered suitable-habitat mask sets approximately a quarter "
    "of the domain to structural (delta) zeros. Nonstationary worlds additionally vary the depth and "
    "temperature coefficients smoothly across space."))
A(P(f"The factorial crosses two stationarity settings {TIMES} three autocorrelation ranges (30, 80, 200 km) "
    f"= six worlds, five per timestep sample sizes (30, 50, 100, 200, 500), three sampling designs (random, "
    f"clustered, stratified) and three timestep counts (1, 3, 5), giving 270 core scenarios, each replicated "
    f"six times (1,620 datasets). Detection probability (0.7, 1.0) is examined as a separate robustness "
    f"analysis. The generator is provided as make_world() and matern_grf() in the marineSim package (S7) and "
    f"is illustrated in Figure 8 of the main text."))

# ---------------- S2 ----------------
A(H("S3 | Hyperparameter sensitivity"))
A(P("The main benchmark fits every method with default hyperparameters to test robustness as drop-in tools. "
    "To confirm that this choice does not bias the comparison, the three random forest methods were refit "
    f"across a 3 {TIMES} 3 grid of the two most influential settings, minimum node size (1, 3, 5) and the "
    "number of features considered per split (the square root of the feature count, 0.5 and 0.8 of it), on "
    "six worlds at n = 50, 100 and 200 with the same data supplied to every setting. Test AUC varied by at "
    f"most 0.015 within each method, the default setting was within 0.002 of the best setting, and the "
    "qualitative ranking stRF > Spatial RF > standard RF held at all nine settings (Figure S1). Tuning would "
    "therefore not change the conclusions while inflating computational cost by about two orders of magnitude."))
A(table([2600, 1700, 1700, 1700, 1660],
        ["Method", "Default AUC", "Best AUC", "Worst AUC", "Range"],
        [["Standard RF", "0.677", "0.677", "0.672", "0.005"],
         ["Spatial RF", "0.797", "0.797", "0.783", "0.014"],
         ["stRF", "0.844", "0.846", "0.832", "0.015"]]))
A(CAP("Table S1. ", "Sensitivity of the random forest family to hyperparameters. AUC is averaged over six "
      "worlds and three sample sizes; the default setting is minimum node size 1 with the square root feature "
      "rule. The range is the spread across all nine settings."))
A(add_image("figures/FigureS1_tuning.png", "FigureS1_tuning.png", cx_fixed=4114800))
A(CAP("Figure S1. ", "Test AUC across the nine hyperparameter settings for each random forest method. Each "
      "point is one setting (averaged over worlds and sample sizes), the diamond marks the default setting and "
      "the bar the mean. The bands do not overlap, so the ranking is invariant to tuning."))

# ---------------- S3 ----------------
A(H("S4 | Full pairwise method contrasts"))
A(P("Table S2 reports mean AUC by sample size for all five methods. The generalised additive model collapses "
    "at the smallest samples (AUC near 0.50 at n = 30) and becomes competitive only from n = 200, whereas "
    "stRF leads at every sample size. Table S3 gives the stRF advantage over each method with a 95% paired "
    "bootstrap interval; every interval excludes zero. The gain of stRF over Spatial RF, which isolates the "
    "temporal features, is 0.026 to 0.046 AUC across sizes. Table S4 shows the effect of sampling design at "
    "n = 100: clustered sampling costs roughly 0.06 to 0.09 AUC relative to stratified for every method, "
    "while stratified is marginally better than random."))
A(table([1560, 1560, 1560, 1560, 1560, 1560],
        ["n per timestep", "GLM", "GAM", "RF", "Spatial RF", "stRF"],
        [["30", "0.624", "0.502", "0.614", "0.700", "0.743"],
         ["50", "0.644", "0.510", "0.639", "0.736", "0.778"],
         ["100", "0.663", "0.573", "0.670", "0.782", "0.828"],
         ["200", "0.682", "0.803", "0.707", "0.834", "0.872"],
         ["500", "0.693", "0.850", "0.752", "0.882", "0.908"]]))
A(CAP("Table S2. ", "Mean test AUC by method and per timestep sample size, pooled over worlds, designs and "
      "timestep counts."))
A(table([1900, 1500, 1500, 1500, 1480, 1480],
        ["Contrast", "n = 30", "n = 50", "n = 100", "n = 200", "n = 500"],
        [[f"stRF {MINUS} GLM", "+0.119", "+0.134", "+0.165", "+0.191", "+0.215"],
         [f"stRF {MINUS} GAM", "+0.241", "+0.268", "+0.254", "+0.069", "+0.058"],
         [f"stRF {MINUS} RF", "+0.129", "+0.139", "+0.157", "+0.165", "+0.156"],
         [f"stRF {MINUS} Spatial RF", "+0.043", "+0.042", "+0.046", "+0.038", "+0.026"]]))
A(CAP("Table S3. ", "Mean AUC advantage of stRF over each method by sample size. All 95% paired bootstrap "
      "intervals (2,000 resamples) exclude zero; the largest half width is 0.014."))
A(table([2400, 2320, 2320, 2320],
        ["Sampling design", "Standard RF", "Spatial RF", "stRF"],
        [["Clustered", "0.627", "0.733", "0.772"],
         ["Random", "0.687", "0.799", "0.847"],
         ["Stratified", "0.697", "0.813", "0.863"]]))
A(CAP("Table S4. ", "Effect of sampling design on mean test AUC at n = 100 (GLM and GAM omitted for brevity; "
      "the pattern is the same)."))

# ---------------- S4 ----------------
A(H("S5 | Probabilistic calibration"))
A(P("Calibration was assessed on six worlds at n = 100 (three replicates, random design) by comparing "
    "predicted occurrence probabilities against the known ground truth. stRF was the best calibrated method "
    "by the Brier score, followed by Spatial RF; the aspatial models were substantially worse and the "
    "generalised additive model slightly over-confident at low probabilities (Table S5, Figure S2). Better "
    "calibration tracks the AUC ranking, indicating that the spatial and temporal features improve the "
    "probability estimates and not only the ranking of sites."))
A(table([3000, 3180, 3180],
        ["Method", "Brier score (mean)", "Brier score (SD)"],
        [["stRF", "0.153", "0.031"],
         ["Spatial RF", "0.171", "0.027"],
         ["Standard RF", "0.217", "0.014"],
         ["GLM", "0.217", "0.013"],
         ["GAM", "0.234", "0.008"]]))
A(CAP("Table S5. ", "Brier score by method (lower is better) on the calibration scenarios. Values are the "
      "mean and standard deviation across the six worlds and three replicates."))
A(add_image("figures/FigureS2_calibration.png", "FigureS2_calibration.png", cx_fixed=4114800))
A(CAP("Figure S2. ", "Reliability diagram. For each method the observed frequency of occurrence is plotted "
      "against the mean predicted probability within probability bins; points on the diagonal indicate perfect "
      "calibration. Brier scores are given in the legend."))

# ---------------- S5 ----------------
A(H("S6 | Extended comparison with the geostatistical model"))
A(P("The geostatistical comparator is a Matern Gaussian process classifier over coordinates and covariates, "
    "a tractable representative of the random field approach used by sdmTMB and VAST. Table S6 reports mean "
    "AUC by sample size for stRF and the comparator, together with convergence and median fit time. The "
    "comparator converged in every fit, including at n = 30. stRF held a small advantage at the smallest "
    "samples (+0.025 AUC at n = 30, +0.011 at n = 50), the two methods were indistinguishable at n = 100, and "
    "the comparator drew marginally ahead at the largest samples (by 0.002 to 0.005 AUC). Fit time for the "
    "Gaussian process grew steeply with sample size (from 0.03 s at n = 30 to about 2 s at n = 500), whereas "
    "stRF scaled approximately linearly. The practical case for stRF is therefore operational rather than a "
    "claim of superior accuracy: it needs no mesh, accommodates many covariates without bespoke "
    "specification, and converges reliably, while a geostatistical model remains preferable when formal "
    "parametric uncertainty or a generative spatial model is required."))
A(table([1560, 1760, 1760, 1560, 2720],
        ["n per timestep", "stRF AUC", "Gaussian process AUC", "Convergence", "Median GP fit time (s)"],
        [["30", "0.775", "0.751", "100%", "0.03"],
         ["50", "0.809", "0.798", "100%", "0.04"],
         ["100", "0.853", "0.847", "100%", "0.08"],
         ["200", "0.881", "0.883", "100%", "0.24"],
         ["500", "0.912", "0.917", "100%", "2.07"]]))
A(CAP("Table S6. ", "stRF versus the Gaussian process (sdmTMB-class) comparator by sample size, pooled over "
      "the six worlds, with convergence rate and median Gaussian process fit time."))

# ---------------- S6 ----------------
A(H("S7 | Case study data sources and harmonisation"))
A(P("The case study uses long term Posidonia sinuosa shoot density monitoring from five marine parks along "
    "the Western Australian coast: Jurien Bay, Marmion, Cockburn Sound, Shoalwater Islands and the Ngari "
    "Capes, spanning 2003 to 2024. Raw quadrat counts (17,935 observations) were aggregated to site by year "
    "means, yielding 748 site-years across 74 unique monitoring sites. The response is the natural logarithm "
    "of one plus shoot density per 0.04 m squared."))
A(P("Each site-year carries depth, water clarity (the light attenuation coefficient kd490), a categorical "
    "habitat rank and three summer marine heatwave descriptors (the mean, maximum and coefficient of "
    "variation of summer intensity). Harmonisation involved unifying site identifiers and coordinates across "
    "monitoring programmes, aggregating quadrats to site-years, and reconciling covariate definitions to a "
    "common set. For the prediction maps (Figure 7) the coastline was taken from the geo-maps "
    "countries-land-100m dataset (Primarosa, 2020), derived from Natural Earth 1:10 m physical vectors and "
    "OpenStreetMap, and covariate surfaces along the coastal corridor were interpolated from the sites by "
    "inverse distance weighting."))

# ---------------- S7 ----------------
A(H("S8 | marineSim R package reference"))
A(P("The methods are provided as an R package, marineSim (version 0.3.0), that mirrors the Python "
    "implementation used for the benchmark. The package depends on ranger for random forests, mgcv for "
    "generalised additive models and FNN for nearest neighbour features. The exported functions are "
    "summarised in Table S7; a full simulation run and the case study are driven by the scripts "
    "01_run_simulation.R and 02_case_study.R, and a reproducible pipeline is provided in pipeline/_targets.R."))
A(table([2600, 6760],
        ["Function", "Purpose"],
        [["matern_grf()", "Simulate a Matern Gaussian random field by spectral synthesis."],
         ["make_world()", "Generate a ground-truth world (drivers, occurrence probability and labels)."],
         ["strf_features()", "Build the stRF feature set: multi scale spatial lags, temporal lag and distance fields."],
         ["fit_model()", "Fit any of GLM, GAM, RF, Spatial RF or stRF with a common interface."],
         ["evaluate_auc()", "Compute the area under the ROC curve for held-out predictions."],
         ["decompose_gain()", "Decompose the AUC gain into machine learning, spatial and temporal components."],
         ["area_of_applicability()", "Flag query points outside the training dissimilarity threshold."],
         ["run_scenario()", "Draw one replicate of a scenario and fit and evaluate all methods."]]))
A(CAP("Table S7. ", "Exported functions of the marineSim R package."))

# ---------------- assemble ----------------
new_doc = preamble + "".join(body) + sectpr
open(DOCPATH, "w").write(new_doc)

# rels: keep non-image, add every registered image
rels = open(f"{U}/word/_rels/document.xml.rels").read()
kept = [m.group(0) for m in re.finditer(r'<Relationship [^>]*/>', rels)
        if "/image" not in m.group(0)]
imgrels = "".join(
    f'<Relationship Id="{rid}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/{media}"/>'
    for rid, media, _ in IMAGES)
newrels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
           '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
           + "".join(kept) + imgrels + '</Relationships>')
open(f"{U}/word/_rels/document.xml.rels", "w").write(newrels)

# media: drop main-figure PNGs, copy every registered image
for f in os.listdir(f"{U}/word/media"):
    if f.lower().endswith(".png"):
        os.remove(f"{U}/word/media/{f}")
for rid, media, src in IMAGES:
    shutil.copy(src, f"{U}/word/media/{media}")

# header text -> Supplementary Material
hp = f"{U}/word/header1.xml"
if os.path.exists(hp):
    h = open(hp).read()
    h = h.replace("Spatio-temporal RF for marine SDM \u2014 Afrifa-Yamoah et al.",
                  "Supplementary Material \u2014 Afrifa-Yamoah et al.")
    open(hp, "w").write(h)

print("SI document built. body blocks:", len(body))
print("kept rels:", len(kept), "| media:", os.listdir(f"{U}/word/media"))
