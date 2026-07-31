#!/usr/bin/env python3
"""
finish.py — complete the manuscript and supplement update.

  1. Swap Figure 2 and Figure 6 images and resize their frames.
  2. Swap the supplement forecast-skill figure; add a new coverage figure.
  3. Insert new main-text Results subsections 5.8 and 5.9.
  4. Renumber the duplicated 5.6 heading.
  5. Update the case study Results numbers.
  6. Correct Table 1 and the residual Table 3 share.
  7. Add the three new references and cite Webster et al. (2024).
"""
import re, shutil, os

MSD = "/home/claude/work/ms"
SUPPD = "/home/claude/work/supp"
RES = "/home/claude/work/repo/results"
EMU = 914400


def txt_para(text, style_pPr):
    """A body paragraph carrying `text`, using the supplied pPr block."""
    rpr = ('<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"'
           ' w:cs="Times New Roman"/></w:rPr>')
    esc = (text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
    return (f'<w:p>{style_pPr}<w:r>{rpr}'
            f'<w:t xml:space="preserve">{esc}</w:t></w:r></w:p>')


def heading_para(text, pPr):
    esc = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    rpr = ('<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"'
           ' w:cs="Times New Roman"/><w:b/><w:bCs/></w:rPr>')
    return f'<w:p>{pPr}<w:r>{rpr}<w:t xml:space="preserve">{esc}</w:t></w:r></w:p>'


def drawing(rid, cx, cy, name):
    return (
        '<w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:r><w:drawing>'
        f'<wp:inline distT="0" distB="0" distL="0" distR="0">'
        f'<wp:extent cx="{cx}" cy="{cy}"/>'
        '<wp:effectExtent l="0" t="0" r="0" b="0"/>'
        f'<wp:docPr id="{9000 + abs(hash(name)) % 900}" name="{name}"/>'
        '<wp:cNvGraphicFramePr><a:graphicFrameLocks '
        'xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" '
        'noChangeAspect="1"/></wp:cNvGraphicFramePr>'
        '<a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
        '<a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">'
        '<pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">'
        f'<pic:nvPicPr><pic:cNvPr id="0" name="{name}"/><pic:cNvPicPr/></pic:nvPicPr>'
        f'<pic:blipFill><a:blip r:embed="{rid}"/><a:stretch><a:fillRect/>'
        '</a:stretch></pic:blipFill>'
        f'<pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{cx}" cy="{cy}"/></a:xfrm>'
        '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr>'
        '</pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>')


def add_image(docdir, png, target_name):
    """Copy png into media/ and register a relationship; return its rId."""
    shutil.copy(png, f"{docdir}/word/media/{target_name}")
    rp = f"{docdir}/word/_rels/document.xml.rels"
    s = open(rp, encoding="utf-8").read()
    ids = [int(x) for x in re.findall(r'Id="rId(\d+)"', s)]
    new = f"rId{max(ids) + 1}"
    rel = (f'<Relationship Id="{new}" Type="http://schemas.openxmlformats.org/'
           f'officeDocument/2006/relationships/image" Target="media/{target_name}"/>')
    s = s.replace("</Relationships>", rel + "</Relationships>")
    open(rp, "w", encoding="utf-8").write(s)
    return new


def resize_before(s, anchor, cx, cy):
    """Resize the last wp:extent and a:ext occurring before `anchor`."""
    i = s.find(anchor)
    if i < 0:
        return s, False
    seg = s[:i]
    m = list(re.finditer(r'<wp:extent cx="\d+" cy="\d+"/>', seg))
    if not m:
        return s, False
    last = m[-1]
    s = s[:last.start()] + f'<wp:extent cx="{cx}" cy="{cy}"/>' + s[last.end():]
    seg2 = s[:s.find(anchor)]
    m2 = list(re.finditer(r'<a:ext cx="\d+" cy="\d+"/>', seg2))
    if m2:
        last2 = m2[-1]
        s = s[:last2.start()] + f'<a:ext cx="{cx}" cy="{cy}"/>' + s[last2.end():]
    return s, True


# ============================================================== MANUSCRIPT
def do_manuscript():
    p = f"{MSD}/word/document.xml"
    s = open(p, encoding="utf-8").read()
    log = []

    # -- 1. swap figure images ------------------------------------------
    shutil.copy(f"{RES}/FIG_decomposition.png", f"{MSD}/word/media/image1.png")
    shutil.copy(f"{RES}/FIG_case_analysis.png", f"{MSD}/word/media/image5.png")
    s, ok = resize_before(s, "Decomposition of the predictive gain",
                          int(6.5 * EMU), int(3.55 * EMU))
    log.append(("Figure 2 frame", ok))
    s, ok = resize_before(s, "Empirical behaviour of the methods",
                          int(6.9 * EMU), int(2.30 * EMU))
    log.append(("Figure 6 frame", ok))

    # -- 2. residual Table 3 share --------------------------------------
    n = s.count("<w:t>9%</w:t>")
    s = s.replace("<w:t>9%</w:t>", "<w:t>13.4%</w:t>")
    log.append(("Table 3 ML share", n == 1))

    # -- 3. Table 1 corrections -----------------------------------------
    for old, new, lab in [
        ("<w:t>0.7, 1.0</w:t>",
         "<w:t>0.7, 1.0 (separate reduced grid; not crossed with the core factorial)</w:t>",
         "Table 1 detection values"),
        ("<w:t xml:space=\"preserve\">GLM, GAM, RF, Spatial RF, stRF (plus a Gaussian process, sdmTMB-class, supplementary)</w:t>",
         "<w:t xml:space=\"preserve\">Core grid: GLM, RF, Spatial RF, stRF. Auxiliary grids: GAM, and a Matern Gaussian process comparator</w:t>",
         "Table 1 methods values"),
        ("Combinations are fully crossed, and each scenario was replicated six times with different random seeds; the final column summarises what each factor is designed to probe.",
         "Ground-truth world, sample size, sampling design and number of timesteps are fully crossed, giving 270 core scenarios, each replicated six times (1,620 datasets). Detection probability, the generalised additive model and the geostatistical comparator were run on separate reduced grids. Across all runs the benchmark comprises 3,792 datasets and 10,068 model fits.",
         "Table 1 caption"),
    ]:
        c = s.count(old)
        if c == 1:
            s = s.replace(old, new)
        log.append((lab, c == 1))

    # -- 4. case study Results numbers ----------------------------------
    for old, new, lab in [
        ("Under random cross-validation every method predicted shoot density well, with Spearman rank correlations of ",
         "Under random ten-fold cross-validation every method predicted shoot density well, with Spearman rank correlations of 0.75 (standard RF), 0.80 (Spatial RF) and 0.81 (stRF). A site-grouped design, in which whole sites but not whole parks are withheld, gives an intermediate picture (0.49, 0.58 and 0.60), so skill degrades gradually rather than at a cliff edge. Earlier reported values were ",
         "5.6 interpolation numbers"),
        ("Only the spatial methods retained skill, with rank correlations of ",
         "Averaged within held-out parks the spatial methods retained limited skill (0.31 for Spatial RF and 0.35 for stRF) against 0.24 for the GLM and 0.21 for the standard forest, so the aspatial collapse reported when predictions are pooled across parks is partly an artefact of pooling between-park with within-park variation. Both summaries are reported here. The pooled values were ",
         "5.6 LORO numbers"),
        ("Absolute calibration nonetheless remained imperfect: the pooled coefficient of determination was close to zero, reflecting the genuine difficulty of predicting into a latitudinal band that was never sampled.",
         "Absolute calibration failed in every case. The coefficient of determination was negative for all methods (-0.35 for stRF, -0.36 for Spatial RF, -1.28 for the GLM and -1.32 for the standard forest), meaning that none improved on simply predicting the training mean. Rank transfer and calibrated transfer are therefore quite different things here, and only the former survives.",
         "5.6 R2 statement"),
    ]:
        c = s.count(old)
        if c == 1:
            s = s.replace(old, new)
        log.append((lab, c == 1))

    # -- 5. renumber duplicate 5.6 --------------------------------------
    c = s.count("5.6 | Case study: Temporal evolution")
    s = s.replace("5.6 | Case study: Temporal evolution",
                  "5.7 | Case study: temporal reconstruction")
    log.append(("renumber 5.6 -> 5.7", c == 1))

    # -- 6. insert new subsections before the Discussion ----------------
    i = s.find("6 | Discussion")
    hstart = s.rfind("<w:p ", 0, i)
    hpPr = re.search(r'<w:pPr>.*?</w:pPr>', s[hstart:hstart + 3000], re.S)
    hpPr = hpPr.group(0) if hpPr else "<w:pPr></w:pPr>"
    bi = s.find("The benchmark provides a practical decision framework")
    bstart = s.rfind("<w:p ", 0, bi)
    bpPr = re.search(r'<w:pPr>.*?</w:pPr>', s[bstart:bstart + 3000], re.S)
    bpPr = bpPr.group(0) if bpPr else "<w:pPr></w:pPr>"

    S58 = (
        "Under rolling-origin evaluation the framework was outperformed by both "
        "naive baselines at every horizon tested (Figure S6). At one year ahead the "
        "rank correlation with observed density was 0.742 for stRF against 0.770 for "
        "persistence, a site's last observed value, and 0.758 for a site-level "
        "climatology; the gap widened with horizon, reaching 0.630 against 0.693 and "
        "0.695 at five years. Skill at predicting departures from a site's own mean "
        "was 0.146 at one year and fell to between 0.014 and 0.078 thereafter. Shoot "
        "density in this system is strongly persistent and the annual thermal "
        "covariates carry little anticipatory information about heatwave-driven "
        "decline, so a site's own history is already close to the best available "
        "predictor. This bounds the claim made for the method: the engineered "
        "features support spatial prediction and retrospective reconstruction within "
        "the area of applicability, not forecasting.")
    S59 = (
        "The shaded band accompanying ensemble predictions is the spread across trees "
        "and is not a calibrated predictive interval. Against a nominal 80%, its "
        "empirical coverage under random cross-validation was 0.743, while a quantile "
        "regression forest achieved 0.832 and split conformal prediction 0.826 "
        "(Figure S7). Under leave-one-region-out none of the three constructions "
        "reached nominal coverage, at 0.628, 0.680 and 0.472 respectively, and 65.3% "
        "of held-out site-years fell outside the area of applicability. Principled "
        "interval construction therefore repairs uncertainty quantification for "
        "interpolation but not for transfer beyond the sampled latitudes. This is the "
        "central practical argument for reporting an applicability diagnostic "
        "alongside predictions: where prediction is unsupported the remedy is not a "
        "wider interval but an explicit refusal to predict.")

    block = (heading_para("5.8 | Case study: forecast skill", hpPr)
             + txt_para(S58, bpPr)
             + heading_para("5.9 | Case study: predictive uncertainty", hpPr)
             + txt_para(S59, bpPr))
    s = s[:hstart] + block + s[hstart:]
    log.append(("insert 5.8 and 5.9", True))

    # -- 7. order-dependence paragraph in 5.1 ---------------------------
    anchor = "Temporal features (stRF minus Spatial RF) added a further"
    ai = s.find(anchor)
    if ai > 0:
        pend = s.find("</w:p>", ai) + 6
        OD = (
            "Because the engineered blocks encode overlapping information, the "
            "increment attributed to any one of them depends on the order in which "
            "blocks enter. Adding the temporal block to a covariate-only forest raises "
            "test AUC by +0.091; adding it to a model that already carries coordinates, "
            "distance fields and spatial lags raises it by +0.007; and the nested "
            "contrast against a coordinate-only spatial forest gives +0.052. The spread "
            "across orderings exceeds an order of magnitude, so no single nested ladder "
            "identifies the contribution. A placebo test makes the consequence concrete. "
            "At a single timestep the temporal features are constant by construction and "
            "can carry no information, yet the nested contrast still attributes +0.047 "
            "AUC to them, whereas the isolated increment returns +0.001 (Figure 2b). We "
            "therefore report Shapley values over all sixteen subsets of the four blocks "
            "as the primary decomposition, and report the added-last increment separately "
            "because it answers the question a practitioner actually faces: whether to "
            "build temporal features given a spatial model already in hand.")
        s = s[:pend] + txt_para(OD, bpPr) + s[pend:]
        log.append(("insert order-dependence paragraph", True))
    else:
        log.append(("insert order-dependence paragraph", False))

    # -- 8. references and Webster citation -----------------------------
    for old, new, lab in [
        ("Wood, S.N. (2017).",
         "Valavi, R., Elith, J., Lahoz-Monfort, J. J., & Guillera-Arroita, G. (2019). "
         "blockCV: An R package for generating spatially or environmentally separated "
         "folds for k-fold cross-validation of species distribution models. Methods in "
         "Ecology and Evolution, 10(2), 225-232. "
         "MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Royle, J. A., & "
         "Langtimm, C. A. (2002). Estimating site occupancy rates when detection "
         "probabilities are less than one. Ecology, 83(8), 2248-2255. "
         "Thorson, J. T., et al. (2025). tinyVAST: R package with an expressive "
         "interface to specify lagged and simultaneous effects in multivariate "
         "spatio-temporal models. Global Ecology and Biogeography, in press. "
         "Wood, S.N. (2017).", "add 3 references"),
        ("the data comprised 748 site years at 61 sites observed between 2003 and 2024.",
         "the data comprised 748 site years at 61 named sites observed between 2003 and "
         "2024 (Webster et al. 2024); the 61 sites carry 74 distinct coordinate pairs "
         "because some were re-positioned between survey rounds, which matters because "
         "the spatial features key off coordinates rather than site labels.",
         "cite Webster and site count"),
    ]:
        c = s.count(old)
        if c == 1:
            s = s.replace(old, new)
        log.append((lab, c == 1))

    open(p, "w", encoding="utf-8").write(s)
    return log


# ============================================================== SUPPLEMENT
def do_supplement():
    p = f"{SUPPD}/word/document.xml"
    s = open(p, encoding="utf-8").read()
    log = []

    shutil.copy(f"{RES}/FIG_forecast_skill.png", f"{SUPPD}/word/media/image6.png")
    s, ok = resize_before(s, "Temporal forecast skill for the stRF case study",
                          int(6.5 * EMU), int(2.90 * EMU))
    log.append(("Figure S6 frame", ok))

    for old, new, lab in [
        ("Under rolling-origin evaluation, stRF matches but does not exceed the naive baselines for predicting density levels (Figure S6a): at one year ahead the three predictors are within 0.02 in Spearman correlation, and at longer horizons stRF falls below a site-level climatology as its lagged density feature reverts to a regional mean.",
         "Under rolling-origin evaluation stRF is outperformed by both naive baselines at every horizon (Figure S6a). At one year ahead the rank correlation with observed density is 0.742 for stRF, 0.770 for persistence and 0.758 for a site-level climatology; by five years the values are 0.630, 0.693 and 0.695. The lagged density feature reverts toward a regional mean as the horizon lengthens, which is why the gap widens.",
         "S forecast levels"),
        ("stRF has essentially no skill at predicting departures from a site\u2019s mean, the correlation between predicted and observed anomaly is +0.01 at one year ahead and becomes negative at longer horizons, and on the site-years that declined by more than 20% it predicted a small increase on average.",
         "stRF has little skill at predicting departures from a site's mean: the correlation between predicted and observed anomaly is 0.146 at one year ahead and falls to between 0.014 and 0.078 at longer horizons.",
         "S forecast anomaly"),
    ]:
        c = s.count(old)
        if c == 1:
            s = s.replace(old, new)
        log.append((lab, c == 1))

    # new coverage figure at the end of S7
    rid = add_image(SUPPD, f"{RES}/FIG_interval_coverage.png",
                    "image_coverage.png")
    anchor = "Temporal forecast skill for the stRF case study"
    ai = s.find(anchor)
    pend = s.find("</w:p>", ai) + 6 if ai > 0 else -1
    if pend > 0:
        bi = s.find("The case study uses long term")
        bstart = s.rfind("<w:p ", 0, bi)
        bpPr = re.search(r'<w:pPr>.*?</w:pPr>', s[bstart:bstart + 3000], re.S)
        bpPr = bpPr.group(0) if bpPr else "<w:pPr></w:pPr>"
        cap = ("Figure S7. Empirical coverage and width of 80% prediction intervals "
               "for the seagrass case study. (a) Coverage of held-out observations "
               "for the across-tree ensemble spread, a quantile regression forest and "
               "split conformal prediction, under random ten-fold cross-validation and "
               "under leave-one-region-out; the dashed line marks nominal coverage. "
               "(b) Mean interval width on the same scale. No construction attains "
               "nominal coverage once prediction moves beyond the sampled latitudes.")
        s = (s[:pend]
             + drawing(rid, int(6.5 * EMU), int(2.90 * EMU), "FigureS7")
             + txt_para(cap, bpPr) + s[pend:])
        log.append(("insert Figure S7", True))
    else:
        log.append(("insert Figure S7", False))

    open(p, "w", encoding="utf-8").write(s)
    return log


if __name__ == "__main__":
    for lab, lg in [("MANUSCRIPT", do_manuscript()), ("SUPPLEMENT", do_supplement())]:
        print(f"\n### {lab}")
        for name, okf in lg:
            print(f"  {'OK  ' if okf else 'MISS'}  {name}")
