import re, shutil
from PIL import Image

U = "si_edit"; P = f"{U}/word/document.xml"
doc = open(P).read()
tpl = open("/tmp/si_img_tpl.xml").read()

def esc(s): return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

def H(t):
    return ('<w:p><w:pPr><w:spacing w:after="120" w:line="276" w:before="240"/>'
            '<w:jc w:val="left"/></w:pPr><w:r><w:rPr><w:b/><w:bCs/>'
            '<w:sz w:val="26"/><w:szCs w:val="26"/></w:rPr>'
            f'<w:t xml:space="preserve">{esc(t)}</w:t></w:r></w:p>')

def Pp(t):
    return ('<w:p><w:pPr><w:spacing w:after="140" w:line="276"/><w:jc w:val="both"/></w:pPr>'
            f'<w:r><w:t xml:space="preserve">{esc(t)}</w:t></w:r></w:p>')

def CAP(label, t):
    return ('<w:p><w:pPr><w:spacing w:after="140" w:line="276"/><w:jc w:val="both"/></w:pPr>'
            f'<w:r><w:rPr><w:b/><w:bCs/></w:rPr><w:t xml:space="preserve">{esc(label)}</w:t></w:r>'
            f'<w:r><w:t xml:space="preserve">{esc(t)}</w:t></w:r></w:p>')

IMAGES = []
def IMG(src, media, rid, docpr, cx=5120640):
    w, h = Image.open(src).size
    cy = round(cx * h / w)
    IMAGES.append((rid, media, src))
    x = (tpl.replace('r:embed="rId200"', f'r:embed="{rid}"')
            .replace('wp:docPr id="200"', f'wp:docPr id="{docpr}"')
            .replace('name="FigureS1_tuning.png"', f'name="{media}"')
            .replace('title="FigureS1_tuning.png"', f'title="{media}"'))
    x = re.sub(r'descr="[^"]*"', f'descr="{media}"', x)
    x = re.sub(r'cx="\d+" cy="\d+"', f'cx="{cx}" cy="{cy}"', x)
    return x

M2 = "\u00b2"; MIN = "\u2212"
body = []
A = body.append
A(H("S9 | Temporal evolution of the case study"))
A(Pp("The case study in the main text maps a single year. This section examines how the fitted "
     "stRF surface evolves across the two decades of monitoring, whether that evolution is supported "
     "by the data, and whether it reflects genuine forecasting skill."))

# ---- Figure S3
A(IMG("figures/FigureT1_spatiotemporal_reconstruction.png", "FigureS3_spatiotemporal.png",
      "rId202", 202))
A(CAP("Figure S3. ",
      "Spatiotemporal behaviour of the stRF case study, 2003 to 2024. (a) Regional mean shoot density: "
      "observed values (points) against the stRF reconstruction (lines) obtained under site-grouped "
      "cross-validation, so each park's trajectory is predicted from the other sites. The five parks "
      "diverge, with Jurien Bay, Cockburn Sound and Shoalwater declining and Marmion and Ngari Capes "
      "stable to rising. (b) Predicted density along the coastal corridor as a function of latitude and "
      "year, the corridor collapsed into 70 latitude bands. The model is fitted once on all site-years "
      "and the coast is then predicted for each year using that year's water clarity and summer "
      "heatwave covariates (interpolated from the year's site observations) and that year's density "
      "memory (the previous-year regional mean). Observed monitoring events are overlaid as circles on "
      "the same colour scale. (c) Area of applicability for the same latitude-year field: the "
      "importance-weighted dissimilarity index of Meyer and Pebesma (2021), with the black contour "
      "marking the 95th-percentile training threshold; cells outside the contour are extrapolation, "
      "where the predictions in (b) are unsupported. North is at the top throughout; park positions are "
      "on the left axis and latitude on the right."))
A(Pp("The reconstruction recovers the spatiotemporal structure of the meadow network. Predicted "
     f"density agrees with the monitoring records across the corridor (Spearman 0.77, mean bias {MIN}0.7 "
     f"shoots per 0.04 m{M2}), and the alongshore pattern of change is reproduced in the correct latitude "
     "bands: the northern and central parks fade over the record while the southern capes hold or "
     "strengthen (Figure S3a, b). The magnitude of change is understated, however. The steepest observed "
     "declines are compressed toward each site's long-term mean, the same shrinkage documented for the "
     "static prediction surface, whereby the forest reproduces roughly 70% of the between-site spread, "
     "so the field in panel (b) fades gently rather than tracing the full observed collapse (for example, "
     "a predicted Jurien Bay decline of about five shoots per 0.04 m" + M2 + " against an observed decline "
     "near twelve)."))
A(Pp("The area of applicability qualifies where this field can be trusted (Figure S3c). We compute the "
     "dissimilarity index exactly as for the static map, in the standardised, importance-weighted space "
     "of the environmental covariates and coordinates, but reference it to the training site-years rather "
     "than to site means so that the year's environment enters the comparison. Because a site's repeated "
     "years share the same coordinates and depth, the highest-importance predictors, their mutual "
     "distances are near zero and would collapse the training distance scale; the threshold is therefore "
     "calibrated on between-site nearest-neighbour distances, the temporal analogue of the site-mean "
     "referencing used for the static analysis. On this basis about a third of the corridor-year "
     "combinations fall outside the area of applicability. The flagged cells are informative rather than "
     "incidental: the inter-park gaps are flagged far more often than the parks themselves (42% against "
     "19% of cells), predictions at a park are better supported in years it was surveyed than in years it "
     "was not (for Shoalwater, a mean index of 1.6 against 2.2), and extrapolation increases in the recent "
     "marine-heatwave years (from roughly a third of the coast in the mid-record to over 40% in 2020 to "
     "2024), when the interpolated summer-temperature covariates move beyond the sampled envelope. Much of "
     "the smooth fill in the inter-park gaps and the recent years in panel (b) therefore lies outside the "
     "area of applicability, so the muted trends there reflect genuine extrapolation and not merely "
     "shrinkage."))
A(Pp("Panel (b) is thus a reconstruction of the fitted field, bounded by panel (c), rather than a "
     "forecast, and it should be read as the method representing the observed dynamics within its support, "
     "not anticipating them. Formal forecast skill is assessed next."))

# ---- Figure S4
A(IMG("figures/FigureT2_forecast_skill.png", "FigureS4_forecast_skill.png", "rId203", 203))
A(CAP("Figure S4. ",
      "Temporal forecast skill for the stRF case study. (a) Skill, the Spearman correlation between "
      "predicted and observed density pooled across origins, as a function of forecast horizon, for stRF "
      "and two naive baselines: persistence (a site's last observed value) and climatology (a site's "
      "historical mean). Horizon 0 is a blocked leave-one-year-out evaluation (temporal interpolation, "
      "which may use later years); horizons 1 to 5 are rolling-origin forecasts trained only on years up "
      "to the origin. (b) Skill at predicting departures from a site's mean (the anomaly), the correlation "
      "between predicted and observed anomaly by horizon; values near or below zero indicate no skill at "
      "anticipating change. Persistence and climatology predict near-zero anomaly by construction and are "
      "omitted from (b)."))
A(Pp("Under rolling-origin evaluation, stRF matches but does not exceed the naive baselines for "
     "predicting density levels (Figure S4a): at one year ahead the three predictors are within 0.02 in "
     "Spearman correlation, and at longer horizons stRF falls below a site-level climatology as its lagged "
     "density feature reverts to a regional mean. The decomposition into levels and change is more "
     "revealing (Figure S4b). stRF has essentially no skill at predicting departures from a site's mean, "
     "the correlation between predicted and observed anomaly is +0.01 at one year ahead and becomes "
     "negative at longer horizons, and on the site-years that declined by more than 20% it predicted a "
     "small increase on average. The marine-heatwave covariates therefore carry little anticipatory "
     "information about the observed declines at the annual scale. This locates the temporal value of the "
     "method in the controlled simulation, where the temporal signal is non-trivial by construction and "
     "the temporal features add measurable skill (Section 4), rather than in forecasting this strongly "
     "persistent system, in which a site's historical mean is already a very strong predictor."))

s9 = "".join(body)
sp = doc.rfind("<w:sectPr")
ps = doc.rfind("<w:p>", 0, sp)   # keep any trailing empty paragraph after our block
doc = doc[:sp] + s9 + doc[sp:]
open(P, "w").write(doc)

# rels + media
rels = open(f"{U}/word/_rels/document.xml.rels").read()
add = "".join(
    f'<Relationship Id="{rid}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/{media}"/>'
    for rid, media, _ in IMAGES)
rels = rels.replace("</Relationships>", add + "</Relationships>")
open(f"{U}/word/_rels/document.xml.rels", "w").write(rels)
for rid, media, src in IMAGES:
    shutil.copy(src, f"{U}/word/media/{media}")

print("S9 added with Figures S3, S4:", [m for _, m, _ in IMAGES])
print("body paragraphs added:", len(body))
