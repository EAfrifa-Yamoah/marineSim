# marineSim — spatiotemporal random forests for marine SDM under data limitation

Reproducible code for the simulation benchmark, applied seagrass case study and
figures behind the manuscript *Spatio-temporal random forests with adaptive local
features for marine species distribution modelling under data limitation*
(submitted to *Methods in Ecology and Evolution*).

The repository ships both a self contained **R package** (`marineSim/`) that
implements the method and the benchmark, and the **Python pipeline**
(`python/`) used to run the full factorial study and render the figures. The two
implementations follow the same design so results can be cross checked.

---

## What the method does

The spatiotemporal random forest (stRF) keeps the nonparametric learner but
replaces the hope that trees will rediscover spatial structure with an explicit
local feature layer:

1. **Multi scale spatial lags** — distance weighted neighbour means of the
   response at several radii (25, 75 and 150 km in the study), matched within
   time. These carry the broad spatial signal that raw coordinates only encode
   implicitly.
2. **Temporal lags** — the previous time point response at the nearest site, plus
   a basin wide annual mean, so the model can borrow strength across years.
3. **Euclidean distance fields** — distance to a small set of spatial anchors,
   giving the forest a smooth global coordinate system.

The benchmark isolates where predictive skill comes from by walking a nested
ladder: GLM (covariates only) → RF (adds algorithmic flexibility) → SpatialRF
(adds coordinates) → stRF (adds the full feature layer). Differencing the rungs
splits the total gain into algorithmic, spatial and temporal contributions.

## Headline findings

- In the data limited regime (tens of sites), **spatial structure delivers about
  two thirds of the achievable gain** over a non spatial baseline, temporal
  features add roughly a quarter and are most useful with at least three time
  points, and pure algorithmic flexibility is marginal once spatial structure is
  available.
- **stRF reaches useful skill (AUC 0.75) at around 50 sites**, where an aspatial
  random forest needs roughly an order of magnitude more.
- A correctly specified geostatistical comparator matches stRF on accuracy; the
  practical advantage of stRF is operational (no mesh, many covariates, scalable).
- In the seagrass case study, aspatial models fail to transfer across latitude
  while spatial methods recover the gradient, and the **area of applicability**
  correctly flags the unmonitored latitudinal gaps as unsupported extrapolation.

---

## Repository layout

```
marineSim/            R package
  R/                  grf, world, features, fit, evaluate, aoa
  tests/testthat/     unit tests
  DESCRIPTION, NAMESPACE, LICENSE
pipeline/
  _targets.R          reproducible {targets} pipeline for the benchmark
python/               Python pipeline used for the manuscript
  sim_engine.py       calibrated ground-truth engine + model ladder
  benchmark.py        main factorial grid runner
  aux_runs.py         GAM / geostatistical / detection subsets
  aggregate.py        decomposition, cluster bootstrap, mixed model, fig tables
  case_study.py       seagrass latitudinal extrapolation + AOA
  make_figures.py     renders the six manuscript figures
  figstyle.py         shared Okabe-Ito publication style
scripts/
  01_run_simulation.R serial R driver for the benchmark
  02_case_study.R     R driver for the seagrass analysis
results/              regenerated outputs (csv, json)
figures/              regenerated figures (png)
requirements.txt      Python dependencies
```

---

## Reproducing the study

### Python pipeline (full factorial study and figures)

```bash
pip install -r requirements.txt
cd python
# main grid (six worlds), then auxiliary runs
for w in 0 1 2 3 4 5; do python benchmark.py $w; done
for w in 0 1 2 3 4 5; do python aux_runs.py gam   $w $((w+1)); done
for w in 0 1 2 3 4 5; do python aux_runs.py geo   $w $((w+1)); done
for w in 0 1 2 3 4 5; do python aux_runs.py det07 $w $((w+1)); done
python aggregate.py        # -> results/derived.json + figure tables
python case_study.py       # -> results/case_*.{csv,json}   (needs the dataset)
python make_figures.py     # -> figures/Figure1..6.png
```

The seagrass case study expects `Bayesiandataset_2025_final.csv` (the *Posidonia
sinuosa* shoot density dataset). Place it where `case_study.py` and
`scripts/02_case_study.R` expect it, or edit the path at the top of each file.

### R package

```r
# install
install.packages(c("ranger", "mgcv", "FNN"))
# from the repository root
install.packages("marineSim", repos = NULL, type = "source")

library(marineSim)
w   <- make_world(world_id = 0, nonstationary = TRUE, spatial_range_km = 80)
res <- run_scenario(w, n = 100, design = "random")
res[, c("GLM", "RF", "SpatialRF", "stRF", "spatial", "temporal")]
```

The full benchmark can be run with the `{targets}` pipeline:

```r
setwd("pipeline"); targets::tar_make(); targets::tar_read(decomposition)
```

---

## A note on reproducibility and scope

The simulation deliberately makes the environmental covariates rough so they
cannot act as a hidden location proxy; this isolates the contribution of genuine
spatial structure from covariate driven prediction. All randomness is seeded per
world. The benchmark is built for the small sample regime that motivates the
method and should not be read as a claim about behaviour with thousands of sites,
where a correctly specified geostatistical model or a richer learner may be
preferable.

## Citation

If you use this code, please cite the manuscript (in review) and this repository.
A `CITATION.cff` will be added on acceptance.

## License

MIT (see `LICENSE`).

## Coastline data (Figure 7)

The Western Australian coastline drawn in Figure 7 comes from the geo-maps
`countries-land-100m` dataset (Primarosa, 2020), derived from Natural Earth
1:10 m physical vectors and OpenStreetMap. The clipped 45 KB `results/coast_wa.json`
is shipped so the maps render out of the box. To regenerate it:

```bash
npm install @geo-maps/countries-land-100m
python python/extract_coastline.py \
  node_modules/@geo-maps/countries-land-100m/map.geo.json results/coast_wa.json
```

In the manuscript workflow the polygons are processed with the sf package in R
(Pebesma, 2018); `extract_coastline.py` provides an equivalent dependency-free
path in Python.
