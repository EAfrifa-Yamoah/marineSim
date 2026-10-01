# marineSim — spatio-temporal random forests for marine SDM under data limitation

Code, results and figures behind the article

> Afrifa-Yamoah, E., Fouedjio, F., Mueller, U., McMahon, K. and Hyndes, G. A.
> *Spatio-temporal random forests with adaptive local features for marine species
> distribution modelling under data limitation.* Methods in Ecology and Evolution
> (submitted).

Everything reported in the article and its Supplement is produced in **R** by the
`marineSim` package and the scripts in `analysis/`. The figure files are rendered
from the stored CSV results by short Python (matplotlib) scripts in
`figures/scripts/`; no number is computed in Python.

## What the method does

stRF is an ordinary random forest fitted on an augmented feature map
Φ = [covariates | C | E | L | T] built from the training sites alone:

| Block | Name | Content |
|---|---|---|
| C | Coordinates | x, y |
| E | Euclidean distance fields | distance to k-means anchors of the training sites (15 in the simulation, 8 in the case study) |
| L | Multi scale spatial lags | mean training response within 25, 75 and 150 km, leave one out at training sites |
| T | Temporal features | previous timestep response and site history mean, read at the nearest training site |

The benchmark attributes the gain of stRF over a covariate only GLM to algorithmic
flexibility (RF on covariates minus GLM) and to each block by an exact Shapley
decomposition over the 16 coalitions of {C, E, L, T}, with a placebo (T = 1), a
same timestep nearest neighbour control (NN0) and a cluster bootstrap over
independent field realisations. A standard RF is the empty coalition, a spatial RF
is {C} and stRF is {C, E, L, T}.

## Repository layout

```
marineSim/              R package: the simulation engine (install this first)
  R/                    grf, world, designs, features, dataset, fit, coalitions, aoa
  tests/testthat/       the verification protocol as unit tests
analysis/               scripts that produce every number in the article, in order
  00_verify_engine.R    PASS/FAIL mechanism checks (run first)
  01_ablation_expanded.R  GLM + 20 coalitions on the expanded grid (resumable, shardable)
  02_decompose.R        Shapley tables T1 to T5 with cluster bootstrap
  03_gam_core.R         spatial GAM on the Stage A datasets
  04_core_benchmark_table.R  five method core benchmark (Figure 3, Table 2)
  05_si_figures_data.R  worlds, tuning sensitivity, calibration (S2 to S5)
  06_s6_sdmtmb.R        sdmTMB comparator (S6)
  07_s6_pooled.R        stRF on pooled timesteps (S6)
  08_case_study.R       Posidonia sinuosa case study (Section 4, S7 to S9)
results/                outputs of the scripts above, as reported
  ablation_expanded.csv       3,779 datasets x (GLM + 20 coalitions)
  expanded_T1..T5_*.csv       decomposition tables
  core_five_methods.csv       Figure 3 / Table 2 data
  si/                         gam_core, figS2..S4 data, tableS1/S5/S6
  case_study/                 cv_regimes, importance, corridor, hovmoller, rolling_origin, intervals, meta.json
figures/
  main/Figure2..8.png   article figures (Figure 1 is a schematic drawn outside R)
  si/FigureS1..S8.png   Supplement figures
  scripts/              matplotlib scripts that draw the figures from results/
data/
  coast_wa.json         simplified Western Australian coastline used for the corridor maps
  raw/                  (not distributed) case study data, see Data availability
```

## Reproducing the study

### 1. Install

R >= 4.1 with `ranger` (the only hard dependency of the package). The auxiliary
scripts also use `mgcv`, `jsonlite`, `sp` and, for Supplement S6, `sdmTMB`.

```r
install.packages(c("ranger", "mgcv", "jsonlite", "sp", "testthat"))
install.packages("marineSim", repos = NULL, type = "source")   # from the repository root
# or: remotes::install_github("GITHUB-USER/marineSim", subdir = "marineSim")
```

The reported results were produced with R 4.3.3, ranger 0.16.0 and mgcv 1.9-1 on
x86_64 Linux; sdmTMB was installed from its GitHub source at the time of the S6 run
and its version string was not captured in the archived outputs (record it with
`packageVersion("sdmTMB")` when S6 is rerun). Forest fits are seeded and single threaded
(`num.threads = 1`), so results are reproducible on the same ranger version;
a different ranger version may change the forest's random stream.

### 2. Verify the engine (two minutes)

```bash
cd analysis
Rscript 00_verify_engine.R        # every line must read PASS
```

The same checks run as `testthat` tests inside the package (`R CMD check`).

### 3. Run the simulation benchmark

Single core, about 4 to 5 hours:

```bash
Rscript 01_ablation_expanded.R    # stages A_designs B_T10 C_detect D_abundance
```

In parallel, one shard per core, writing separate files that are merged afterwards:

```bash
for k in 0 1 2 3; do ABL_SHARD=$k/4 Rscript 01_ablation_expanded.R & done; wait
Rscript -e 'f <- list.files("../results", "^ablation_expanded_shard.*csv$", full.names = TRUE)
  d <- do.call(rbind, lapply(f, read.csv, check.names = FALSE))
  key <- c("response","detection","design","world_cfg","realisation","n","n_time","rep")
  write.csv(d[!duplicated(d[key]), ], "../results/ablation_expanded.csv", row.names = FALSE)'
```

Both forms are resumable (completed dataset keys are skipped) and give identical
rows, because every dataset seed is a deterministic function of the world seed,
sample size, timesteps, design, replicate, detection and response type
(`dataset_seed()`).

Then:

```bash
Rscript 02_decompose.R            # tables T1 to T5            (minutes)
Rscript 03_gam_core.R             # spatial GAM                 (~1 h)
Rscript 04_core_benchmark_table.R # Figure 3 data, Table 2 columns
Rscript 05_si_figures_data.R      # S2 to S5 data              (~30 min)
Rscript 07_s6_pooled.R            # S6 pooled stRF             (minutes)
Rscript 06_s6_sdmtmb.R            # S6 sdmTMB (needs sdmTMB)   (~2 h)
```

### 4. Run the case study

Place the site level monitoring data at `data/raw/Bayesiandataset_2025_final.csv`
(see Data availability), then

```bash
Rscript 08_case_study.R           # ~5 min on one core
```

### 5. Draw the figures

```bash
cd ../figures/scripts
python3 make_figures_2_4_5.py; python3 make_figure_3.py
python3 make_figures_S2_S3_S4.py; python3 make_figure_S8.py; python3 make_figure_S1.py
python3 make_figures_case.py      # Figures 6 to 8 and S5 to S7
```

Requires `matplotlib`, `pandas` and `numpy`. Running these on the stored results
regenerates the committed PNG files pixel for pixel.

## Verification performed before release

* `R CMD check` passes (one NOTE: sdmTMB, a suggested package, was not installed on the checking machine).
* A 131 dataset slice of Stage A rerun through the installed package reproduces the
  archived `ablation_expanded.csv` rows exactly (21 score columns, maximum absolute
  difference 0).
* `02_decompose.R` and `04_core_benchmark_table.R` reproduce the archived tables
  exactly from `ablation_expanded.csv`.
* The case study rerun from `08_case_study.R` reproduces every archived CSV in `results/case_study/` exactly.
* Twelve GAM fits rerun through `03_gam_core.R` match `si/gam_core.csv` to 1e-9.
* All fifteen figure files regenerate pixel identical from `results/`.

## Data availability

The simulation is fully synthetic and regenerates from the seeds in `world_grid()`.
The case study uses long term *Posidonia sinuosa* monitoring data from five Western
Australian marine parks that are not ours to redistribute; they are available from
the data custodians on request (see the article's Data Availability Statement).
`data/raw/` is excluded by `.gitignore`.

## Citation

See `CITATION.cff`. Please cite both the article and the archived software release
(Zenodo DOI to be added on release).

## Licence

MIT, see `LICENSE`.
