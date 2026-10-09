# =============================================================================
# 03_gam_core.R — spatial GAM (mgcv) on the Stage A datasets
# =============================================================================
# s(covariates) + te(x, y), REML, binomial; same worlds, designs and dataset seeds
# as Stage A of 01_ablation_expanded.R so that rows merge one to one. Writes
# results/si/gam_core.csv. Usage: Rscript 03_gam_core.R   (~1 h on one core)
# =============================================================================
suppressPackageStartupMessages({library(marineSim); library(mgcv)})
OUT_DIR <- Sys.getenv("ABL_OUT_DIR", "../results")
OUT_G <- file.path(OUT_DIR, "si", "gam_core.csv"); dir.create(dirname(OUT_G), showWarnings = FALSE, recursive = TRUE)
rows <- list(); DESIGNS <- c("random", "clustered", "stratified")
if (nzchar(Sys.getenv("ABL_DESIGNS"))) DESIGNS <- strsplit(Sys.getenv("ABL_DESIGNS"), ",")[[1]]
if (file.exists(OUT_G) && nzchar(Sys.getenv("ABL_DESIGNS"))) {   # keep rows of the designs not being rerun
  prev <- read.csv(OUT_G, stringsAsFactors = FALSE); rows <- split(prev[!prev$design %in% DESIGNS, ], seq_len(sum(!prev$design %in% DESIGNS))) }
for (wg in world_grid(0:5)) {
  w <- make_world(wg$stat, wg$range_km, n_time = 5, seed = wg$seed)
  for (dsg in DESIGNS) for (n in c(30, 50, 100, 200, 500)) for (nt in c(1, 3, 5)) {
    di <- match(dsg, c("random", "clustered", "stratified"))
    dseed <- dataset_seed(wg$seed, n, nt, di)
    ds <- build_dataset(w, n, nt, dsg, dseed); if (is.null(ds)) next
    dtr <- data.frame(y = ds$y_tr, ds$Xc_tr, X = ds$btr$COORD[, 1], Y = ds$btr$COORD[, 2])
    dte <- data.frame(ds$Xc_te, X = ds$bte$COORD[, 1], Y = ds$bte$COORD[, 2])
    kk <- if (n <= 50) 4 else 5
    a <- tryCatch({ g <- gam(y ~ s(depth, k = kk) + s(subst, k = kk) + s(wave, k = kk) + s(sst, k = kk) + te(X, Y, k = c(kk, kk)),
                              data = dtr, family = binomial(), method = "REML")
                    auc(ds$y_te, predict(g, dte, type = "response")) }, error = function(e) NA)
    rows[[length(rows) + 1]] <- data.frame(world_cfg = wg$cfg, realisation = wg$real, design = dsg, n = n, n_time = nt, GAM = a)
  }
  write.csv(do.call(rbind, rows), OUT_G, row.names = FALSE); cat("GAM:", wg$cfg, "r", wg$real, "|", length(rows), "\n")
}
cat("GAM COMPLETE\n")
