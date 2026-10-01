# =============================================================================
# 05_si_figures_data.R — data behind the simulation SI figures and tables
# =============================================================================
#   Figure S2: six occurrence probability worlds + sampling designs + n gradient
#   Figure S3 / Table S1: hyperparameter sensitivity (3 x 3 grid, RF family)
#   Figure S4 / Table S5: probabilistic calibration (reliability + Brier), incl. GAM
# Writes CSVs to results/si/; the figures are drawn from these files by
# figures/scripts/plot_si_figs.py. Usage: Rscript 05_si_figures_data.R (~30 min)
# =============================================================================
suppressPackageStartupMessages({library(marineSim); library(mgcv)})
OUT_DIR <- Sys.getenv("ABL_OUT_DIR", "../results")
OUT_SI <- file.path(OUT_DIR, "si"); dir.create(OUT_SI, showWarnings = FALSE, recursive = TRUE)
worlds6 <- world_grid(0)             # realisation 0 of each of the six configurations

# ---------------------------------------------------------------- Figure S2 data
rows <- list()
for (wg in worlds6) {
  w <- make_world(wg$stat, wg$range_km, n_time = 5, seed = wg$seed)
  P <- w$p[[5]]
  rows[[length(rows) + 1]] <- data.frame(cfg = wg$cfg, stationary = wg$stat, range_km = wg$range_km,
                                         cell = seq_len(GRID^2), p = as.vector(P))
}
write.csv(do.call(rbind, rows), file.path(OUT_SI, "figS2_worlds.csv"), row.names = FALSE)
# representative world for the design / n panels: non stationary, 80 km
wrep <- worlds6[[which(sapply(worlds6, function(g) g$cfg == "nonstat_r80"))]]
w <- make_world(wrep$stat, wrep$range_km, n_time = 5, seed = wrep$seed)
des <- list()
set.seed(wrep$seed + 11)
for (dsg in c("random", "clustered", "stratified")) {
  s <- sample_sites(100, dsg); xy <- cell_xy(s)
  des[[length(des) + 1]] <- data.frame(panel = paste0("design_", dsg), x = xy[, 1], y = xy[, 2], y_obs = as.vector(w$y[[5]])[s])
}
for (n in c(30, 100, 500)) {
  s <- sample_sites(n, "random"); xy <- cell_xy(s)
  des[[length(des) + 1]] <- data.frame(panel = paste0("n_", n), x = xy[, 1], y = xy[, 2], y_obs = as.vector(w$y[[5]])[s])
}
write.csv(do.call(rbind, des), file.path(OUT_SI, "figS2_designs.csv"), row.names = FALSE)
write.csv(data.frame(cell = seq_len(GRID^2), p = as.vector(w$p[[5]])), file.path(OUT_SI, "figS2_repworld.csv"), row.names = FALSE)
cat("Figure S2 data written\n")

# ---------------------------------------------------------------- tuning helper
fit_auc <- function(ds, subset, mtry_rule, mns, seed) {
  Xtr <- do.call(cbind, c(list(ds$Xc_tr), ds$btr[subset])); Xte <- do.call(cbind, c(list(ds$Xc_te), ds$bte[subset]))
  colnames(Xtr) <- colnames(Xte) <- paste0("f", seq_len(ncol(Xtr))); p <- ncol(Xtr)
  mtry <- switch(mtry_rule, sqrt = max(1, floor(sqrt(p))), half = max(1, floor(0.5 * p)), p8 = max(1, floor(0.8 * p)))
  dat <- data.frame(y = factor(ds$y_tr, levels = c(0, 1)), Xtr)
  fit <- ranger(y ~ ., data = dat, num.trees = RF_TREES, mtry = mtry, min.node.size = mns,
                probability = TRUE, seed = seed, num.threads = 1)
  auc(ds$y_te, predict(fit, data.frame(Xte))$predictions[, "1"])
}
METHODS <- list("Standard RF" = character(0), "Spatial RF" = "COORD", "stRF" = BLOCKS)

# ---------------------------------------------------------------- Figure S3 / Table S1
tun <- list()
for (wg in worlds6) {
  w <- make_world(wg$stat, wg$range_km, n_time = 5, seed = wg$seed)
  for (n in c(50, 100, 200)) {
    ds <- build_dataset(w, n, 5, "random", seed = wg$seed * 7 + n)
    if (is.null(ds)) next
    for (m in names(METHODS)) for (mns in c(1, 3, 5)) for (mr in c("sqrt", "half", "p8"))
      tun[[length(tun) + 1]] <- data.frame(cfg = wg$cfg, n = n, method = m, min_node_size = mns,
                                           mtry_rule = mr, auc = fit_auc(ds, METHODS[[m]], mr, mns, seed = 7))
  }
  cat("tuning:", wg$cfg, "done\n")
}
tun <- do.call(rbind, tun)
write.csv(tun, file.path(OUT_SI, "figS3_tuning.csv"), row.names = FALSE)
agg <- aggregate(auc ~ method + min_node_size + mtry_rule, tun, mean)
tab <- do.call(rbind, lapply(names(METHODS), function(m) {
  a <- agg[agg$method == m, ]; dflt <- a$auc[a$min_node_size == 1 & a$mtry_rule == "sqrt"]
  data.frame(method = m, default_auc = dflt, best_auc = max(a$auc), worst_auc = min(a$auc), range = max(a$auc) - min(a$auc))
}))
write.csv(tab, file.path(OUT_SI, "tableS1_tuning.csv"), row.names = FALSE); print(tab, digits = 3)

# ---------------------------------------------------------------- Figure S4 / Table S5
cal <- list()
for (wg in worlds6) {
  w <- make_world(wg$stat, wg$range_km, n_time = 5, seed = wg$seed)
  for (rep in 0:2) {
    ds <- build_dataset(w, 100, 5, "random", seed = wg$seed * 13 + rep)
    if (is.null(ds)) next
    dtr <- data.frame(y = ds$y_tr, ds$Xc_tr); dte <- data.frame(ds$Xc_te)
    preds <- list()
    preds[["GLM"]] <- predict(suppressWarnings(glm(y ~ ., data = dtr, family = binomial())), dte, type = "response")
    gm <- suppressWarnings(gam(y ~ s(depth, k = 5) + s(subst, k = 5) + s(wave, k = 5) + s(sst, k = 5),
                               data = dtr, family = binomial(), method = "REML"))
    preds[["GAM"]] <- as.vector(predict(gm, dte, type = "response"))
    for (m in names(METHODS)) {
      Xtr <- do.call(cbind, c(list(ds$Xc_tr), ds$btr[METHODS[[m]]])); Xte <- do.call(cbind, c(list(ds$Xc_te), ds$bte[METHODS[[m]]]))
      colnames(Xtr) <- colnames(Xte) <- paste0("f", seq_len(ncol(Xtr)))
      fit <- ranger(y ~ ., data = data.frame(y = factor(ds$y_tr, levels = c(0, 1)), Xtr), num.trees = RF_TREES,
                    mtry = floor(sqrt(ncol(Xtr))), min.node.size = 1, probability = TRUE, seed = 7, num.threads = 1)
      preds[[m]] <- predict(fit, data.frame(Xte))$predictions[, "1"]
    }
    for (m in names(preds))
      cal[[length(cal) + 1]] <- data.frame(cfg = wg$cfg, rep = rep, method = m, p_hat = preds[[m]], y = ds$y_te, p_true = ds$p_te)
  }
  cat("calibration:", wg$cfg, "done\n")
}
cal <- do.call(rbind, cal)
write.csv(cal, file.path(OUT_SI, "figS4_calibration.csv"), row.names = FALSE)
cal$sq <- (cal$p_hat - cal$y)^2
bri <- aggregate(sq ~ method + cfg + rep, cal, mean)
tab5 <- do.call(rbind, lapply(unique(bri$method), function(m) {
  b <- bri$sq[bri$method == m]; data.frame(method = m, brier = mean(b), sd = sd(b)) }))
tab5 <- tab5[order(tab5$brier), ]
write.csv(tab5, file.path(OUT_SI, "tableS5_brier.csv"), row.names = FALSE); print(tab5, digits = 3)
cat("SI figure data COMPLETE\n")
