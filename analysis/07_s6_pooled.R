# =============================================================================
# 07_s6_pooled.R — Supplement S6: stRF fitted on all five timesteps' rows (pooled),
# predicting the final timestep, like for like with spatiotemporal sdmTMB.
# Writes results/si/tableS6_strf_pooled.csv. Usage: Rscript 07_s6_pooled.R
# =============================================================================
suppressPackageStartupMessages(library(marineSim))
OUT_DIR <- Sys.getenv("ABL_OUT_DIR", "../results")
build_pooled <- function(w, ds) {
  xy <- ds$btr$COORD; n <- nrow(xy); K <- min(N_ANCHORS, n)
  anchors <- kmeans(xy, centers = K, nstart = 3, iter.max = 50)$centers
  E <- edf(xy, anchors); gm_all <- mean(unlist(ds$obs))
  X <- list(); Y <- list()
  for (t in 1:5) {
    yt <- ds$obs[[t]]; gm <- mean(yt)
    L <- spatial_lags(xy, yt, xy, TRUE)
    if (t >= 2) { prev <- ds$obs[[t - 1]]; hist <- Reduce("+", ds$obs[1:(t - 1)]) / (t - 1)
      j <- nearest_train(xy, xy); TEMP <- cbind(prev[j], hist[j]) } else TEMP <- matrix(gm, n, 2)
    X[[t]] <- cbind(covariates(w, ds$site, t), xy, E, L, TEMP, t = t); Y[[t]] <- yt
  }
  Xtr <- do.call(rbind, X); ytr <- unlist(Y)
  Xte <- cbind(ds$Xc_te, ds$bte$COORD, ds$bte$EDF, ds$bte$LAG, ds$bte$TEMP, t = 5)
  colnames(Xtr) <- colnames(Xte) <- paste0("f", seq_len(ncol(Xtr)))
  list(Xtr = Xtr, ytr = ytr, Xte = Xte)
}
rows <- list()
for (wg in world_grid(0)) {
  w <- make_world(wg$stat, wg$range_km, n_time = 5, seed = wg$seed)
  for (n in c(30, 50, 100, 200, 500)) for (rep in 0:2) {
    ds <- build_dataset(w, n, 5, "random", seed = wg$seed * 17 + n * 3 + rep); if (is.null(ds)) next
    P <- build_pooled(w, ds); p <- ncol(P$Xtr)
    fit <- ranger(y ~ ., data = data.frame(y = factor(P$ytr, levels = c(0, 1)), P$Xtr), num.trees = RF_TREES,
                  mtry = floor(sqrt(p)), min.node.size = 1, probability = TRUE, seed = 7, num.threads = 1)
    a <- auc(ds$y_te, predict(fit, data.frame(P$Xte))$predictions[, "1"])
    rows[[length(rows) + 1]] <- data.frame(cfg = wg$cfg, n = n, rep = rep, stRF_pooled = a)
  }
  cat("pooled:", wg$cfg, "done\n")
}
write.csv(do.call(rbind, rows), file.path(OUT_DIR, "si", "tableS6_strf_pooled.csv"), row.names = FALSE); cat("POOLED COMPLETE\n")
