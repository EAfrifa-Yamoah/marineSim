# =============================================================================
# 02b_temporal_coverage.R — temporal value against spatial coverage (Figure 4c)
# =============================================================================
# For every binary, full detection dataset at T = 5 in Stage A, rebuilds the
# dataset (no model fitting) to measure how far the held out cells lie from their
# nearest training site, and relates the temporal Shapley value to that distance.
# Writes results/temporal_coverage_T5.csv (one row per dataset) and
# results/expanded_T8_temporal_by_coverage.csv (binned means with cluster
# bootstrap intervals, by design). Also reports how many clustered sites are
# random fill once the three cluster footprints saturate (Section 3.2).
# Usage: Rscript 02b_temporal_coverage.R   (~5 min; needs ablation_expanded.csv)
# =============================================================================
suppressPackageStartupMessages(library(marineSim))
OUT_DIR <- Sys.getenv("ABL_OUT_DIR", "../results")
d <- read.csv(file.path(OUT_DIR, "ablation_expanded.csv"), check.names = FALSE, stringsAsFactors = FALSE)
d <- add_shapley(d); B <- d[d$response == "binary" & d$detection == 1 & d$n_time == 5 & d$stage == "A_designs", ]
dsg_i <- c(random = 1, clustered = 2, stratified = 3); wg_all <- world_grid(0:5)
key <- sapply(wg_all, function(w) paste(w$cfg, w$real)); wcache <- list(); nn_med <- numeric(nrow(B))
for (i in seq_len(nrow(B))) {
  k <- paste(B$world_cfg[i], B$realisation[i]); wg <- wg_all[[match(k, key)]]
  if (is.null(wcache[[k]])) wcache[[k]] <- make_world(wg$stat, wg$range_km, n_time = 5, seed = wg$seed)
  ds <- build_dataset(wcache[[k]], B$n[i], 5, B$design[i], dataset_seed(wg$seed, B$n[i], 5, dsg_i[[B$design[i]]]))
  nn_med[i] <- median(apply(dist_km(ds$bte$COORD, ds$btr$COORD), 1, min))
}
B$nn_median_km <- nn_med
cov <- B[, c("world_cfg", "realisation", "field", "design", "n", "nn_median_km", "phi_COORD", "phi_EDF", "phi_LAG", "phi_TEMP",
             "spatial", "alg", "total", "temp_marg_withNN0")]
write.csv(cov, file.path(OUT_DIR, "temporal_coverage_T5.csv"), row.names = FALSE)

br <- c(5, 10, 15, 20, 30, 45, 70, 160); B$bin <- cut(B$nn_median_km, br)
t8 <- do.call(rbind, lapply(split(B, list(B$bin, B$design), drop = TRUE), function(x) if (nrow(x) >= 10)
  data.frame(design = x$design[1], bin = as.character(x$bin[1]), x_km = exp(mean(log(x$nn_median_km))), t(boot_ci(x, "phi_TEMP")), n_datasets = nrow(x))))
t8 <- t8[order(t8$design, t8$x_km), ]
write.csv(t8, file.path(OUT_DIR, "expanded_T8_temporal_by_coverage.csv"), row.names = FALSE)
cat("Temporal Shapley (T = 5) by median distance from a test cell to its nearest site\n"); print(t8, digits = 3, row.names = FALSE)

# Random fill in the clustered design: cells drawn uniformly once the clusters saturate
set.seed(1); fill <- NULL
for (n in c(30, 50, 100, 200, 500)) for (r in 1:50) {
  g <- GRID; ctr <- matrix(runif(6, 0.15, 0.85) * g, 3, 2); per <- rep(n %/% 3, 3); per[seq_len(n %% 3)] <- per[seq_len(n %% 3)] + 1
  cells <- integer(0)
  for (k in 1:3) { rr <- round(ctr[k, 1] + rnorm(per[k], 0, g * 0.06)); cc <- round(ctr[k, 2] + rnorm(per[k], 0, g * 0.06)); cells <- c(cells, cell_idx(rr, cc)) }
  fill <- rbind(fill, data.frame(n = n, fill_pct = 100 * (n - length(unique(cells))) / n))
}
f <- aggregate(fill_pct ~ n, fill, mean); cat("\nClustered design: mean percentage of sites placed at random after de-duplication\n"); print(f, digits = 3, row.names = FALSE)
write.csv(f, file.path(OUT_DIR, "clustered_random_fill.csv"), row.names = FALSE)
