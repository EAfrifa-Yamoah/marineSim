#!/usr/bin/env Rscript
# Reproduce the simulation benchmark in R (serial reference implementation).
# For the full factorial grid the {targets} pipeline in pipeline/_targets.R is
# preferred; this script is a transparent, dependency light driver.
suppressPackageStartupMessages({
  library(marineSim)
})

set.seed(20260601)
world_grid <- expand.grid(nonstationary = c(FALSE, TRUE),
                          spatial_range_km = c(30, 80, 200))
designs <- c("random", "clustered", "stratified")
sample_sizes <- c(30L, 50L, 100L, 200L, 500L)
n_rep <- 6L

rows <- list(); k <- 1
for (i in seq_len(nrow(world_grid))) {
  w <- make_world(world_id = i - 1L,
                  nonstationary = world_grid$nonstationary[i],
                  spatial_range_km = world_grid$spatial_range_km[i],
                  n_time = 5L)
  for (n in sample_sizes) for (d in designs) for (r in seq_len(n_rep)) {
    res <- run_scenario(w, n = n, design = d)
    res$world_id <- i - 1L; res$rep <- r
    rows[[k]] <- res; k <- k + 1
    message(sprintf("world %d  n=%d  design=%s  rep=%d  stRF AUC=%.3f",
                    i - 1L, n, d, r, res$stRF))
  }
}
benchmark <- do.call(rbind, rows)
dir.create("results", showWarnings = FALSE)
write.csv(benchmark, "results/benchmark_R.csv", row.names = FALSE)

# headline decomposition (cluster averaged over worlds)
agg <- aggregate(cbind(algorithmic, spatial, temporal, total) ~ world_id,
                 data = benchmark, FUN = mean)
print(round(colMeans(agg[, -1]), 3))
