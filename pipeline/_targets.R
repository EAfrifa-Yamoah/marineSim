# Reproducible {targets} pipeline for the stRF simulation benchmark.
# Run with: targets::tar_make()   (from within the pipeline/ directory)
library(targets)
tar_option_set(packages = c("marineSim"))

# the six ground-truth worlds: {stationary, non-stationary} x range {30, 80, 200} km
world_grid <- expand.grid(
  nonstationary = c(FALSE, TRUE),
  spatial_range_km = c(30, 80, 200)
)

list(
  tar_target(designs, c("random", "clustered", "stratified")),
  tar_target(sample_sizes, c(30L, 50L, 100L, 200L, 500L)),
  tar_target(reps, seq_len(6L)),

  tar_target(worlds, lapply(seq_len(nrow(world_grid)), function(i)
    make_world(world_id = i - 1L,
               nonstationary = world_grid$nonstationary[i],
               spatial_range_km = world_grid$spatial_range_km[i],
               n_time = 5L))),

  tar_target(results, {
    out <- list(); k <- 1
    for (wi in seq_along(worlds)) for (n in sample_sizes)
      for (d in designs) for (r in reps) {
        res <- run_scenario(worlds[[wi]], n = n, design = d)
        res$world_id <- wi - 1L; res$rep <- r
        out[[k]] <- res; k <- k + 1
      }
    do.call(rbind, out)
  }),

  tar_target(decomposition, {
    agg <- aggregate(cbind(algorithmic, spatial, temporal, total) ~ world_id,
                     data = results, FUN = mean)
    colMeans(agg[, c("algorithmic", "spatial", "temporal", "total")])
  })
)
