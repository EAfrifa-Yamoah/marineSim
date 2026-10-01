# =============================================================================
# 01_ablation_expanded.R — feature block ablation over the expanded grid
# =============================================================================
# Fits, for every dataset, the GLM and the 20 coalitions of required_subsets()
# (16 for the exact Shapley decomposition plus four NN0 controls) and appends one
# row per dataset to results/ablation_expanded.csv. Resumable: completed keys are
# skipped, so an interrupted run is restarted with the same command.
#
# Usage (from analysis/):
#   Rscript 01_ablation_expanded.R                 # all four stages, ~4-5 h on one core
#   Rscript 01_ablation_expanded.R A_designs       # one stage
#   ABL_SHARD=0/4 Rscript 01_ablation_expanded.R   # one of four shards (see README)
#   ABL_BUDGET=120 Rscript 01_ablation_expanded.R  # stop after 120 s (smoke test)
#
# Environment: ABL_OUT_DIR (default ../results), ABL_SHARD (k/K), ABL_BUDGET (s).
# =============================================================================
suppressPackageStartupMessages(library(marineSim))

OUT_DIR <- Sys.getenv("ABL_OUT_DIR", "../results")
SHARD   <- Sys.getenv("ABL_SHARD", "")
BUDGET  <- as.numeric(Sys.getenv("ABL_BUDGET", "0"))
sfx <- if (nzchar(SHARD)) paste0("_shard", strsplit(SHARD, "/")[[1]][1]) else ""
OUT <- file.path(OUT_DIR, paste0("ablation_expanded", sfx, ".csv"))
KEY <- c("response", "detection", "design", "world_cfg", "realisation", "n", "n_time", "rep")

T_START <- Sys.time()
run_stage <- function(stage, response, detection, designs, realisations, timesteps,
                      sizes = c(30, 50, 100, 200, 500), reps = 1) {
  subs <- required_subsets()
  done <- character(0)
  if (file.exists(OUT)) { prev <- read.csv(OUT, stringsAsFactors = FALSE, check.names = FALSE)
    done <- do.call(paste, c(prev[KEY], sep = "|")) }
  grid <- world_grid(realisations)
  if (nzchar(SHARD)) { kn <- as.integer(strsplit(SHARD, "/")[[1]])
    grid <- grid[(seq_along(grid) - 1) %% kn[2] == kn[1]] }
  nrow_done <- 0L
  for (wg in grid) {
    w <- make_world(wg$stat, wg$range_km, n_time = max(timesteps), seed = wg$seed)
    for (di in seq_along(designs)) for (n in sizes) for (nt in timesteps) for (rep in 0:(reps - 1)) {
      key <- paste(response, detection, designs[di], wg$cfg, wg$real, n, nt, rep, sep = "|")
      if (key %in% done) next
      if (BUDGET > 0 && as.numeric(difftime(Sys.time(), T_START, units = "secs")) > BUDGET) {
        cat(sprintf("[%s] budget reached; rows %d (resumable)\n", stage, nrow_done)); return(invisible()) }
      dseed <- dataset_seed(wg$seed, n, nt, di, rep, detection, response)
      ds <- build_dataset(w, n, nt, designs[di], dseed, detection, response)
      if (is.null(ds)) next
      rec <- data.frame(stage = stage, response = response, detection = detection,
                        design = designs[di], world_cfg = wg$cfg, world_id = wg$wid,
                        realisation = wg$real, stationary = wg$stat, range_km = wg$range_km,
                        n = n, n_time = nt, rep = rep, GLM = glm_score(ds), stringsAsFactors = FALSE)
      for (S in subs) rec[[skey(S)]] <- score(ds, S, seed = dseed %% 10000L)
      write.table(rec, OUT, sep = ",", row.names = FALSE, col.names = !file.exists(OUT),
                  append = file.exists(OUT))
      nrow_done <- nrow_done + 1L
    }
    cat(sprintf("[%s] world %d %s r%d done | rows %d | %.0fs\n", stage, wg$wid, wg$cfg,
                wg$real, nrow_done, as.numeric(difftime(Sys.time(), T_START, units = "secs"))))
  }
  cat(sprintf("[%s] COMPLETE rows=%d\n", stage, nrow_done))
}

STAGES <- list(
  A_designs   = list(response = "binary",    detection = 1.0, designs = c("random", "clustered", "stratified"),
                     realisations = 0:5, timesteps = c(1, 3, 5)),
  B_T10       = list(response = "binary",    detection = 1.0, designs = c("random", "clustered", "stratified"),
                     realisations = 0:5, timesteps = 10),
  C_detect    = list(response = "binary",    detection = 0.7, designs = c("random", "clustered", "stratified"),
                     realisations = 0:2, timesteps = c(1, 3, 5)),
  D_abundance = list(response = "abundance", detection = 1.0, designs = c("random", "clustered", "stratified"),
                     realisations = 0:2, timesteps = c(1, 3, 5)))

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE); if (!length(args)) args <- names(STAGES)
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  for (s in args) do.call(run_stage, c(list(stage = s), STAGES[[s]]))
}
