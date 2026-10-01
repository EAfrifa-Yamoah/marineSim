# =============================================================================
# 06_s6_sdmtmb.R — Supplement S6: stRF versus sdmTMB (spatial and spatiotemporal
# AR(1) Matérn random field GLMMs), realisation 0 of the six worlds, random design.
# Writes results/si/tableS6_sdmtmb.csv. Requires sdmTMB (>= 0.6) and INLA free
# meshes via fmesher. Usage: Rscript 06_s6_sdmtmb.R (~2 h; timings depend on the
# TMB build flags, see README)
# =============================================================================
suppressPackageStartupMessages({library(marineSim); library(sdmTMB)})
OUT_DIR <- Sys.getenv("ABL_OUT_DIR", "../results")
OUT_S6 <- file.path(OUT_DIR, "si"); dir.create(OUT_S6, showWarnings = FALSE, recursive = TRUE)
CUTOFF <- 20
rows <- list()
for (wg in world_grid(0)) {
  w <- make_world(wg$stat, wg$range_km, n_time = 5, seed = wg$seed)
  for (n in c(30, 50, 100, 200, 500)) for (rep in 0:2) {
    ds <- build_dataset(w, n, 5, "random", seed = wg$seed * 17 + n * 3 + rep)
    if (is.null(ds)) next
    xy <- ds$btr$COORD; nd <- data.frame(ds$Xc_te, X = ds$bte$COORD[, 1], Y = ds$bte$COORD[, 2])
    # stRF and GLM
    t0 <- Sys.time(); a_strf <- score(ds, BLOCKS, 7); t_strf <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    a_glm <- glm_score(ds)
    # sdmTMB spatial: same final timestep rows
    d1 <- data.frame(y = ds$y_tr, ds$Xc_tr, X = xy[, 1], Y = xy[, 2])
    mesh <- make_mesh(d1, c("X", "Y"), cutoff = CUTOFF)
    t0 <- Sys.time()
    f1 <- tryCatch(sdmTMB(y ~ depth + I(depth^2) + subst + wave + sst, data = d1, mesh = mesh,
                          family = binomial(), spatial = "on"), error = function(e) NULL)
    t_sp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ok1 <- !is.null(f1) && isTRUE(sanity(f1, silent = TRUE)$all_ok)
    a_sp <- if (!is.null(f1)) auc(ds$y_te, plogis(predict(f1, newdata = nd)$est)) else NA
    # sdmTMB spatiotemporal AR(1): all five timesteps as rows, predict final timestep
    d5 <- do.call(rbind, lapply(1:5, function(t) data.frame(y = ds$obs[[t]], covariates(w, ds$site, t),
                                                              X = xy[, 1], Y = xy[, 2], t = t)))
    mesh5 <- make_mesh(d5, c("X", "Y"), cutoff = CUTOFF)
    t0 <- Sys.time()
    f5 <- tryCatch(sdmTMB(y ~ depth + I(depth^2) + subst + wave + sst, data = d5, mesh = mesh5, time = "t",
                          family = binomial(), spatial = "on", spatiotemporal = "ar1"), error = function(e) NULL)
    t_st <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ok5 <- !is.null(f5) && isTRUE(sanity(f5, silent = TRUE)$all_ok)
    a_st <- if (!is.null(f5)) auc(ds$y_te, plogis(predict(f5, newdata = cbind(nd, t = 5))$est)) else NA
    rows[[length(rows) + 1]] <- data.frame(cfg = wg$cfg, n = n, rep = rep, mesh_vertices = mesh$mesh$n,
      GLM = a_glm, stRF = a_strf, sdmTMB_spatial = a_sp, sdmTMB_spatiotemporal = a_st,
      conv_spatial = ok1, conv_st = ok5, t_stRF = t_strf, t_spatial = t_sp, t_st = t_st)
    write.csv(do.call(rbind, rows), file.path(OUT_S6, "tableS6_sdmtmb.csv"), row.names = FALSE)
  }
  cat("S6:", wg$cfg, "done |", length(rows), "datasets |", format(Sys.time(), "%H:%M:%S"), "\n")
}
cat("S6 COMPLETE\n")
