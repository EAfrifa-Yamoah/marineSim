# =============================================================================
# 02_decompose.R — Shapley decomposition tables T1 to T5 from ablation_expanded.csv
# =============================================================================
# Exact Shapley over the 16 coalitions of {COORD, EDF, LAG, TEMP}; 95% intervals
# from a cluster bootstrap over independent field realisations (B = 2000).
# Writes results/expanded_T1_decomposition_by_design.csv, T2_temporal_by_T,
# T3_detection, T4_abundance, T5_NN0_by_design. Usage: Rscript 02_decompose.R
# =============================================================================
suppressPackageStartupMessages(library(marineSim))
OUT_DIR <- Sys.getenv("ABL_OUT_DIR", "../results")

main <- function() {
  d <- read.csv(file.path(OUT_DIR, "ablation_expanded.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  d <- add_shapley(d)
  cat(sprintf("loaded %d datasets; stages: %s\n", nrow(d), paste(unique(d$stage), collapse = ", ")))
  out <- list()
  A <- d[d$response == "binary" & d$detection == 1 & d$n_time %in% c(1, 3, 5), ]
  t1 <- list(decomp_table(A, "pooled (all designs)"))
  for (dsg in c("random", "clustered", "stratified"))
    if (any(A$design == dsg)) t1[[length(t1) + 1]] <- decomp_table(A[A$design == dsg, ], dsg)
  out$T1_decomposition_by_design <- do.call(rbind, t1)

  B <- d[d$response == "binary" & d$detection == 1, ]
  out$T2_temporal_by_T <- do.call(rbind, lapply(sort(unique(B$n_time)), function(T) {
    s <- B[B$n_time == T, ]
    do.call(rbind, lapply(c(Shapley = "phi_TEMP", `Added last` = "temp_added_last",
                            `Added last, NN0 controlled` = "temp_marg_withNN0"), function(col) {
      ci <- boot_ci(s, col); data.frame(n_time = T, value = ci[["value"]], ci_lo = ci[["ci_lo"]],
                                       ci_hi = ci[["ci_hi"]], n_datasets = nrow(s)) })) -> r
    r$quantity <- rep(c("Shapley", "Added last", "Added last, NN0 controlled"), each = 1); r }))

  Cc <- d[d$response == "binary" & d$n_time %in% c(1, 3, 5) & d$realisation <= 2, ]
  out$T3_detection <- do.call(rbind, lapply(sort(unique(Cc$detection)), function(det)
    decomp_table(Cc[Cc$detection == det, ], paste("detection", det))))

  D <- d[d$response == "abundance", ]
  if (nrow(D)) {
    t4 <- list(decomp_table(D, "abundance pooled"))
    for (dsg in c("random", "clustered", "stratified"))
      if (any(D$design == dsg)) t4[[length(t4) + 1]] <- decomp_table(D[D$design == dsg, ], paste("abundance", dsg))
    out$T4_abundance <- do.call(rbind, t4)
  }

  out$T5_NN0_by_design <- do.call(rbind, lapply(c("random", "clustered", "stratified"), function(dsg) {
    s <- A[A$design == dsg & A$n_time >= 3, ]; if (!nrow(s)) return(NULL)
    rbind(data.frame(design = dsg, quantity = "Temporal marginal, no control", t(boot_ci(s, "temp_added_last")), n_datasets = nrow(s)),
          data.frame(design = dsg, quantity = "Temporal marginal, NN0 controlled", t(boot_ci(s, "temp_marg_withNN0")), n_datasets = nrow(s)))
  }))

  for (nm in names(out)) {
    write.csv(out[[nm]], file.path(OUT_DIR, paste0("expanded_", nm, ".csv")), row.names = FALSE)
    cat("\n", strrep("=", 78), "\n", nm, "\n", strrep("=", 78), "\n", sep = "")
    print(out[[nm]], digits = 4, row.names = FALSE)
  }
  invisible(out)
}
if (sys.nframe() == 0) main()
