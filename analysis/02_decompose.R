# =============================================================================
# 02_decompose.R — Shapley decomposition tables T1 to T5 from ablation_expanded.csv
# =============================================================================
# Exact Shapley over the 16 coalitions of {COORD, EDF, LAG, TEMP}; 95% intervals
# from a cluster bootstrap over independent field realisations (B = 2000).
# Writes results/expanded_T1_decomposition_by_design.csv, T2_temporal_by_T,
# T3_detection, T4_abundance, T5_NN0_by_design, T6_temporal_by_T_n (temporal
# Shapley, its share of the total gain and the NN0 controlled increment within
# each timesteps x sample size cell) and T7_temporal_by_n_design (T >= 3).
# Usage: Rscript 02_decompose.R
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

  # T6: temporal value resolved by timesteps x sample size (pooled over designs and worlds)
  share_ci <- function(s, B = 2000, seed = 1) {   # ratio of means, cluster bootstrap over fields
    set.seed(seed); f <- split(s, s$field); G <- length(f)
    m <- sapply(f, function(x) c(mean(x$phi_TEMP), mean(x$total)))
    dr <- replicate(B, { k <- sample.int(G, G, TRUE); 100 * mean(m[1, k]) / mean(m[2, k]) })
    c(value = 100 * mean(s$phi_TEMP) / mean(s$total), ci_lo = unname(quantile(dr, 0.025)), ci_hi = unname(quantile(dr, 0.975)))
  }
  out$T6_temporal_by_T_n <- do.call(rbind, lapply(split(B, list(B$n_time, B$n)), function(s) {
    rbind(data.frame(n_time = s$n_time[1], n = s$n[1], quantity = "Shapley", t(boot_ci(s, "phi_TEMP")), n_datasets = nrow(s)),
          data.frame(n_time = s$n_time[1], n = s$n[1], quantity = "Added last, NN0 controlled", t(boot_ci(s, "temp_marg_withNN0")), n_datasets = nrow(s)),
          if (s$n_time[1] > 1) data.frame(n_time = s$n_time[1], n = s$n[1], quantity = "Share of total gain (%)", t(share_ci(s)), n_datasets = nrow(s))) }))
  out$T6_temporal_by_T_n <- out$T6_temporal_by_T_n[order(out$T6_temporal_by_T_n$quantity, out$T6_temporal_by_T_n$n, out$T6_temporal_by_T_n$n_time), ]
  # T7: temporal Shapley by sample size x design, pooled over T in {3, 5, 10}
  B3 <- B[B$n_time >= 3, ]
  out$T7_temporal_by_n_design <- do.call(rbind, lapply(split(B3, list(B3$n, B3$design)), function(s)
    data.frame(n = s$n[1], design = s$design[1], t(boot_ci(s, "phi_TEMP")), n_datasets = nrow(s))))

  for (nm in names(out)) {
    write.csv(out[[nm]], file.path(OUT_DIR, paste0("expanded_", nm, ".csv")), row.names = FALSE)
    cat("\n", strrep("=", 78), "\n", nm, "\n", strrep("=", 78), "\n", sep = "")
    print(out[[nm]], digits = 4, row.names = FALSE)
  }
  invisible(out)
}
if (sys.nframe() == 0) main()
