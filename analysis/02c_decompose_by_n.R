# =============================================================================
# 02c_decompose_by_n.R — Shapley decomposition along the sample size gradient
# =============================================================================
# Exact Shapley over the 16 coalitions of {COORD, EDF, LAG, TEMP} within each
# per timestep sample size n (main grid: binary response, detection 1, spans
# 1, 3 and 5), pooled over designs (T9) and within each design (T10). 95%
# intervals from a cluster bootstrap over field realisations (B = 2000), as in
# 02_decompose.R. Writes results/expanded_T9_decomposition_by_n.csv and
# results/expanded_T10_decomposition_by_n_design.csv.
# Usage: Rscript 02c_decompose_by_n.R
# =============================================================================
suppressPackageStartupMessages(library(marineSim))
OUT_DIR <- Sys.getenv("ABL_OUT_DIR", "../results")

main <- function() {
  d <- read.csv(file.path(OUT_DIR, "ablation_expanded.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  d <- add_shapley(d)
  A <- d[d$response == "binary" & d$detection == 1 & d$n_time %in% c(1, 3, 5), ]
  ns <- sort(unique(A$n))
  t9 <- do.call(rbind, lapply(ns, function(n) {
    r <- decomp_table(A[A$n == n, ], paste0("n = ", n)); r$n <- n; r }))
  t10 <- do.call(rbind, lapply(c("random", "stratified", "clustered"), function(dsg)
    do.call(rbind, lapply(ns, function(n) {
      r <- decomp_table(A[A$n == n & A$design == dsg, ], paste0(dsg, ", n = ", n)); r$design <- dsg; r$n <- n; r }))))
  write.csv(t9, file.path(OUT_DIR, "expanded_T9_decomposition_by_n.csv"), row.names = FALSE)
  write.csv(t10, file.path(OUT_DIR, "expanded_T10_decomposition_by_n_design.csv"), row.names = FALSE)
  cat(sprintf("T9: %d rows; T10: %d rows\n", nrow(t9), nrow(t10)))
  print(t9[t9$component %in% c("Algorithmic flexibility", "Spatial total", "Temporal", "Total gain over GLM"),
           c("subset", "component", "value", "ci_lo", "ci_hi", "share_pct")], row.names = FALSE, digits = 3)
}
main()
