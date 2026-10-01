# =============================================================================
# 04_core_benchmark_table.R — the five method core benchmark (Figure 3, Table 2)
# =============================================================================
# Extracts GLM, standard RF (v_BASE), spatial RF (v_COORD) and stRF
# (v_COORD+EDF+LAG+TEMP) from Stage A of the ablation and merges the GAM scores.
# Writes results/core_five_methods.csv and prints the Table 2 empirical columns:
# the smallest n at which median test AUC reaches 0.75, and the percentage of the
# 45 design x size x span scenarios in which each method had the highest mean AUC.
# Usage: Rscript 04_core_benchmark_table.R
# =============================================================================
OUT_DIR <- Sys.getenv("ABL_OUT_DIR", "../results")
d <- read.csv(file.path(OUT_DIR, "ablation_expanded.csv"), check.names = FALSE, stringsAsFactors = FALSE)
A <- d[d$stage == "A_designs", ]
g <- read.csv(file.path(OUT_DIR, "si", "gam_core.csv"), stringsAsFactors = FALSE)
core <- data.frame(world_cfg = A$world_cfg, realisation = A$realisation, design = A$design, n = A$n, n_time = A$n_time,
                   GLM = A$GLM, RF = A[["v_BASE"]], SRF = A[["v_COORD"]], STRF = A[["v_COORD+EDF+LAG+TEMP"]])
core <- merge(core, g, by = c("world_cfg", "realisation", "design", "n", "n_time"))
core <- core[order(core$world_cfg, core$realisation, core$design, core$n, core$n_time), ]
M <- c("GLM", "GAM", "RF", "SRF", "STRF")
write.csv(core, file.path(OUT_DIR, "core_five_methods.csv"), row.names = FALSE)
cat("merged datasets:", nrow(core), "\n")
cat("mean AUC:\n"); print(round(colMeans(core[M]), 3))
med <- aggregate(core[M], by = list(n = core$n), FUN = median)
sc <- aggregate(core[M], by = list(design = core$design, n = core$n, n_time = core$n_time), FUN = mean)
best <- table(factor(M[apply(sc[M], 1, which.max)], levels = M))
cat("\nTable 2 columns\n")
for (m in M) {
  reach <- med$n[med[[m]] >= 0.75]
  cat(sprintf("  %-5s n for AUC >= 0.75: %-12s best in %.1f%% of %d scenarios\n", m,
              if (length(reach)) as.character(min(reach)) else "Not reached", 100 * best[[m]] / nrow(sc), nrow(sc)))
}
