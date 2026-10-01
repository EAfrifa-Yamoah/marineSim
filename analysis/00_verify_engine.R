# =============================================================================
# 00_verify_engine.R — run before the full grid and before trusting any number.
# =============================================================================
# Prints PASS/FAIL for each mechanism check of RUN_EXPANDED_R.md. The same checks
# run as the package's testthat suite (marineSim/tests/testthat/test-engine.R);
# this script is the human readable form. Usage: Rscript 00_verify_engine.R
# =============================================================================
suppressPackageStartupMessages(library(marineSim))
ok <- function(cond, msg) cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))

cat("== 1. World mechanisms match Section 3.1 ==\n")
w <- make_world(stationary = FALSE, range_km = 80, n_time = 5, seed = 11)
ok(abs(sd(as.vector(w$depth)) - 1) < 0.02 && abs(mean(w$depth)) < 0.02, "covariates are standardised (mean 0, sd 1)")
ok(abs(mean(w$mu[[1]] == 0) - 0.25) < 0.02, sprintf("habitat mask gives ~25%% structural zeros (%.3f)", mean(w$mu[[1]] == 0)))
sst_change <- sd(as.vector(w$sst[[5]] - w$sst[[1]])); ok(sst_change > 0.5, "SST evolves through time (AR(1))")
# only SST varies in time: eta_5 - eta_1 must equal b_sst * (sst_5 - sst_1) exactly
resid <- (w$eta[[5]] - w$eta[[1]]) - w$b_sst * (w$sst[[5]] - w$sst[[1]])
ok(max(abs(resid)) < 1e-9, "persistent site effect and latent field are static (only SST changes eta over time)")
ok(abs(sd(as.vector(w$b_depth)) - 0.20) < 0.04, sprintf("non stationary coefficient perturbation sd ~0.20 (%.3f)", sd(as.vector(w$b_depth))))
wn <- make_world(stationary = TRUE, range_km = 80, n_time = 3, seed = 11)
ok(TRUE, "stationary and non stationary worlds both construct")

cat("\n== 2. Abundance response is ecological ==\n")
yc <- as.vector(w$ycount[[5]])
ok(mean(yc) > 3 && mean(yc) < 15 && max(yc) < 2000, sprintf("counts plausible: mean %.1f, median %.0f, p95 %.0f, max %d, zeros %.2f",
   mean(yc), median(yc), quantile(yc, 0.95), max(yc), mean(yc == 0)))

cat("\n== 3. Sampling designs ==\n")
set.seed(1); for (dsg in c("random", "clustered", "stratified")) {
  s <- sample_sites(100, dsg); ok(length(unique(s)) == 100 && all(s >= 1 & s <= GRID^2), sprintf("%s: 100 unique valid cells", dsg)) }
set.seed(2); sc <- cell_xy(sample_sites(100, "clustered"))
ok(sd(sc[, 1]) < sd(cell_xy(sample_sites(100, "random"))[, 1]), "clustered design is more spatially concentrated than random")

cat("\n== 4. Determinism (same seed -> identical dataset and scores) ==\n")
d1 <- build_dataset(w, 100, 5, "random", seed = 7); d2 <- build_dataset(w, 100, 5, "random", seed = 7)
ok(identical(d1$y_tr, d2$y_tr) && identical(d1$btr$LAG, d2$btr$LAG) && identical(d1$bte$TEMP, d2$bte$TEMP), "dataset build is deterministic")
ok(abs(score(d1, BLOCKS, 7) - score(d2, BLOCKS, 7)) < 1e-12, "forest score is deterministic given seed")

cat("\n== 5. Leakage: test features depend on training responses only ==\n")
# Test responses are never an input to any block builder. Non vacuity: perturbing a
# TRAINING response must move the test lag features.
lag_a <- spatial_lags(d1$btr$COORD, d1$y_tr, d1$bte$COORD, FALSE)
y_pert <- d1$y_tr; y_pert[1] <- 1 - y_pert[1]
lag_b <- spatial_lags(d1$btr$COORD, y_pert, d1$bte$COORD, FALSE)
ok(any(lag_a != lag_b), "non vacuity: perturbing a training response moves test lag features")
ok(!("y_te" %in% names(formals(spatial_lags))) && !("y_te" %in% names(formals(nearest_train))), "block builders take no held out response argument")

cat("\n== 6. Placebo: at T = 1 the temporal block is constant ==\n")
d0 <- build_dataset(w, 100, 1, "random", seed = 9)
ok(all(d0$btr$TEMP == d0$btr$TEMP[1, 1]) && all(d0$bte$TEMP == d0$btr$TEMP[1, 1]), "TEMP features constant at T = 1 (carry no information)")

cat("\n== 7. Sanity of scores ==\n")
g <- glm_score(d1); base <- score(d1, character(0), 7); full <- score(d1, BLOCKS, 7)
ok(g > 0.5 && base > 0.5 && full > 0.5 && full <= 1, sprintf("AUCs in range: GLM %.3f  RF %.3f  stRF %.3f", g, base, full))
ok(full > base, "full coalition beats covariate only forest on this dataset")
da <- build_dataset(w, 100, 5, "stratified", seed = 7, response = "abundance")
ok(!is.na(score(da, BLOCKS, 7)) && !is.na(glm_score(da)), sprintf("abundance path scores: GLM rho %.3f  stRF rho %.3f", glm_score(da), score(da, BLOCKS, 7)))
d_full <- build_dataset(w, 100, 5, "clustered", seed = 7, detection = 1.0)
d_thin <- build_dataset(w, 100, 5, "clustered", seed = 7, detection = 0.7)
ok(all(d_thin$y_tr <= d_full$y_tr) && sum(d_thin$y_tr) < sum(d_full$y_tr) && identical(d_thin$btr$COORD, d_full$btr$COORD),
   sprintf("detection 0.7 only removes positives at the same training sites (%d -> %d positives)", sum(d_full$y_tr), sum(d_thin$y_tr)))

cat("\n== 8. Reference pattern from the Python engine (qualitative, random design, T in {1,3,5}) ==\n")
cat("   Expected once the R grid runs: distance fields the largest block (~35% share);\n")
cat("   spatial total ~65-70%; algorithmic ~13-17%; temporal Shapley ~0 at T=1, rising with T;\n")
cat("   clustered design: total gain lower and spatial lag share roughly halved.\n")
cat("   Numbers will differ (different random stream); orderings and signs should agree.\n")
