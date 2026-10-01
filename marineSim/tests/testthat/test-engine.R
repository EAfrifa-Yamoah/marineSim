# The verification protocol of RUN_EXPANDED_R.md, as unit tests. Every mechanism
# below is a stated property of Section 3.1 to 3.3 of the article.

w <- make_world(stationary = FALSE, range_km = 80, n_time = 5, seed = 11)

test_that("world mechanisms match Section 3.1", {
  expect_lt(abs(sd(as.vector(w$depth)) - 1), 0.02)
  expect_lt(abs(mean(w$depth)), 0.02)
  expect_lt(abs(mean(w$mu[[1]] == 0) - 0.25), 0.02)          # habitat mask
  expect_gt(sd(as.vector(w$sst[[5]] - w$sst[[1]])), 0.5)     # SST evolves
  resid <- (w$eta[[5]] - w$eta[[1]]) - w$b_sst * (w$sst[[5]] - w$sst[[1]])
  expect_lt(max(abs(resid)), 1e-9)                            # only SST changes eta
  expect_lt(abs(sd(as.vector(w$b_depth)) - 0.20), 0.04)      # non stationary sd
  ws <- make_world(stationary = TRUE, range_km = 80, n_time = 3, seed = 11)
  expect_true(all(ws$b_depth == -0.55))
})

test_that("abundance counts are ecological", {
  yc <- as.vector(w$ycount[[5]])
  expect_gt(mean(yc), 3); expect_lt(mean(yc), 15); expect_lt(max(yc), 2000)
})

test_that("sampling designs draw valid distinct cells", {
  set.seed(1)
  for (dsg in c("random", "clustered", "stratified")) {
    s <- sample_sites(100, dsg)
    expect_equal(length(unique(s)), 100)
    expect_true(all(s >= 1 & s <= GRID^2))
  }
  set.seed(2)
  expect_lt(sd(cell_xy(sample_sites(100, "clustered"))[, 1]),
            sd(cell_xy(sample_sites(100, "random"))[, 1]))
})

test_that("dataset build and forest score are deterministic", {
  d1 <- build_dataset(w, 100, 5, "random", seed = 7)
  d2 <- build_dataset(w, 100, 5, "random", seed = 7)
  expect_identical(d1$y_tr, d2$y_tr)
  expect_identical(d1$btr$LAG, d2$btr$LAG)
  expect_identical(d1$bte$TEMP, d2$bte$TEMP)
  expect_lt(abs(score(d1, BLOCKS, 7) - score(d2, BLOCKS, 7)), 1e-12)
})

test_that("block builders never see a held out response", {
  d1 <- build_dataset(w, 100, 5, "random", seed = 7)
  lag_a <- spatial_lags(d1$btr$COORD, d1$y_tr, d1$bte$COORD, FALSE)
  y_pert <- d1$y_tr; y_pert[1] <- 1 - y_pert[1]
  lag_b <- spatial_lags(d1$btr$COORD, y_pert, d1$bte$COORD, FALSE)
  expect_true(any(lag_a != lag_b))                            # non vacuity
  expect_false("y_te" %in% names(formals(spatial_lags)))
  expect_false("y_te" %in% names(formals(nearest_train)))
})

test_that("temporal block is constant at T = 1 (placebo)", {
  d0 <- build_dataset(w, 100, 1, "random", seed = 9)
  expect_true(all(d0$btr$TEMP == d0$btr$TEMP[1, 1]))
  expect_true(all(d0$bte$TEMP == d0$btr$TEMP[1, 1]))
})

test_that("scores are in range and the full coalition beats the covariate only forest", {
  d1 <- build_dataset(w, 100, 5, "random", seed = 7)
  g <- glm_score(d1); base <- score(d1, character(0), 7); full <- score(d1, BLOCKS, 7)
  expect_true(g > 0.5 && base > 0.5 && full > 0.5 && full <= 1)
  expect_gt(full, base)
  da <- build_dataset(w, 100, 5, "stratified", seed = 7, response = "abundance")
  expect_false(is.na(score(da, BLOCKS, 7)))
  expect_false(is.na(glm_score(da)))
})

test_that("detection thinning removes positives at the same training sites only", {
  d_full <- build_dataset(w, 100, 5, "clustered", seed = 7, detection = 1.0)
  d_thin <- build_dataset(w, 100, 5, "clustered", seed = 7, detection = 0.7)
  expect_true(all(d_thin$y_tr <= d_full$y_tr))
  expect_lt(sum(d_thin$y_tr), sum(d_full$y_tr))
  expect_identical(d_thin$btr$COORD, d_full$btr$COORD)
})

test_that("Shapley values sum to the full coalition gain over the covariate only forest", {
  set.seed(3)
  row <- list(); for (S in required_subsets()[1:16]) row[[skey(S)]] <- runif(1)
  phi <- shapley_row(row)
  expect_equal(unname(sum(phi)), row[[skey(BLOCKS)]] - row[[skey(character(0))]], tolerance = 1e-12)
  expect_equal(length(required_subsets()), 20)
})

test_that("area of applicability flags clear extrapolation", {
  set.seed(5)
  tr <- matrix(rnorm(200), 100, 2)
  q <- rbind(matrix(rnorm(20), 10, 2), matrix(rnorm(10, mean = 8), 5, 2))
  a <- area_of_applicability(tr, q)
  expect_true(all(!a$inside[11:15]))
  expect_gt(mean(a$inside[1:10]), 0.5)
})
