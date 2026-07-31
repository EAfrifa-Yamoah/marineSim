test_that("matern_grf returns a standardised field of the requested size", {
  g <- matern_grf(grid_n = 32, range_km = 60, seed = 1)
  expect_equal(dim(g$field), c(32, 32))
  expect_lt(abs(mean(g$field)), 1e-6)
  expect_lt(abs(stats::sd(as.vector(g$field)) - 1), 0.05)
})

test_that("evaluate_auc matches a known ordering", {
  y <- c(0, 0, 1, 1)
  p <- c(0.1, 0.2, 0.3, 0.9)
  expect_equal(evaluate_auc(y, p), 1)
  expect_true(is.na(evaluate_auc(rep(1, 4), p)))
})

test_that("decompose_gain sums to the total", {
  auc <- c(GLM = 0.60, RF = 0.62, SpatialRF = 0.73, stRF = 0.77)
  d <- decompose_gain(auc)
  expect_equal(unname(d["algorithmic"] + d["spatial"] + d["temporal"]),
               unname(d["total"]), tolerance = 1e-10)
})

test_that("strf_features returns one row per query point with named columns", {
  set.seed(3)
  df <- data.frame(x = runif(40, 0, 400), y = runif(40, 0, 400),
                   time = rep(1:2, 20), depth = rnorm(40), subst = rnorm(40),
                   wave = rnorm(40), sst = rnorm(40), y = rbinom(40, 1, 0.5))
  X <- strf_features(df, df)
  expect_equal(nrow(X), nrow(df))
  expect_true(all(c("lag25", "tlag", "amean", "x", "y") %in% colnames(X)))
})

test_that("area_of_applicability flags clear extrapolation", {
  set.seed(5)
  tr <- matrix(rnorm(200), 100, 2)
  q <- rbind(matrix(rnorm(20), 10, 2), matrix(rnorm(10, mean = 8), 5, 2))
  a <- area_of_applicability(tr, q)
  expect_true(all(a$inside[1:10] | !a$inside[1:10]))  # defined
  expect_true(mean(a$inside[11:15]) < mean(a$inside[1:10]))
})
