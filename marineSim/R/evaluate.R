#' Area under the ROC curve
#'
#' A dependency free AUC via the Mann-Whitney U statistic.
#'
#' @param y Binary outcome (0/1).
#' @param p Predicted probabilities.
#' @return Scalar AUC, or NA if only one class is present.
#' @export
evaluate_auc <- function(y, p) {
  pos <- which(y == 1); neg <- which(y == 0)
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  r <- rank(p)
  (sum(r[pos]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg))
}

#' Decompose the predictive gain along the model ladder
#'
#' Given the AUC of the four nested models on a common test set, splits the total
#' gain over the GLM baseline into an algorithmic component (RF minus GLM), a
#' spatial component (SpatialRF minus RF) and a temporal component (stRF minus
#' SpatialRF).
#'
#' @param auc Named numeric vector with names GLM, RF, SpatialRF, stRF.
#' @return Named numeric vector: algorithmic, spatial, temporal, total.
#' @export
decompose_gain <- function(auc) {
  c(algorithmic = unname(auc["RF"] - auc["GLM"]),
    spatial = unname(auc["SpatialRF"] - auc["RF"]),
    temporal = unname(auc["stRF"] - auc["SpatialRF"]),
    total = unname(auc["stRF"] - auc["GLM"]))
}

#' Run a single benchmark scenario end to end
#'
#' Draws a training and an independent test sample from a world, fits all models
#' and returns their test AUC together with the gain decomposition.
#'
#' @param world A world object from \code{make_world}.
#' @param n Training sample size (number of sites).
#' @param design Sampling design for the training sample.
#' @param detection Detection probability.
#' @param test_n Size of the independent random test sample.
#' @param methods Character vector of methods to fit.
#' @return A one row data frame of AUC values plus decomposition columns.
#' @export
run_scenario <- function(world, n = 100L, design = "random", detection = 1.0,
                         test_n = 600L,
                         methods = c("GLM", "GAM", "RF", "SpatialRF", "stRF")) {
  train <- world$sample_sites(n, design, detection)
  test <- world$sample_sites(test_n, "random", 1.0)
  auc <- sapply(methods, function(m) {
    p <- tryCatch(fit_model(m, train, test), error = function(e) rep(NA, nrow(test)))
    evaluate_auc(test$y, p)
  })
  dec <- decompose_gain(auc[c("GLM", "RF", "SpatialRF", "stRF")])
  out <- as.data.frame(as.list(c(auc, dec)))
  out$n <- n; out$design <- design; out$detection <- detection
  out
}
