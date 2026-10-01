#' Area of applicability (Meyer and Pebesma 2021)
#'
#' Dissimilarity index of each query point: the minimum distance in importance
#' weighted standardised predictor space to the training reference, divided by
#' the mean nearest neighbour distance within the reference. Points whose index
#' exceeds the \code{quantile_thresh} quantile of the within reference indices
#' are outside the area of applicability. The reference should be one row per
#' site (or per site year with a between site calibration, as in Supplement S9)
#' so that repeated visits do not collapse the distance scale.
#'
#' @param train_pred Matrix or data frame of reference predictors.
#' @param query_pred Matrix or data frame of query predictors (same columns).
#' @param importance Optional non negative weights, one per predictor; the
#'   square root is applied, as in the article.
#' @param quantile_thresh Quantile defining the threshold (0.95).
#' @return A list with \code{DI}, \code{threshold} and logical \code{inside}.
#' @references Meyer, H. and Pebesma, E. (2021) Predicting into unknown space?
#'   Estimating the area of applicability of spatial prediction models.
#'   Methods in Ecology and Evolution 12, 1620 to 1633.
#' @export
area_of_applicability <- function(train_pred, query_pred, importance = NULL,
                                  quantile_thresh = 0.95) {
  Xtr <- as.matrix(train_pred); Xq <- as.matrix(query_pred)
  mu <- colMeans(Xtr); sdv <- apply(Xtr, 2, stats::sd) + 1e-12
  Ztr <- sweep(sweep(Xtr, 2, mu), 2, sdv, "/")
  Zq <- sweep(sweep(Xq, 2, mu), 2, sdv, "/")
  if (!is.null(importance)) {
    w <- sqrt(pmax(importance, 0))
    Ztr <- Ztr * rep(w, each = nrow(Ztr)); Zq <- Zq * rep(w, each = nrow(Zq))
  }
  Dtr <- dist_km(Ztr, Ztr); diag(Dtr) <- Inf
  nn_tr <- apply(Dtr, 1, min); dbar <- mean(nn_tr)
  nn_q <- apply(dist_km(Zq, Ztr), 1, min)
  DI <- nn_q / dbar
  thresh <- as.numeric(stats::quantile(nn_tr / dbar, quantile_thresh))
  list(DI = DI, threshold = thresh, inside = DI <= thresh)
}
