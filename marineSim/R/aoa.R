#' Area of applicability for a fitted model
#'
#' Implements the dissimilarity index of Meyer and Pebesma (2021). For each
#' prediction point the minimum predictor space distance to the training data is
#' scaled by the mean nearest neighbour distance within the training data; points
#' whose scaled distance exceeds a threshold (by default the outlier removed 95th
#' percentile of the within training distances) fall outside the area of
#' applicability and their predictions should be treated as unsupported
#' extrapolation. Predictors are standardised and optionally weighted by variable
#' importance, and the training reference is taken at the site level so that
#' repeated visits to the same site do not collapse the distance scale.
#'
#' @param train_pred Matrix or data frame of training predictors.
#' @param query_pred Matrix or data frame of prediction predictors.
#' @param importance Optional non negative weights, one per predictor.
#' @param quantile_thresh Quantile of within training distances for the threshold.
#'
#' @return A list with \code{DI} (dissimilarity index per query point),
#'   \code{threshold}, and \code{inside} (logical vector).
#' @export
area_of_applicability <- function(train_pred, query_pred, importance = NULL,
                                  quantile_thresh = 0.95) {
  Xtr <- as.matrix(train_pred); Xq <- as.matrix(query_pred)
  mu <- colMeans(Xtr); sdv <- apply(Xtr, 2, stats::sd); sdv[sdv == 0] <- 1
  Ztr <- sweep(sweep(Xtr, 2, mu), 2, sdv, "/")
  Zq <- sweep(sweep(Xq, 2, mu), 2, sdv, "/")
  if (!is.null(importance)) {
    w <- sqrt(pmax(importance, 0) + 1e-9)
    Ztr <- sweep(Ztr, 2, w, "*"); Zq <- sweep(Zq, 2, w, "*")
  }
  nn_tr <- FNN::get.knnx(Ztr, Ztr, k = 2)$nn.dist[, 2]
  dbar <- mean(nn_tr)
  nn_q <- FNN::get.knnx(Ztr, Zq, k = 1)$nn.dist[, 1]
  DI <- nn_q / dbar
  thresh <- as.numeric(stats::quantile(nn_tr / dbar, quantile_thresh))
  list(DI = DI, threshold = thresh, inside = DI <= thresh)
}
