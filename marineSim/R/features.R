#' Euclidean distance matrix between two point sets
#' @param A,B Numeric matrices with coordinates in columns (km).
#' @return A \code{nrow(A)} by \code{nrow(B)} matrix of distances.
#' @export
dist_km <- function(A, B) {                    # rows of A vs rows of B
  outer(rowSums(A^2), rowSums(B^2), "+") - 2 * tcrossprod(A, B) -> d2
  sqrt(pmax(d2, 0))
}

#' Multi scale spatial lag block (L)
#'
#' For each query point and each radius in \code{RADII_KM}, the mean training
#' response within that radius, excluding the point itself when the query set
#' is the training set (leave one out), and falling back to the training mean
#' when no neighbour lies within the radius. Only training responses enter;
#' held out responses are never an argument.
#'
#' @param train_xy,train_y Training coordinates (km) and responses.
#' @param query_xy Query coordinates (km).
#' @param is_train \code{TRUE} when \code{query_xy} is \code{train_xy}, so that
#'   self neighbours are excluded.
#' @return A matrix with one column per radius.
#' @export
spatial_lags <- function(train_xy, train_y, query_xy, is_train) {
  D <- dist_km(query_xy, train_xy); gm <- mean(train_y)
  out <- matrix(NA_real_, nrow(query_xy), length(RADII_KM))
  for (k in seq_along(RADII_KM)) {
    mask <- D <= RADII_KM[k]
    if (is_train) diag(mask) <- FALSE
    cnt <- rowSums(mask); s <- as.vector(mask %*% train_y)
    out[, k] <- ifelse(cnt > 0, s / pmax(cnt, 1), gm)
  }
  out
}

#' Index of the nearest training site for each query point
#' @param train_xy Training coordinates (km).
#' @param query_xy Query coordinates (km).
#' @param exclude_self Exclude the zero distance match (query set is the training set).
#' @return Integer vector of training row indices.
#' @export
nearest_train <- function(train_xy, query_xy, exclude_self = FALSE) {
  D <- dist_km(query_xy, train_xy)
  if (exclude_self) diag(D) <- Inf
  apply(D, 1, which.min)
}

#' Euclidean distance field block (E)
#' @param query_xy Query coordinates (km).
#' @param anchors Anchor coordinates, usually k-means centres of the training sites.
#' @return A matrix with one column per anchor.
#' @export
edf <- function(query_xy, anchors) dist_km(query_xy, anchors)
