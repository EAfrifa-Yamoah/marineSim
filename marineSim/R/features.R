#' Build spatiotemporal random forest features
#'
#' Adds the three feature families that define the stRF representation on top of
#' the raw environmental covariates and coordinates: multi scale spatial lags
#' (distance weighted neighbour means of the response at several radii), temporal
#' lags (the previous time point response at the same location, plus a basin wide
#' annual mean), and a Euclidean distance field (distance to a small set of
#' spatial anchors). All features are computed for a query set from a training
#' set so the same routine serves both in sample fitting and out of sample
#' prediction, including spatial or temporal extrapolation.
#'
#' @param train Data frame with columns x, y, time, y (response) and covariates.
#' @param query Data frame to build features for (defaults to \code{train}).
#' @param radii Numeric vector of spatial lag radii in the coordinate units.
#' @param n_anchor Number of distance field anchors.
#' @param covars Character vector of covariate column names to retain.
#'
#' @return A numeric matrix of features aligned to the rows of \code{query}.
#' @export
strf_features <- function(train, query = train,
                          radii = c(25, 75, 150), n_anchor = 8L,
                          covars = c("depth", "subst", "wave", "sst")) {
  stopifnot(all(c("x", "y", "time", "y") %in% names(train)))
  Xc <- as.matrix(query[, covars, drop = FALSE])

  # ---- multi scale spatial lags (per radius, matched within time) ----
  lag <- matrix(0, nrow(query), length(radii))
  for (i in seq_len(nrow(query))) {
    d <- sqrt((train$x - query$x[i])^2 + (train$y - query$y[i])^2)
    same_t <- train$time == query$time[i]
    for (j in seq_along(radii)) {
      sel <- same_t & d <= radii[j] & d > 1e-6
      if (sum(sel) < 2) sel <- d <= radii[j] & d > 1e-6
      if (sum(sel) >= 2) {
        w <- 1 / (d[sel] + 1)
        lag[i, j] <- sum(w * train$y[sel]) / sum(w)
      } else {
        k <- min(5, sum(d > 1e-6))
        nn <- order(d)[seq_len(k)]
        lag[i, j] <- mean(train$y[nn])
      }
    }
  }

  # ---- temporal lag: previous time point at the nearest training site ----
  tlag <- numeric(nrow(query))
  for (i in seq_len(nrow(query))) {
    prev_t <- train$time == (query$time[i] - 1)
    if (any(prev_t)) {
      d <- sqrt((train$x[prev_t] - query$x[i])^2 + (train$y[prev_t] - query$y[i])^2)
      tlag[i] <- train$y[prev_t][which.min(d)]
    } else {
      tlag[i] <- mean(train$y)
    }
  }
  annual <- tapply(train$y, train$time, mean)
  amean <- as.numeric(annual[as.character(query$time)])
  amean[is.na(amean)] <- mean(train$y)

  # ---- Euclidean distance field to spatial anchors ----
  anc_idx <- sample(seq_len(nrow(train)), min(n_anchor, nrow(train)))
  anchors <- cbind(train$x[anc_idx], train$y[anc_idx])
  edf <- matrix(0, nrow(query), nrow(anchors))
  for (a in seq_len(nrow(anchors))) {
    edf[, a] <- sqrt((query$x - anchors[a, 1])^2 + (query$y - anchors[a, 2])^2)
  }

  out <- cbind(Xc, lag, tlag, amean, edf, query$x, query$y)
  colnames(out) <- c(covars,
                     paste0("lag", radii), "tlag", "amean",
                     paste0("edf", seq_len(ncol(edf))), "x", "y")
  out
}
