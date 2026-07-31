#' Fit one of the benchmark models and predict on new data
#'
#' Thin wrappers around the five models compared in the benchmark. The ladder is
#' deliberate: \code{GLM} uses covariates only; \code{RF} adds nonparametric
#' flexibility but still no coordinates; \code{SpatialRF} adds raw coordinates;
#' \code{stRF} adds the full spatiotemporal feature layer. \code{GAM} is a spline
#' baseline with a spatial tensor product. Isolating each rung lets the gain be
#' decomposed into algorithmic, spatial and temporal contributions.
#'
#' @param method One of "GLM", "GAM", "RF", "SpatialRF", "stRF".
#' @param train Training data frame (with response column \code{y}).
#' @param test Test data frame.
#' @param covars Covariate column names.
#' @param num_trees Number of trees for the random forest learners.
#'
#' @return Numeric vector of predicted occurrence probabilities for \code{test}.
#' @export
fit_model <- function(method = c("GLM", "GAM", "RF", "SpatialRF", "stRF"),
                      train, test, covars = c("depth", "subst", "wave", "sst"),
                      num_trees = 150L) {
  method <- match.arg(method)
  ytr <- train$y

  if (method == "GLM") {
    f <- stats::as.formula(paste("y ~", paste(covars, collapse = " + ")))
    m <- stats::glm(f, data = train, family = stats::binomial())
    return(stats::predict(m, newdata = test, type = "response"))
  }

  if (method == "GAM") {
    sm <- paste(sprintf("s(%s)", covars), collapse = " + ")
    f <- stats::as.formula(paste("y ~ te(x, y) +", sm))
    m <- mgcv::gam(f, data = train, family = stats::binomial(), method = "REML")
    p <- as.numeric(stats::predict(m, newdata = test, type = "response"))
    return(pmin(pmax(p, 1e-4), 1 - 1e-4))
  }

  build <- function(df) {
    if (method == "RF") {
      as.matrix(df[, covars, drop = FALSE])
    } else if (method == "SpatialRF") {
      as.matrix(df[, c(covars, "x", "y"), drop = FALSE])
    } else {
      NULL  # stRF handled separately
    }
  }

  if (method == "stRF") {
    Xtr <- strf_features(train, train, covars = covars)
    Xte <- strf_features(train, test, covars = covars)
  } else {
    Xtr <- build(train); Xte <- build(test)
  }

  dat <- data.frame(Xtr, y = factor(ytr, levels = c(0, 1)))
  rf <- ranger::ranger(y ~ ., data = dat, num.trees = num_trees,
                       probability = TRUE, min.node.size = 5,
                       respect.unordered.factors = "order")
  pr <- stats::predict(rf, data = data.frame(Xte))$predictions
  as.numeric(pr[, "1"])
}
