#' Area under the ROC curve by the rank statistic
#' @param y Binary responses (0/1).
#' @param p Scores.
#' @return AUC in \code{[0, 1]}.
#' @export
auc <- function(y, p) { r <- rank(p); n1 <- sum(y == 1); n0 <- length(y) - n1
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }

#' Fit and score a random forest on covariates plus a coalition of blocks
#'
#' The forest backbone of the article: \pkg{ranger} with \code{RF_TREES} trees,
#' minimum node size 1, bootstrap with replacement, \code{mtry = floor(sqrt(p))}
#' for occurrence (probability forest, scored by AUC) and \code{mtry = p} on a
#' \code{log1p} response for abundance (scored by Spearman correlation with the
#' held out counts). The empty coalition is the standard RF, \code{"COORD"} alone
#' is the spatial RF, and \code{BLOCKS} is stRF.
#'
#' @param ds A dataset from \code{\link{build_dataset}}.
#' @param subset Character vector of block names to append to the covariates.
#' @param seed Forest seed.
#' @return A single score (AUC or Spearman rho).
#' @export
score <- function(ds, subset, seed) {
  Xtr <- do.call(cbind, c(list(ds$Xc_tr), ds$btr[subset]))
  Xte <- do.call(cbind, c(list(ds$Xc_te), ds$bte[subset]))
  colnames(Xtr) <- colnames(Xte) <- paste0("f", seq_len(ncol(Xtr)))
  p <- ncol(Xtr)
  if (ds$response == "binary") {
    dat <- data.frame(y = factor(ds$y_tr, levels = c(0, 1)), Xtr)
    fit <- ranger::ranger(y ~ ., data = dat, num.trees = RF_TREES, mtry = floor(sqrt(p)),
                          min.node.size = 1, probability = TRUE, seed = seed, num.threads = 1)
    pr <- stats::predict(fit, data.frame(Xte))$predictions[, "1"]
    return(auc(ds$y_te, pr))
  }
  dat <- data.frame(y = log1p(ds$y_tr), Xtr)
  fit <- ranger::ranger(y ~ ., data = dat, num.trees = RF_TREES, mtry = p, min.node.size = 1,
                        seed = seed, num.threads = 1)
  stats::cor(stats::predict(fit, data.frame(Xte))$predictions, ds$y_te, method = "spearman")
}

#' Score the covariate only GLM baseline
#'
#' Unpenalised logistic regression for occurrence (AUC) or Poisson regression
#' for abundance (Spearman rho), on the four covariates only.
#'
#' @inheritParams score
#' @return A single score.
#' @export
glm_score <- function(ds) {
  dtr <- data.frame(y = ds$y_tr, ds$Xc_tr); dte <- data.frame(ds$Xc_te)
  if (ds$response == "binary") {
    m <- suppressWarnings(stats::glm(y ~ ., data = dtr, family = stats::binomial()))
    return(auc(ds$y_te, stats::predict(m, dte, type = "response")))
  }
  m <- suppressWarnings(stats::glm(y ~ ., data = dtr, family = stats::poisson()))
  stats::cor(stats::predict(m, dte, type = "response"), ds$y_te, method = "spearman")
}
