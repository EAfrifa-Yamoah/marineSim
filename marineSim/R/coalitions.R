#' Feature-block coalitions and Shapley decomposition of predictive skill
#'
#' The published decomposition differenced a nested ladder of models
#' (GLM -> RF -> SpatialRF -> stRF). Because the engineered feature blocks
#' encode overlapping information, the increment attributed to any one block
#' depends on the order in which blocks enter, and that dependence spans an
#' order of magnitude for the temporal block. These functions fit every
#' coalition of the four blocks on a common random forest backbone and return
#' an order-free Shapley attribution.
#'
#' Blocks:
#'   COORD  raw coordinates (x, y)
#'   EDF    Euclidean distance fields to k-means anchors
#'   LAG    multi-scale spatial lags of the response
#'   TEMP   previous-timestep response at the nearest site, and site history mean
#'
#' Environmental covariates are always present and are not part of the
#' attribution.
#'
#' NOTE: these functions mirror python/ablation.py and python/decompose.py.
#' They were written alongside the Python implementation but have not been
#' executed, because the revision analyses were run through the Python engine.
#' Run marineSim's testthat suite before relying on them.
#'
#' @name coalitions
NULL

BLOCKS <- c("COORD", "EDF", "LAG", "TEMP")


#' Enumerate all subsets of a character vector
#'
#' @param blocks character vector of block names
#' @return list of character vectors, including the empty set
#' @export
all_coalitions <- function(blocks = BLOCKS) {
  out <- list(character(0))
  for (k in seq_along(blocks)) {
    cmb <- utils::combn(blocks, k, simplify = FALSE)
    out <- c(out, cmb)
  }
  out
}


#' Assemble a design matrix from covariates plus a chosen set of blocks
#'
#' @param covariates numeric matrix of environmental covariates
#' @param blocks named list of numeric matrices, one per block
#' @param subset character vector naming the blocks to include
#' @return numeric matrix
#' @export
assemble_design <- function(covariates, blocks, subset) {
  stopifnot(is.matrix(covariates))
  parts <- list(covariates)
  for (b in sort(subset)) {
    if (is.null(blocks[[b]])) {
      stop(sprintf("block '%s' not supplied", b))
    }
    parts[[length(parts) + 1L]] <- as.matrix(blocks[[b]])
  }
  do.call(cbind, parts)
}


#' Evaluate every feature-block coalition on a common forest backbone
#'
#' @param covariates_train,covariates_test covariate matrices
#' @param blocks_train,blocks_test named lists of block matrices
#' @param y_train training response
#' @param y_test test response, used only for scoring
#' @param metric function(observed, predicted) returning a scalar skill score;
#'   defaults to AUC for a binary response
#' @param num.trees,min.node.size passed to ranger
#' @param seed random seed
#' @return data.frame with one row per coalition
#' @export
coalition_skill <- function(covariates_train, covariates_test,
                            blocks_train, blocks_test,
                            y_train, y_test,
                            metric = auc_binary,
                            num.trees = 500, min.node.size = 1,
                            seed = 1L) {
  subsets <- all_coalitions()
  res <- vector("list", length(subsets))
  for (i in seq_along(subsets)) {
    S <- subsets[[i]]
    Xtr <- assemble_design(covariates_train, blocks_train, S)
    Xte <- assemble_design(covariates_test, blocks_test, S)
    df <- data.frame(.y = y_train, Xtr)
    set.seed(seed)
    fit <- ranger::ranger(
      dependent.variable.name = ".y", data = df,
      num.trees = num.trees, min.node.size = min.node.size,
      probability = is.factor(y_train), seed = seed
    )
    pr <- stats::predict(fit, data = data.frame(Xte))$predictions
    if (is.matrix(pr)) pr <- pr[, ncol(pr)]
    res[[i]] <- data.frame(
      coalition = if (length(S)) paste(sort(S), collapse = "+") else "BASE",
      k = length(S),
      skill = metric(y_test, pr),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, res)
}


#' Shapley value of each feature block
#'
#' @param skill data.frame from \code{coalition_skill}
#' @param blocks character vector of block names
#' @return named numeric vector of Shapley values
#' @export
shapley_blocks <- function(skill, blocks = BLOCKS) {
  v <- stats::setNames(skill$skill, skill$coalition)
  key <- function(S) if (length(S)) paste(sort(S), collapse = "+") else "BASE"
  k <- length(blocks)
  phi <- stats::setNames(numeric(k), blocks)
  for (b in blocks) {
    others <- setdiff(blocks, b)
    total <- 0
    for (r in 0:length(others)) {
      cmb <- if (r == 0) list(character(0)) else
        utils::combn(others, r, simplify = FALSE)
      for (S in cmb) {
        w <- factorial(length(S)) * factorial(k - length(S) - 1) / factorial(k)
        total <- total + w * (v[[key(c(S, b))]] - v[[key(S)]])
      }
    }
    phi[[b]] <- total
  }
  phi
}


#' Sequential (order-dependent) increments along a nested ladder
#'
#' Retained for comparison with the Shapley attribution. The spread between
#' orderings is itself a reportable quantity.
#'
#' @param skill data.frame from \code{coalition_skill}
#' @param order character vector giving the order in which blocks enter
#' @return named numeric vector of increments
#' @export
sequential_blocks <- function(skill, order = BLOCKS) {
  v <- stats::setNames(skill$skill, skill$coalition)
  key <- function(S) if (length(S)) paste(sort(S), collapse = "+") else "BASE"
  cur <- character(0)
  out <- stats::setNames(numeric(length(order)), order)
  for (b in order) {
    out[[b]] <- v[[key(c(cur, b))]] - v[[key(cur)]]
    cur <- c(cur, b)
  }
  out
}


#' Placebo check for the temporal block
#'
#' With a single timestep the temporal features are constant by construction
#' and can carry no information, so any honest estimator must return
#' approximately zero. The nested contrast fails this check; the isolated
#' increment passes it.
#'
#' @param skill data.frame from \code{coalition_skill} fitted at one timestep
#' @return list with the nested contrast and the isolated increment
#' @export
temporal_placebo <- function(skill) {
  v <- stats::setNames(skill$skill, skill$coalition)
  list(
    nested_contrast = unname(v[["EDF+LAG+TEMP"]] - v[["COORD"]]),
    isolated_increment = unname(v[["COORD+EDF+LAG+TEMP"]] - v[["COORD+EDF+LAG"]])
  )
}


#' Assert that held-out responses cannot influence any predictor
#'
#' Implements the reviewer's invariance test: perturbing the test response must
#' leave every training and query feature unchanged, while perturbing the
#' training response must move the response-derived blocks (otherwise the test
#' is vacuous).
#'
#' @param builder function(y_train, y_test) returning a named list with
#'   \code{train} and \code{test} block lists
#' @param y_train,y_test responses
#' @return list with \code{no_leakage} and \code{non_vacuous}
#' @export
assert_no_leakage <- function(builder, y_train, y_test) {
  base <- builder(y_train, y_test)
  perturbed_test <- builder(y_train, sample(y_test) + 1)
  same <- all(mapply(function(a, b) isTRUE(all.equal(a, b)),
                     unlist(base, recursive = FALSE),
                     unlist(perturbed_test, recursive = FALSE)))
  perturbed_train <- builder(y_train + 1, y_test)
  moved <- !isTRUE(all.equal(base$test$LAG, perturbed_train$test$LAG))
  list(no_leakage = isTRUE(same), non_vacuous = isTRUE(moved))
}


#' Area under the ROC curve for a binary response
#'
#' @param obs binary observations
#' @param pred predicted scores
#' @export
auc_binary <- function(obs, pred) {
  obs <- as.numeric(as.character(obs))
  r <- rank(pred)
  n1 <- sum(obs == 1); n0 <- sum(obs == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[obs == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
