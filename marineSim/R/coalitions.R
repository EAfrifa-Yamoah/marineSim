#' Coalitions fitted per dataset
#'
#' The 16 subsets of \code{BLOCKS} needed for an exact Shapley decomposition,
#' plus four control coalitions involving the same timestep nearest neighbour
#' response (NN0): \{NN0\}, \{NN0, TEMP\}, \{COORD, EDF, LAG, NN0\} and
#' \{COORD, EDF, LAG, NN0, TEMP\}.
#'
#' @return A list of character vectors.
#' @export
required_subsets <- function() {
  subs <- list(character(0))
  for (k in 1:4) subs <- c(subs, utils::combn(BLOCKS, k, simplify = FALSE))
  c(subs, list("NN0", c("NN0", "TEMP"), c("COORD", "EDF", "LAG", "NN0"),
               c("COORD", "EDF", "LAG", "NN0", "TEMP")))
}

#' Column name of a coalition's score
#' @param S Character vector of block names (empty for the covariate only forest).
#' @return A string such as \code{"v_COORD+LAG"} or \code{"v_BASE"}.
#' @export
skey <- function(S) paste0("v_", if (length(S)) paste(sort(S), collapse = "+") else "BASE")

#' Exact Shapley attribution of one dataset's coalition scores
#'
#' @param row A named list (or one row data frame coerced to list) holding the
#'   16 coalition scores named by \code{\link{skey}}.
#' @return A named numeric vector of Shapley values, one per block.
#' @export
shapley_row <- function(row) {
  k <- length(BLOCKS); phi <- stats::setNames(numeric(k), BLOCKS)
  for (b in BLOCKS) {
    others <- setdiff(BLOCKS, b); tot <- 0
    for (r in 0:length(others)) {
      Ss <- if (r == 0) list(character(0)) else utils::combn(others, r, simplify = FALSE)
      for (S in Ss) {
        w <- factorial(length(S)) * factorial(k - length(S) - 1) / factorial(k)
        tot <- tot + w * (row[[skey(c(S, b))]] - row[[skey(S)]])
      }
    }
    phi[b] <- tot
  }
  phi
}

#' Add Shapley values and derived quantities to a benchmark table
#'
#' Appends, for each dataset, the block Shapley values (\code{phi_COORD} etc.),
#' the algorithmic term (covariate only forest minus GLM), the total gain of the
#' full coalition over the GLM, the temporal increment when added last, the same
#' increment with the NN0 control present, the spatial total, and a
#' \code{field} identifier (world configuration by realisation) used as the
#' bootstrap cluster.
#'
#' @param d A data frame as written by the benchmark runner, read with
#'   \code{check.names = FALSE} so that \code{"v_COORD+EDF"} style names survive.
#' @return The augmented data frame.
#' @export
add_shapley <- function(d) {
  d$field <- paste0(d$world_cfg, "_r", d$realisation)
  sh <- t(apply(d, 1, function(r) shapley_row(as.list(sapply(r, function(v) suppressWarnings(as.numeric(v)))))))
  for (b in BLOCKS) d[[paste0("phi_", b)]] <- sh[, b]
  d$alg <- d[[skey(character(0))]] - d$GLM
  d$total <- d[[skey(BLOCKS)]] - d$GLM
  d$temp_added_last <- d[[skey(BLOCKS)]] - d[[skey(c("COORD", "EDF", "LAG"))]]
  d$temp_marg_withNN0 <- d[[skey(c("COORD", "EDF", "LAG", "NN0", "TEMP"))]] -
    d[[skey(c("COORD", "EDF", "LAG", "NN0"))]]
  d$spatial <- d$phi_COORD + d$phi_EDF + d$phi_LAG
  d
}

#' Cluster bootstrap interval for a mean
#'
#' Resamples the independent field realisations (column \code{field}) with
#' replacement \code{B} times and returns the mean with a 95 percent percentile
#' interval. The point estimate is the dataset level mean.
#'
#' @param d Data frame with a \code{field} column.
#' @param col Name of the column to summarise.
#' @param B Bootstrap replicates.
#' @param seed Seed.
#' @return A named vector \code{value}, \code{ci_lo}, \code{ci_hi}.
#' @export
boot_ci <- function(d, col, B = 2000, seed = 1) {
  set.seed(seed)
  m <- tapply(d[[col]], d$field, mean); G <- length(m)
  draws <- replicate(B, mean(m[sample.int(G, G, replace = TRUE)]))
  c(value = mean(d[[col]]), ci_lo = unname(stats::quantile(draws, 0.025)),
    ci_hi = unname(stats::quantile(draws, 0.975)))
}

#' Decomposition table for a subset of datasets
#'
#' @param d Data frame from \code{\link{add_shapley}}.
#' @param label Label for the subset.
#' @return A data frame with one row per component (algorithmic, the four blocks,
#'   spatial total, total gain), its bootstrap interval and share of the total.
#' @export
decomp_table <- function(d, label) {
  LAB <- c(COORD = "Coordinates", EDF = "Distance fields", LAG = "Spatial lags", TEMP = "Temporal")
  tot <- mean(d$total)
  comps <- c("Algorithmic flexibility" = "alg", stats::setNames(paste0("phi_", BLOCKS), LAB[BLOCKS]),
             "Spatial total" = "spatial", "Total gain over GLM" = "total")
  do.call(rbind, lapply(names(comps), function(nm) {
    ci <- boot_ci(d, comps[[nm]])
    data.frame(subset = label, component = nm, value = ci[["value"]], ci_lo = ci[["ci_lo"]],
               ci_hi = ci[["ci_hi"]], share_pct = 100 * ci[["value"]] / tot,
               n_datasets = nrow(d), n_fields = length(unique(d$field)))
  }))
}
