#' Draw survey sites under a sampling design (Section 3.2)
#'
#' Returns \code{n} distinct grid cells.
#' \describe{
#'   \item{\code{"random"}}{Simple random sampling of cells.}
#'   \item{\code{"clustered"}}{Preferential sampling from three Gaussian clusters:
#'     the cluster centres are uniform on the central 70 percent of each axis and
#'     the spread is 6 percent of the grid side (about 24 km). Cells are drawn
#'     without replacement with probability proportional to the mixture intensity
#'     at the cell centre, so that exactly \code{n} distinct cells are returned and
#'     the design remains clustered at every \code{n}.}
#'   \item{\code{"stratified"}}{The grid is partitioned into
#'     \eqn{\lceil\sqrt{n}\rceil^2} equal blocks, \code{n} blocks are chosen without
#'     replacement, and one cell is drawn uniformly from each chosen block.}
#'   \item{\code{"clustered_legacy"}, \code{"stratified_legacy"}}{The generators of
#'     marineSim 3.0.x, retained so that those results can be reproduced. They
#'     drew continuous positions, rounded them to cells, discarded duplicates and
#'     topped the sample up with uniformly random cells; for the clustered design
#'     the top up reached 25 percent of sites at \code{n = 200} and 45 percent at
#'     \code{n = 500}.}
#' }
#'
#' @param n Number of sites.
#' @param design One of \code{"random"}, \code{"clustered"}, \code{"stratified"},
#'   \code{"clustered_legacy"}, \code{"stratified_legacy"}.
#' @return Integer vector of \code{n} distinct cell indices.
#' @export
sample_sites <- function(n, design) {
  g <- GRID
  if (design == "random") return(sample.int(g * g, n))
  if (design == "clustered") {
    ctr <- matrix(stats::runif(6, 0.15, 0.85) * g, 3, 2)           # (row, col) centres in cell units
    idx <- seq_len(g * g); r <- (idx - 1) %% g + 0.5; cc <- (idx - 1) %/% g + 0.5
    s2 <- 2 * (0.06 * g)^2
    w <- exp(-((r - ctr[1, 1])^2 + (cc - ctr[1, 2])^2) / s2) +
         exp(-((r - ctr[2, 1])^2 + (cc - ctr[2, 2])^2) / s2) +
         exp(-((r - ctr[3, 1])^2 + (cc - ctr[3, 2])^2) / s2)
    return(sample.int(g * g, n, prob = w))
  }
  if (design == "stratified") {
    s <- ceiling(sqrt(n))
    idx <- seq_len(g * g); r <- (idx - 1) %% g + 0.5; cc <- (idx - 1) %/% g + 0.5
    block <- (pmin(floor(r / g * s), s - 1)) * s + pmin(floor(cc / g * s), s - 1) + 1   # block of each cell centre
    chosen <- sample.int(s * s, n)
    return(vapply(chosen, function(b) { cells <- idx[block == b]; if (length(cells) == 1) cells else sample(cells, 1) }, integer(1)))
  }
  if (design == "clustered_legacy") {
    ctr <- matrix(stats::runif(6, 0.15, 0.85) * g, 3, 2)
    per <- rep(n %/% 3, 3); per[seq_len(n %% 3)] <- per[seq_len(n %% 3)] + 1
    cells <- integer(0)
    for (k in 1:3) {
      r  <- round(ctr[k, 1] + stats::rnorm(per[k], 0, g * 0.06))
      cc <- round(ctr[k, 2] + stats::rnorm(per[k], 0, g * 0.06))
      cells <- c(cells, cell_idx(r, cc))
    }
    cells <- unique(cells)
    while (length(cells) < n) cells <- unique(c(cells, sample.int(g * g, n - length(cells))))
    return(cells[1:n])
  }
  if (design == "stratified_legacy") {
    s <- ceiling(sqrt(n)); edges <- seq(0, g, length.out = s + 1)
    strata <- expand.grid(i = 1:s, j = 1:s)
    strata <- strata[sample.int(nrow(strata), n), ]
    r  <- floor(stats::runif(n, edges[strata$i], edges[strata$i + 1])) + 1
    cc <- floor(stats::runif(n, edges[strata$j], edges[strata$j + 1])) + 1
    cells <- unique(cell_idx(r, cc))
    while (length(cells) < n) cells <- unique(c(cells, sample.int(g * g, n - length(cells))))
    return(cells[1:n])
  }
  stop("unknown design: ", design)
}
