#' Draw survey sites under a sampling design (Section 3.2)
#'
#' Returns \code{n} distinct grid cells. \code{"random"} is simple random sampling
#' of cells. \code{"clustered"} places three Gaussian clusters whose centres are
#' uniform on the central 70 percent of each axis and whose spread is 6 percent
#' of the grid, filling any shortfall from duplicate cells at random.
#' \code{"stratified"} partitions the grid into \eqn{\lceil\sqrt{n}\rceil^2} equal
#' blocks, chooses \code{n} of them without replacement and draws one cell
#' uniformly within each.
#'
#' @param n Number of sites.
#' @param design One of \code{"random"}, \code{"clustered"}, \code{"stratified"}.
#' @return Integer vector of \code{n} distinct cell indices.
#' @export
sample_sites <- function(n, design) {
  g <- GRID
  if (design == "random") return(sample.int(g * g, n))
  if (design == "clustered") {
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
  s <- ceiling(sqrt(n)); edges <- seq(0, g, length.out = s + 1)
  strata <- expand.grid(i = 1:s, j = 1:s)
  strata <- strata[sample.int(nrow(strata), n), ]
  r  <- floor(stats::runif(n, edges[strata$i], edges[strata$i + 1])) + 1
  cc <- floor(stats::runif(n, edges[strata$j], edges[strata$j + 1])) + 1
  cells <- unique(cell_idx(r, cc))
  while (length(cells) < n) cells <- unique(c(cells, sample.int(g * g, n - length(cells))))
  cells[1:n]
}
