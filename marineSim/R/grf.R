#' Discrete Fourier frequencies
#'
#' Equivalent of \code{numpy.fft.fftfreq(n)} with unit sample spacing.
#' @param n Integer length.
#' @return Numeric vector of length \code{n}.
#' @keywords internal
#' @export
fftfreq <- function(n) c(0:floor((n - 1) / 2), (-floor(n / 2)):-1) / n

#' Standardised Matérn Gaussian random field by spectral synthesis
#'
#' Simulates a stationary Gaussian random field with a Matérn covariance of
#' smoothness \code{nu} on a \code{g} by \code{g} grid, by filtering complex white
#' noise with the square root of the Matérn spectral density on a doubled grid
#' (to suppress periodic wrap around) and keeping the first \code{g} rows and
#' columns. The field is returned centred and scaled to unit standard deviation.
#'
#' The function consumes \code{2 * (2g)^2} normal draws from the current random
#' stream, so the sequence of calls inside \code{\link{make_world}} fixes the
#' world for a given seed.
#'
#' @param g Grid size (cells per side).
#' @param range_cells Correlation range expressed in cells.
#' @param nu Matérn smoothness (1.5 throughout the article).
#' @return A \code{g} by \code{g} numeric matrix with mean 0 and sd 1.
#' @export
matern_grf <- function(g, range_cells, nu = 1.5) {
  M <- 2L * g
  ky <- fftfreq(M); kx <- fftfreq(M)
  k2 <- outer(ky^2, kx^2, "+")
  alpha <- sqrt(2 * nu) / max(range_cells, 1e-6)
  spec <- (alpha^2 + (2 * pi)^2 * k2)^(-(nu + 1))
  spec[1, 1] <- 0
  noise <- matrix(stats::rnorm(M * M), M, M) + 1i * matrix(stats::rnorm(M * M), M, M)
  field <- Re(stats::fft(stats::fft(noise) * sqrt(spec), inverse = TRUE)) / (M * M)
  field <- field[1:g, 1:g]
  (field - mean(field)) / (stats::sd(as.vector(field)) + 1e-12)
}
