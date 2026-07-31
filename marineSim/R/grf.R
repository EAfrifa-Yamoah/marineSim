#' Simulate a Matern Gaussian random field on a regular grid
#'
#' Generates a stationary, isotropic Gaussian random field with a Matern
#' covariance on a regular grid using spectral synthesis (a fast Fourier
#' transform of a filtered white noise field). This mirrors the generator used
#' in the Python simulation engine and avoids any heavy geostatistical
#' dependency.
#'
#' @param grid_n Integer grid side length (the field is grid_n by grid_n).
#' @param range_km Effective correlation range in kilometres.
#' @param extent_km Physical side length of the domain in kilometres.
#' @param nu Matern smoothness parameter (1.5 by default).
#' @param sigma Marginal standard deviation of the field.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with elements \code{field} (a grid_n by grid_n matrix),
#'   \code{x} and \code{y} (grid coordinate vectors in kilometres).
#' @export
matern_grf <- function(grid_n = 64, range_km = 80, extent_km = 400,
                       nu = 1.5, sigma = 1.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  # wave numbers
  k <- c(0:(grid_n / 2), (grid_n / 2 - 1):1) * (2 * pi / extent_km)
  kx <- matrix(k, grid_n, grid_n)
  ky <- t(kx)
  k2 <- kx^2 + ky^2
  # Matern spectral density (up to a constant); alpha = nu + d/2 with d = 2
  alpha <- nu + 1
  scale_len <- range_km / sqrt(8 * nu)
  spectrum <- (1 + (scale_len^2) * k2)^(-alpha)
  spectrum[1, 1] <- 0  # remove the mean component
  amp <- sqrt(spectrum)
  # complex white noise in the spectral domain
  white <- matrix(stats::rnorm(grid_n * grid_n), grid_n, grid_n) +
    1i * matrix(stats::rnorm(grid_n * grid_n), grid_n, grid_n)
  field <- Re(stats::fft(amp * white, inverse = TRUE)) / (grid_n)
  field <- (field - mean(field)) / stats::sd(as.vector(field)) * sigma
  x <- seq(0, extent_km, length.out = grid_n)
  list(field = field, x = x, y = x)
}

#' Bilinear sample of a gridded field at arbitrary coordinates
#' @keywords internal
sample_field <- function(grf, xq, yq) {
  gx <- grf$x; gy <- grf$y; F <- grf$field
  n <- length(gx); span <- gx[n] - gx[1]
  fx <- (xq - gx[1]) / span * (n - 1) + 1
  fy <- (yq - gy[1]) / span * (n - 1) + 1
  fx <- pmin(pmax(fx, 1), n); fy <- pmin(pmax(fy, 1), n)
  x0 <- floor(fx); y0 <- floor(fy)
  x1 <- pmin(x0 + 1, n); y1 <- pmin(y0 + 1, n)
  dx <- fx - x0; dy <- fy - y0
  idx <- function(i, j) F[cbind(i, j)]
  v <- idx(x0, y0) * (1 - dx) * (1 - dy) + idx(x1, y0) * dx * (1 - dy) +
    idx(x0, y1) * (1 - dx) * dy + idx(x1, y1) * dx * dy
  v
}
