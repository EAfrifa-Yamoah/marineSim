#' Ground truth ecosystem generator (Section 3.1)
#'
#' Builds one simulated world: four standardised Matérn covariate fields
#' (depth, substrate and wave, static, with ranges 45, 35 and 60 km; sea surface
#' temperature as an AR(1) process in time with \eqn{\rho = 0.6} and spatially
#' correlated innovations of range 45 km), a latent spatial field at the world's
#' range with sd 2.55, a static persistent site effect (26 km, sd 1.30),
#' optionally non stationary depth and SST coefficients perturbed by
#' \eqn{0.20 \times} a 250 km field, and a habitat mask of structural zeros
#' covering the lower quartile of a 120 km field. The linear predictor is
#' \deqn{\eta_t = 0.10 + \beta_d z_d - 0.22 z_d^2 + 0.45 z_s - 0.30 z_w + 0.22 z_d z_w + \beta_T z_{T,t} + \xi + \pi,}
#' with \eqn{\beta_d = -0.55} and \eqn{\beta_T = -0.45} in the stationary case.
#' Occurrence is Bernoulli with probability \eqn{\mathrm{mask} \cdot \mathrm{logit}^{-1}(\eta_t)};
#' abundance is Poisson with mean \eqn{\mathrm{mask} \cdot \exp(\log 3 + \eta_t / 2)}.
#'
#' Fields are stored as \code{GRID} by \code{GRID} matrices. Cell index \code{i}
#' in \code{1:GRID^2} is column major: \code{row = (i - 1) \%\% GRID + 1},
#' \code{col = (i - 1) \%/\% GRID + 1}, and the cell centre is at
#' \code{x = (col - 0.5) * CELL_KM}, \code{y = (row - 0.5) * CELL_KM}.
#'
#' @param stationary Logical; \code{FALSE} gives spatially varying depth and SST coefficients.
#' @param range_km Correlation range of the latent spatial field (30, 80 or 200 in the article).
#' @param n_time Number of timesteps to simulate.
#' @param seed Integer seed; see \code{\link{world_grid}} for the article's seeds.
#' @return A list with the covariate fields (\code{depth}, \code{subst}, \code{wave},
#'   \code{sst} as a list over time), the truth surfaces per timestep (\code{eta},
#'   \code{p}, \code{y}, \code{mu}, \code{ycount}), the coefficient fields
#'   (\code{b_depth}, \code{b_sst}) and the latent components (\code{spatial_re},
#'   \code{persist}).
#' @export
make_world <- function(stationary, range_km, n_time = 5L, seed = 0L) {
  set.seed(seed); g <- GRID
  depth <- matern_grf(g, 45 / CELL_KM); subst <- matern_grf(g, 35 / CELL_KM)
  wave  <- matern_grf(g, 60 / CELL_KM)
  sst <- vector("list", n_time); sst[[1]] <- matern_grf(g, 45 / CELL_KM); rho <- 0.6
  if (n_time > 1) for (t in 2:n_time)
    sst[[t]] <- rho * sst[[t - 1]] + sqrt(1 - rho^2) * matern_grf(g, 45 / CELL_KM)
  spatial_re <- 2.55 * matern_grf(g, range_km / CELL_KM)
  persist <- 1.30 * matern_grf(g, 26 / CELL_KM)          # static across time
  if (stationary) {
    b_depth <- matrix(-0.55, g, g); b_sst <- matrix(-0.45, g, g)
  } else {
    b_depth <- -0.55 + 0.20 * matern_grf(g, 250 / CELL_KM)
    b_sst   <- -0.45 + 0.20 * matern_grf(g, 250 / CELL_KM)
  }
  hab <- matern_grf(g, 120 / CELL_KM)
  struct <- (hab > stats::quantile(hab, 0.25)) * 1
  eta <- p <- y <- mu <- ycount <- vector("list", n_time)
  for (t in seq_len(n_time)) {
    eta[[t]] <- 0.10 + b_depth * depth - 0.22 * depth^2 + 0.45 * subst - 0.30 * wave +
      0.22 * depth * wave + b_sst * sst[[t]] + spatial_re + persist
    p[[t]] <- struct * stats::plogis(eta[[t]])
    y[[t]] <- matrix(stats::rbinom(g * g, 1, as.vector(p[[t]])), g, g)
    mu[[t]] <- struct * exp(C0_COUNT + ETA_SCALE * eta[[t]])
    ycount[[t]] <- matrix(stats::rpois(g * g, as.vector(mu[[t]])), g, g)
  }
  list(g = g, n_time = n_time, depth = depth, subst = subst, wave = wave, sst = sst,
       p = p, y = y, mu = mu, ycount = ycount, eta = eta, b_sst = b_sst, b_depth = b_depth,
       persist = persist, spatial_re = spatial_re)
}

#' Cell centres in kilometres
#' @param idx Integer cell indices in \code{1:GRID^2}.
#' @return A two column matrix of \code{x}, \code{y} coordinates (km).
#' @export
cell_xy <- function(idx) {
  g <- GRID; r <- (idx - 1) %% g + 1; cc <- (idx - 1) %/% g + 1
  cbind(x = (cc - 0.5) * CELL_KM, y = (r - 0.5) * CELL_KM)
}

#' Cell index from row and column, clamped to the grid
#' @param r,cc Integer row and column vectors.
#' @return Integer cell indices.
#' @export
cell_idx <- function(r, cc) (pmin(pmax(cc, 1), GRID) - 1) * GRID + pmin(pmax(r, 1), GRID)

#' Covariate matrix at a set of cells and a timestep
#' @param w A world from \code{\link{make_world}}.
#' @param idx Integer cell indices.
#' @param t Timestep.
#' @return A matrix with columns \code{depth}, \code{subst}, \code{wave}, \code{sst}.
#' @export
covariates <- function(w, idx, t)
  cbind(depth = w$depth[idx], subst = w$subst[idx], wave = w$wave[idx], sst = w$sst[[t]][idx])

#' The six world configurations of the article and their seeds
#'
#' Two stationarity regimes crossed with three latent ranges (30, 80, 200 km).
#' The world seed is \code{BASE_SEED + 1000 * realisation + world_id}, where
#' \code{world_id} runs 0 to 5 in the order stationary 30/80/200, then non
#' stationary 30/80/200.
#'
#' @param realisations Integer vector of field realisation indices (0 to 5 in the article).
#' @return A list of configuration lists with elements \code{wid}, \code{cfg},
#'   \code{real}, \code{stat}, \code{range_km}, \code{seed}.
#' @export
world_grid <- function(realisations) {
  out <- list(); wid <- 0L
  for (stat in c(TRUE, FALSE)) for (rng_km in c(30, 80, 200)) {
    cfg <- sprintf("%s_r%d", if (stat) "stat" else "nonstat", rng_km)
    for (r in realisations)
      out[[length(out) + 1]] <- list(wid = wid, cfg = cfg, real = r, stat = stat,
                                     range_km = rng_km, seed = BASE_SEED + 1000L * r + wid)
    wid <- wid + 1L
  }
  out
}
