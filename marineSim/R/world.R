#' Build a ground truth world for the simulation benchmark
#'
#' Constructs a spatially and temporally structured occurrence process on a
#' continuous domain. Environmental covariates are deliberately rough (short
#' correlation range) so that a learner cannot use them as an implicit location
#' proxy; the smooth spatial signal lives in a separate latent field plus, in
#' the non stationary case, a spatially varying covariate effect. A persistent
#' temporal field with annual innovations introduces genuine temporal structure.
#'
#' @param world_id Integer world identifier (used to set the seed).
#' @param nonstationary Logical; if TRUE the depth effect varies over space.
#' @param spatial_range_km Correlation range of the latent spatial field.
#' @param n_time Number of time points.
#' @param grid_n Grid resolution for the underlying fields.
#' @param extent_km Domain side length in kilometres.
#' @param base_seed Integer seed offset.
#'
#' @return A list with the truth surfaces and a function \code{sample_sites}
#'   that draws observations under a chosen design, sample size and detection
#'   probability.
#' @export
make_world <- function(world_id = 0L, nonstationary = FALSE,
                       spatial_range_km = 80, n_time = 5L,
                       grid_n = 64L, extent_km = 400,
                       base_seed = 20260601L) {
  seed <- base_seed + world_id
  set.seed(seed)

  # rough environmental covariates (short ranges -> not a location proxy)
  depth <- matern_grf(grid_n, 45, extent_km, sigma = 1, seed = seed + 11)$field
  subst <- matern_grf(grid_n, 35, extent_km, sigma = 1, seed = seed + 12)$field
  wave  <- matern_grf(grid_n, 60, extent_km, sigma = 1, seed = seed + 13)$field
  sst   <- matern_grf(grid_n, 45, extent_km, sigma = 1, seed = seed + 14)$field

  # smooth latent spatial field (the structure spatial methods can exploit)
  re_field <- matern_grf(grid_n, spatial_range_km, extent_km,
                         sigma = 2.55, seed = seed + 21)$field

  # spatially varying coefficient for the non stationary regime
  svc <- if (nonstationary) {
    0.9 * matern_grf(grid_n, spatial_range_km, extent_km, sigma = 1,
                     seed = seed + 22)$field
  } else {
    matrix(0, grid_n, grid_n)
  }

  # persistent temporal fields with annual innovations
  persist_range <- 26
  t_fields <- vector("list", n_time)
  prev <- matern_grf(grid_n, persist_range, extent_km, sigma = 1.30,
                     seed = seed + 31)$field
  for (tt in seq_len(n_time)) {
    innov <- matern_grf(grid_n, persist_range, extent_km, sigma = 1.30,
                        seed = seed + 40 + tt)$field
    prev <- 0.75 * prev + sqrt(1 - 0.75^2) * innov
    t_fields[[tt]] <- prev
  }

  b <- list(b0 = 0.10, depth = -0.55, depth2 = -0.22, subst = 0.45,
            wave = -0.30, sst = -0.45, int = 0.22)

  linpred_grid <- function(tt) {
    eff_depth <- b$depth + svc
    lp <- b$b0 + eff_depth * depth + b$depth2 * depth^2 + b$subst * subst +
      b$wave * wave + b$sst * sst + b$int * (depth * wave) +
      re_field + t_fields[[tt]]
    lp
  }

  grf_like <- list(x = seq(0, extent_km, length.out = grid_n),
                   y = seq(0, extent_km, length.out = grid_n))

  sample_sites <- function(n = 100L, design = c("random", "clustered", "stratified"),
                           detection = 1.0, times = NULL) {
    design <- match.arg(design)
    if (is.null(times)) times <- seq_len(n_time)
    co <- draw_coords(n, design, extent_km)
    rows <- list(); r <- 1
    for (tt in times) {
      lpg <- list(field = linpred_grid(tt), x = grf_like$x, y = grf_like$y)
      lp <- sample_field(lpg, co$x, co$y)
      p <- 1 / (1 + exp(-lp))
      y_true <- stats::rbinom(length(p), 1, p)
      y_obs <- ifelse(y_true == 1 & stats::runif(length(p)) > detection, 0, y_true)
      gv <- function(f) sample_field(list(field = f, x = grf_like$x, y = grf_like$y),
                                     co$x, co$y)
      rows[[r]] <- data.frame(x = co$x, y = co$y, time = tt,
                              depth = gv(depth), subst = gv(subst),
                              wave = gv(wave), sst = gv(sst),
                              y = y_obs)
      r <- r + 1
    }
    do.call(rbind, rows)
  }

  list(seed = seed, nonstationary = nonstationary,
       spatial_range_km = spatial_range_km, n_time = n_time,
       extent_km = extent_km, sample_sites = sample_sites)
}

#' Draw site coordinates under a sampling design
#' @keywords internal
draw_coords <- function(n, design, extent_km) {
  if (design == "random") {
    x <- stats::runif(n, 0, extent_km); y <- stats::runif(n, 0, extent_km)
  } else if (design == "clustered") {
    nc <- max(2, round(n / 12))
    cx <- stats::runif(nc, 0, extent_km); cy <- stats::runif(nc, 0, extent_km)
    k <- sample(seq_len(nc), n, replace = TRUE)
    x <- pmin(pmax(cx[k] + stats::rnorm(n, 0, extent_km * 0.04), 0), extent_km)
    y <- pmin(pmax(cy[k] + stats::rnorm(n, 0, extent_km * 0.04), 0), extent_km)
  } else { # stratified over a coarse grid
    g <- ceiling(sqrt(n)); step <- extent_km / g
    gx <- rep(seq_len(g), times = g); gy <- rep(seq_len(g), each = g)
    sel <- sample(seq_along(gx), n)
    x <- (gx[sel] - 0.5) * step + stats::runif(n, -step / 2, step / 2)
    y <- (gy[sel] - 0.5) * step + stats::runif(n, -step / 2, step / 2)
    x <- pmin(pmax(x, 0), extent_km); y <- pmin(pmax(y, 0), extent_km)
  }
  list(x = x, y = y)
}
