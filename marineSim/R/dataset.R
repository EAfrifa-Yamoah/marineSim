#' Build one benchmark dataset with all feature blocks (Sections 3.2 and 3.3)
#'
#' Draws \code{n} sites under \code{design}, observes the response at every
#' timestep (optionally thinned at the training sites by a detection probability),
#' holds out \code{N_TEST} cells not used as sites, and constructs the four feature
#' blocks at the final timestep for both sets: coordinates (COORD), distance fields
#' to \code{min(N_ANCHORS, n)} k-means anchors of the sites (EDF), multi scale
#' spatial lags (LAG) and temporal features (TEMP: previous timestep response and
#' site history mean, read at the nearest training site). A same timestep nearest
#' neighbour response (NN0) is built as a control. When \code{n_time == 1} the
#' temporal block is constant (the placebo condition).
#'
#' Returns \code{NULL} when the training or test response is degenerate.
#'
#' @param w A world from \code{\link{make_world}}.
#' @param n Sites per timestep.
#' @param n_time Number of timesteps observed; the model is fitted at the last.
#' @param design Sampling design, see \code{\link{sample_sites}}.
#' @param seed Dataset seed.
#' @param detection Detection probability applied to training observations only.
#' @param response \code{"binary"} (occurrence) or \code{"abundance"} (Poisson counts).
#' @return A list with covariate matrices \code{Xc_tr}, \code{Xc_te}, responses
#'   \code{y_tr}, \code{y_te}, the true test probability or mean \code{p_te},
#'   block lists \code{btr}, \code{bte}, and the site, test and observation records.
#' @export
build_dataset <- function(w, n, n_time, design, seed, detection = 1, response = "binary") {
  set.seed(seed); te <- n_time
  site <- sample_sites(n, design); site_xy <- cell_xy(site)
  truth <- if (response == "binary") w$y else w$ycount
  obs <- vector("list", n_time)
  for (t in seq_len(n_time)) {
    yt <- as.numeric(truth[[t]][site])
    if (detection < 1) yt <- if (response == "binary") yt * (stats::runif(n) < detection)
                             else stats::rbinom(n, yt, detection)
    obs[[t]] <- yt
  }
  y_tr <- obs[[te]]
  if (response == "binary") { if (sum(y_tr) %in% c(0, n)) return(NULL) }
  else if (stats::sd(y_tr) == 0) return(NULL)
  pool <- setdiff(seq_len(GRID * GRID), site)
  test <- sample(pool, min(N_TEST, length(pool))); test_xy <- cell_xy(test)
  y_te <- as.numeric(truth[[te]][test])
  if (response == "binary") { if (sum(y_te) %in% c(0, length(y_te))) return(NULL) }
  else if (stats::sd(y_te) == 0) return(NULL)
  Xc_tr <- covariates(w, site, te); Xc_te <- covariates(w, test, te); gm <- mean(y_tr)
  btr <- bte <- list()
  btr$COORD <- site_xy; bte$COORD <- test_xy
  K <- min(N_ANCHORS, n)
  anchors <- stats::kmeans(site_xy, centers = K, nstart = 3, iter.max = 50)$centers
  btr$EDF <- edf(site_xy, anchors); bte$EDF <- edf(test_xy, anchors)
  btr$LAG <- spatial_lags(site_xy, y_tr, site_xy, TRUE)
  bte$LAG <- spatial_lags(site_xy, y_tr, test_xy, FALSE)
  if (te >= 2) {
    prev <- obs[[te - 1]]
    hist <- Reduce("+", obs[1:(te - 1)]) / (te - 1)
    j_tr <- nearest_train(site_xy, site_xy); j_te <- nearest_train(site_xy, test_xy)
    btr$TEMP <- cbind(prev[j_tr], hist[j_tr]); bte$TEMP <- cbind(prev[j_te], hist[j_te])
  } else {
    btr$TEMP <- matrix(gm, n, 2); bte$TEMP <- matrix(gm, length(test), 2)
  }
  btr$NN0 <- matrix(y_tr[nearest_train(site_xy, site_xy, exclude_self = TRUE)], ncol = 1)
  bte$NN0 <- matrix(y_tr[nearest_train(site_xy, test_xy)], ncol = 1)
  p_te <- as.numeric((if (response == "binary") w$p else w$mu)[[te]][test])   # true test probability / mean
  list(Xc_tr = Xc_tr, Xc_te = Xc_te, y_tr = y_tr, y_te = y_te, p_te = p_te, btr = btr, bte = bte,
       response = response, site = site, test = test, obs = obs)
}

#' Dataset seed used by the benchmark runner
#'
#' Combines the world seed with the sample size, number of timesteps, design
#' index, replicate, detection probability and response type, so that sharded
#' and single core runs produce identical rows.
#'
#' @param world_seed Seed of the world (see \code{\link{world_grid}}).
#' @param n,n_time Sample size and timesteps.
#' @param design_index 1 for random, 2 for clustered, 3 for stratified.
#' @param rep Replicate index (0 based).
#' @param detection Detection probability.
#' @param response \code{"binary"} or \code{"abundance"}.
#' @return An integer seed.
#' @export
dataset_seed <- function(world_seed, n, n_time, design_index, rep = 0, detection = 1,
                         response = "binary") {
  (world_seed * 31 + n * 131 + n_time * 7 + rep + 1009 * (design_index - 1) +
     as.integer(1000 * detection) + if (response == "abundance") 50000 else 0) %% 2147483647
}
