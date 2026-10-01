# =============================================================================
# 08_case_study.R — Posidonia sinuosa case study (Section 4 and Supplement S7 to S9)
# Mirrors the corrected Python pipeline (case_study_corrected*.py, case_maps.py,
# make_temporal_hovmoller.py) choice for choice: aggregation and imputation, the four
# feature blocks (IDW same-year lags with fallbacks, k = 8 k-means distance fields,
# previous-year regional mean), ranger with 400 trees / min node 3 / mtry = 0.6 p, log1p
# response, three validation regimes, permutation and drop-column importance, coastal
# corridor with leaf raw-scale prediction and 80% bands, static and temporal AOA,
# rolling-origin forecasts, and interval coverage.
# =============================================================================
suppressPackageStartupMessages({library(ranger); library(mgcv); library(jsonlite); library(sp)})
OUT <- Sys.getenv("CASE_OUT", "../results/case_study"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
DATA <- Sys.getenv("CASE_DATA", "../data/raw/Bayesiandataset_2025_final.csv")
COVARS <- c("depth", "kd490", "htyperank", "meansummermIST", "maxsummermIST", "CVsummermIST")
TVAR <- c("kd490", "meansummermIST", "maxsummermIST", "CVsummermIST")
RADII <- c(25, 75, 150); REGIONS <- c("JBMP", "MMP", "CSMC", "SIMP", "NCMP")
BLOCKS <- c("COORD", "EDF", "LAG", "TEMP"); FULL <- BLOCKS
NTREE <- 400; MINNODE <- 3; MTRYF <- 0.6; SEED <- 7

# ----------------------------------------------------------------- data
load_case <- function() {
  raw <- read.csv(DATA)
  agg <- aggregate(cbind(dens = rawcounts, latitude, longitude, depth, kd490, htyperank,
                         meansummermIST, maxsummermIST, CVsummermIST) ~ Region_5 + Site + YearSeagrassSampling,
                   data = raw, FUN = mean, na.action = na.pass)
  names(agg)[1:3] <- c("region", "site", "year")
  for (c in COVARS) {
    m <- ave(agg[[c]], agg$region, FUN = function(v) mean(v, na.rm = TRUE))
    agg[[c]][is.na(agg[[c]])] <- m[is.na(agg[[c]])]
    agg[[c]][is.na(agg[[c]])] <- mean(agg[[c]], na.rm = TRUE)
  }
  agg <- agg[!is.na(agg$dens), ]; agg$site_uid <- paste0(agg$region, ":", agg$site)
  agg$year <- as.integer(agg$year); rownames(agg) <- NULL; agg
}
hav <- function(lon1, lat1, lon2, lat2) {
  R <- 6371; p1 <- lat1 * pi / 180; p2 <- lat2 * pi / 180
  a <- sin((lat2 - lat1) * pi / 360)^2 + cos(p1) * cos(p2) * sin((lon2 - lon1) * pi / 360)^2
  2 * R * asin(sqrt(pmin(1, a)))
}
# ----------------------------------------------------------------- feature blocks
lag_block <- function(tr, q) {
  f <- matrix(0, nrow(q), length(RADII))
  for (i in seq_len(nrow(q))) {
    d <- hav(q$longitude[i], q$latitude[i], tr$longitude, tr$latitude); same <- tr$year == q$year[i]
    for (j in seq_along(RADII)) {
      done <- FALSE
      for (mask in list(same, rep(TRUE, length(d)))) {
        sel <- mask & d <= RADII[j] & d > 1e-6
        if (sum(sel) >= 2) { w <- 1 / (d[sel] + 1); f[i, j] <- sum(w * tr$dens[sel]) / sum(w); done <- TRUE; break }
      }
      if (!done) { k <- min(5, sum(d > 1e-6)); f[i, j] <- mean(tr$dens[order(d)[1:k]]) }
    }
  }
  f
}
edf_block <- function(tr, q, anchors = NULL) {
  if (is.null(anchors)) { set.seed(0); anchors <- kmeans(cbind(tr$longitude, tr$latitude), centers = min(8, nrow(tr)), nstart = 5, iter.max = 300)$centers }
  D <- sapply(seq_len(nrow(anchors)), function(a) hav(q$longitude, q$latitude, anchors[a, 1], anchors[a, 2]))
  list(D = matrix(D, nrow = nrow(q)), anchors = anchors)
}
temp_block <- function(tr, q) {
  ry <- aggregate(dens ~ region + year, tr, mean); gm <- mean(tr$dens)
  key <- paste(ry$region, ry$year); rm <- tapply(tr$dens, tr$region, mean)
  v <- numeric(nrow(q))
  for (i in seq_len(nrow(q))) {
    k <- paste(q$region[i], q$year[i] - 1); j <- match(k, key)
    v[i] <- if (!is.na(j)) ry$dens[j] else if (!is.na(rm[q$region[i]])) rm[q$region[i]] else gm
  }
  matrix(v, ncol = 1)
}
build_blocks <- function(tr, q) {
  e <- edf_block(tr, tr); eq <- edf_block(tr, q, e$anchors)
  list(tr = list(COORD = cbind(tr$longitude, tr$latitude), EDF = e$D, LAG = lag_block(tr, tr), TEMP = temp_block(tr, tr)),
       q  = list(COORD = cbind(q$longitude, q$latitude), EDF = eq$D, LAG = lag_block(tr, q), TEMP = temp_block(tr, q)))
}
design <- function(tr, q, b, subset) {
  s <- sort(subset)
  X <- do.call(cbind, c(list(as.matrix(tr[COVARS])), b$tr[s])); Q <- do.call(cbind, c(list(as.matrix(q[COVARS])), b$q[s]))
  nm <- c(COVARS, unlist(lapply(s, function(k) paste0(k, seq_len(ncol(b$tr[[k]]))))))
  colnames(X) <- colnames(Q) <- nm; list(X = X, Q = Q)
}
rf_fit <- function(X, ylog, importance = "none", quantreg = FALSE, keep = FALSE) {
  ranger(x = X, y = ylog, num.trees = NTREE, min.node.size = MINNODE, mtry = max(1, floor(MTRYF * ncol(X))),
         seed = SEED, num.threads = 1, importance = importance, quantreg = quantreg, keep.inbag = keep)
}
metrics <- function(y, yhat) { ok <- is.finite(yhat); y <- y[ok]; yhat <- yhat[ok]
  c(R2 = 1 - sum((y - yhat)^2) / sum((y - mean(y))^2), Spearman = cor(y, yhat, method = "spearman"), RMSE = sqrt(mean((y - yhat)^2))) }
fit_glm <- function(tr, q) { d <- data.frame(y = log1p(tr$dens), scale(tr[COVARS])); m <- lm(y ~ ., d)
  qq <- data.frame(scale(q[COVARS], center = attr(scale(tr[COVARS]), "scaled:center"), scale = attr(scale(tr[COVARS]), "scaled:scale")))
  expm1(predict(m, qq)) }
fit_gam <- function(tr, q) {
  d <- data.frame(y = log1p(tr$dens), tr[COVARS], lon = tr$longitude, lat = tr$latitude)
  m <- suppressWarnings(gam(y ~ s(depth, k = 5) + s(kd490, k = 5) + htyperank + s(meansummermIST, k = 5) +
                              s(maxsummermIST, k = 5) + s(CVsummermIST, k = 5) + te(lon, lat, k = c(4, 4)), data = d, method = "REML"))
  p <- predict(m, data.frame(q[COVARS], lon = q$longitude, lat = q$latitude))
  expm1(pmin(pmax(p, min(d$y)), max(d$y)))          # clipped to the training range (Section 5.6)
}
leaf_raw_mean <- function(m, X, raw, Q) {          # raw-scale conditional mean from the leaves
  Ltr <- predict(m, X, type = "terminalNodes")$predictions; Lq <- predict(m, Q, type = "terminalNodes")$predictions
  out <- numeric(nrow(Q))
  for (b in seq_len(ncol(Ltr))) { lm_ <- tapply(raw, Ltr[, b], mean); out <- out + lm_[as.character(Lq[, b])] }
  as.numeric(out / ncol(Ltr))
}

# ----------------------------------------------------------------- D. validation regimes
cv_regimes <- function(agg) {
  ladder <- list(RF = character(0), SpatialRF = "COORD", stRF = FULL); rows <- list()
  run <- function(splits, label, per_region = FALSE) {
    store <- list(); obs <- c(); reg <- c(); preds <- setNames(vector("list", 5), c(names(ladder), "GLM", "GAM"))
    for (sp_ in splits) {
      tr <- agg[sp_$tr, ]; te <- agg[sp_$te, ]; b <- build_blocks(tr, te); ylog <- log1p(tr$dens)
      for (nm in names(ladder)) { D <- design(tr, te, b, ladder[[nm]]); preds[[nm]] <- c(preds[[nm]], expm1(predict(rf_fit(D$X, ylog), D$Q)$predictions)) }
      preds$GLM <- c(preds$GLM, fit_glm(tr, te)); preds$GAM <- c(preds$GAM, fit_gam(tr, te))
      obs <- c(obs, te$dens); reg <- c(reg, te$region)
    }
    for (nm in names(preds)) {
      m <- metrics(obs, preds[[nm]])
      within <- if (per_region) mean(sapply(REGIONS, function(r) cor(obs[reg == r], preds[[nm]][reg == r], method = "spearman"))) else NA
      rows[[length(rows) + 1]] <<- data.frame(regime = label, method = nm, R2 = m[["R2"]], Spearman_pooled = m[["Spearman"]], RMSE = m[["RMSE"]], Spearman_within_park = within)
      if (per_region) for (r in REGIONS) { mm <- metrics(obs[reg == r], preds[[nm]][reg == r])
        rows[[length(rows) + 1]] <<- data.frame(regime = paste0("LORO:", r), method = nm, R2 = mm[["R2"]], Spearman_pooled = mm[["Spearman"]], RMSE = mm[["RMSE"]], Spearman_within_park = NA) }
    }
  }
  set.seed(1); f <- sample(rep(1:10, length.out = nrow(agg)))
  run(lapply(1:10, function(k) list(tr = which(f != k), te = which(f == k))), "random 10-fold")
  sites <- unique(agg$site_uid); set.seed(1); sf <- setNames(sample(rep(1:10, length.out = length(sites))), sites); g <- sf[agg$site_uid]
  run(lapply(1:10, function(k) list(tr = which(g != k), te = which(g == k))), "site-grouped")
  run(lapply(REGIONS, function(r) list(tr = which(agg$region != r), te = which(agg$region == r))), "LORO", per_region = TRUE)
  do.call(rbind, rows)
}

# ----------------------------------------------------------------- importance and drop-column
importance_study <- function(agg) {
  b <- build_blocks(agg, agg); D <- design(agg, agg, b, FULL)
  m <- rf_fit(D$X, log1p(agg$dens), importance = "permutation")
  imp <- m$variable.importance; data.frame(feature = names(imp), importance = as.numeric(imp), share = as.numeric(imp) / sum(pmax(imp, 0)))
}
dropcol_study <- function(agg) {
  set.seed(1); f <- sample(rep(1:10, length.out = nrow(agg)))
  groups <- list(full = character(0), longitude = "COORD1", latitude = "COORD2", both_coords = c("COORD1", "COORD2"),
                 spatial_block = c("COORD", "EDF", "LAG"), depth = "depth", kd490 = "kd490", htyperank = "htyperank",
                 meansummermIST = "meansummermIST", maxsummermIST = "maxsummermIST", CVsummermIST = "CVsummermIST", temporal = "TEMP")
  res <- list()
  for (g in names(groups)) {
    preds <- c(); obs <- c()
    for (k in 1:10) { tr <- agg[f != k, ]; te <- agg[f == k, ]; b <- build_blocks(tr, te); D <- design(tr, te, b, FULL)
      drop <- groups[[g]]; keep <- !(colnames(D$X) %in% drop) & !grepl(paste0("^(", paste(drop, collapse = "|"), ")\\d*$"), colnames(D$X))
      if (!length(drop)) keep <- rep(TRUE, ncol(D$X))
      preds <- c(preds, expm1(predict(rf_fit(D$X[, keep, drop = FALSE], log1p(tr$dens)), D$Q[, keep, drop = FALSE])$predictions)); obs <- c(obs, te$dens) }
    res[[g]] <- metrics(obs, preds)[["Spearman"]]
  }
  data.frame(removed = names(res), Spearman = unlist(res), loss = res$full - unlist(res))
}

# ----------------------------------------------------------------- corridor, static AOA, transect
corridor <- function(agg) {
  site <- aggregate(cbind(latitude, longitude, dens, depth, kd490, htyperank, meansummermIST, maxsummermIST, CVsummermIST) ~ region + site_uid, agg, mean)
  o <- order(site$latitude); sp_ <- lowess(site$latitude[o], site$longitude[o], f = 0.5)
  spine <- function(lat) approx(sp_$x, sp_$y, lat, rule = 2)$y
  lats <- seq(min(site$latitude) - 0.05, max(site$latitude) + 0.05, by = 0.022); cells <- list()
  for (la in lats) { c0 <- spine(la); for (lo in seq(c0 - 0.24, c0 + 0.24 + 1e-9, by = 0.022)) cells[[length(cells) + 1]] <- c(lo, la) }
  g <- as.data.frame(do.call(rbind, cells)); names(g) <- c("longitude", "latitude")
  coast <- fromJSON(Sys.getenv("CASE_COAST", "../data/coast_wa.json"), simplifyVector = FALSE)
  onland <- rep(FALSE, nrow(g))
  for (poly in coast) { P <- do.call(rbind, lapply(poly, unlist)); if (nrow(P) < 3) next
    if (max(P[, 1]) < min(g$longitude) - 0.5 || min(P[, 1]) > max(g$longitude) + 0.5 || max(P[, 2]) < min(g$latitude) - 0.5 || min(P[, 2]) > max(g$latitude) + 0.5) next
    onland <- onland | point.in.polygon(g$longitude, g$latitude, P[, 1], P[, 2]) > 0 }
  g <- g[!onland, ]
  nn <- sapply(seq_len(nrow(g)), function(i) which.min(hav(g$longitude[i], g$latitude[i], site$longitude, site$latitude)))
  g$region <- site$region[nn]; g$year <- max(agg$year)
  idw <- function(col, k = 4) sapply(seq_len(nrow(g)), function(i) { d <- hav(g$longitude[i], g$latitude[i], site$longitude, site$latitude); j <- order(d)[1:k]; w <- 1 / (d[j]^2 + 1e-6); sum(w * site[[col]][j]) / sum(w) })
  for (c in COVARS) g[[c]] <- idw(c)
  # stRF on all data; leaf raw mean; 80% band via quantile forest
  b <- build_blocks(agg, g); D <- design(agg, g, b, FULL); ylog <- log1p(agg$dens)
  m <- rf_fit(D$X, ylog, importance = "impurity", quantreg = TRUE)
  g$pred <- leaf_raw_mean(m, D$X, agg$dens, D$Q)
  qq <- predict(m, D$Q, type = "quantiles", quantiles = c(0.1, 0.9))$predictions; g$lo <- expm1(qq[, 1]); g$hi <- expm1(qq[, 2])
  # static AOA on site means (Meyer & Pebesma), importance weighted
  aoa_vars <- c(COVARS, "longitude", "latitude"); ref <- as.matrix(site[aoa_vars]); mu <- colMeans(ref); sd_ <- apply(ref, 2, sd) + 1e-12
  imp_m <- rf_fit(ref, log1p(site$dens), importance = "impurity"); w <- sqrt(pmax(imp_m$variable.importance[aoa_vars], 0))
  Z <- sweep(sweep(ref, 2, mu), 2, sd_, "/") * rep(w, each = nrow(ref)); Zq <- sweep(sweep(as.matrix(g[aoa_vars]), 2, mu), 2, sd_, "/") * rep(w, each = nrow(g))
  Dm <- as.matrix(dist(Z)); dbar <- mean(Dm[upper.tri(Dm)]); diag(Dm) <- Inf; thr <- quantile(apply(Dm, 1, min) / dbar, 0.95)
  dq <- sapply(seq_len(nrow(g)), function(i) min(sqrt(colSums((t(Z) - Zq[i, ])^2)))) / dbar
  g$DI <- dq; g$inside <- dq <= thr
  list(grid = g, site = site, DI_thresh = as.numeric(thr), model = m, X = D$X, band = qq)
}

# ----------------------------------------------------------------- temporal reconstruction + temporal AOA
temporal_field <- function(agg, cor_) {
  g <- cor_$grid; site <- cor_$site; years <- min(agg$year):max(agg$year)
  # E1: site-grouped 5-fold OOF regional trajectories
  sites <- unique(agg$site_uid); set.seed(5); sf <- setNames(sample(rep(1:5, length.out = length(sites))), sites); f <- sf[agg$site_uid]
  oof <- numeric(nrow(agg))
  for (k in 1:5) { tr <- agg[f != k, ]; te <- agg[f == k, ]; b <- build_blocks(tr, te); D <- design(tr, te, b, FULL)
    m <- rf_fit(D$X, log1p(tr$dens)); oof[f == k] <- leaf_raw_mean(m, D$X, tr$dens, D$Q) }
  traj <- aggregate(cbind(obs = dens, pred = oof) ~ region + year, data.frame(agg, oof = oof), mean)
  oof_stats <- c(Spearman = cor(agg$dens, oof, method = "spearman"), bias = mean(oof - agg$dens))
  # E4: per-year corridor prediction; temporal AOA referenced to site-years, between-site calibration
  b <- build_blocks(agg, agg); D <- design(agg, agg, b, FULL); m <- rf_fit(D$X, log1p(agg$dens))
  aoa_vars <- c(COVARS, "longitude", "latitude"); ref <- as.matrix(agg[aoa_vars]); mu <- colMeans(ref); sd_ <- apply(ref, 2, sd) + 1e-12
  imp_m <- rf_fit(ref, log1p(agg$dens), importance = "impurity"); w <- sqrt(pmax(imp_m$variable.importance[aoa_vars], 0))
  Z <- sweep(sweep(ref, 2, mu), 2, sd_, "/") * rep(w, each = nrow(ref))
  Dm <- as.matrix(dist(Z)); same <- outer(agg$site_uid, agg$site_uid, "=="); Dm[same] <- Inf
  nn <- apply(Dm, 1, min); dbar <- mean(nn); thr <- quantile(nn / dbar, 0.95)
  per <- split(agg, agg$site_uid); fields <- list(); dis <- list()
  for (Y in years) {
    q <- g[c("longitude", "latitude", "region", "depth", "htyperank")]; q$year <- Y
    for (c in TVAR) { sv <- sapply(site$site_uid, function(s) { sy <- per[[s]]; j <- if (Y %in% sy$year) which(sy$year == Y)[1] else which.min(abs(sy$year - Y)); sy[[c]][j] })
      q[[c]] <- sapply(seq_len(nrow(q)), function(i) { d <- hav(q$longitude[i], q$latitude[i], site$longitude, site$latitude); j <- order(d)[1:6]; wt <- 1 / (d[j]^2 + 1e-6); sum(wt * sv[j]) / sum(wt) }) }
    bq <- list(tr = b$tr, q = list(COORD = cbind(q$longitude, q$latitude), EDF = edf_block(agg, q, edf_block(agg, agg)$anchors)$D, LAG = lag_block(agg, q), TEMP = temp_block(agg, q)))
    Dq <- design(agg, q, bq, FULL); fields[[as.character(Y)]] <- leaf_raw_mean(m, D$X, agg$dens, Dq$Q)
    Zq <- sweep(sweep(as.matrix(q[aoa_vars]), 2, mu), 2, sd_, "/") * rep(w, each = nrow(q))
    dis[[as.character(Y)]] <- sapply(seq_len(nrow(q)), function(i) min(sqrt(colSums((t(Z) - Zq[i, ])^2)))) / dbar
  }
  nb <- 70; edges <- seq(min(g$latitude), max(g$latitude), length.out = nb + 1); bidx <- pmin(pmax(findInterval(g$latitude, edges), 1), nb)
  M <- Mdi <- matrix(NA, nb, length(years)); dep <- -g$depth
  for (yi in seq_along(years)) for (bb in 1:nb) { s <- bidx == bb; if (any(s)) { M[bb, yi] <- mean(fields[[yi]][s]); Mdi[bb, yi] <- mean(dis[[yi]][s]) } }
  depth_prof <- data.frame(band = 1:nb, lat = (edges[-1] + edges[-(nb + 1)]) / 2, dep_mean = tapply(dep, factor(bidx, 1:nb), mean), dep_lo = tapply(dep, factor(bidx, 1:nb), min), dep_hi = tapply(dep, factor(bidx, 1:nb), max))
  list(traj = traj, oof_stats = oof_stats, M = M, Mdi = Mdi, years = years, lat_centres = depth_prof$lat, DI_thresh = as.numeric(thr), depth_prof = depth_prof,
       outside_frac = mean(unlist(dis) > thr), fields = fields, dis = dis)
}

# ----------------------------------------------------------------- rolling origin
rolling_origin <- function(agg, horizons = c(1, 2, 3, 5)) {
  yrs <- sort(unique(agg$year)); origins <- yrs[sapply(yrs, function(y) sum(agg$year <= y) >= 200) & yrs <= yrs[length(yrs) - 1]]; rows <- list()
  for (o in origins) { tr <- agg[agg$year <= o, ]; clim <- tapply(tr$dens, tr$site_uid, mean); last <- tapply(tr$dens[order(tr$year)], tr$site_uid[order(tr$year)], function(v) v[length(v)])
    for (h in horizons) { te <- agg[agg$year == o + h, ]; te <- te[te$site_uid %in% names(clim), ]; if (nrow(te) < 15) next
      b <- build_blocks(tr, te); D <- design(tr, te, b, FULL); p_strf <- expm1(predict(rf_fit(D$X, log1p(tr$dens)), D$Q)$predictions)
      base <- clim[te$site_uid]; y <- te$dens
      for (nm in c("stRF", "climatology", "persistence")) { p <- switch(nm, stRF = p_strf, climatology = clim[te$site_uid], persistence = last[te$site_uid]); ok <- is.finite(p)
        rows[[length(rows) + 1]] <- data.frame(origin = o, horizon = h, method = nm, n = sum(ok), rho_level = cor(y[ok], p[ok], method = "spearman"),
          rho_anomaly = if (nm == "stRF") cor(y[ok] - base[ok], p[ok] - base[ok], method = "spearman") else NA, rmse = sqrt(mean((y[ok] - p[ok])^2))) } } }
  do.call(rbind, rows)
}

# ----------------------------------------------------------------- intervals + AOA coverage
interval_study <- function(agg, alpha = 0.8) {
  rows <- list(); set.seed(2); f5 <- sample(rep(1:5, length.out = nrow(agg)))
  splits <- c(lapply(1:5, function(k) list(lab = "interpolation (random 5-fold)", tr = which(f5 != k), te = which(f5 == k))),
              lapply(REGIONS, function(r) list(lab = "extrapolation (LORO)", tr = which(agg$region != r), te = which(agg$region == r))))
  for (s in splits) { tr <- agg[s$tr, ]; te <- agg[s$te, ]; if (nrow(te) < 10) next
    b <- build_blocks(tr, te); D <- design(tr, te, b, FULL); ylog <- log1p(tr$dens)
    m <- rf_fit(D$X, ylog, importance = "impurity", quantreg = TRUE)
    pt <- predict(m, D$Q, predict.all = TRUE)$predictions          # per-tree
    lo_e <- expm1(apply(pt, 1, quantile, 0.1)); hi_e <- expm1(apply(pt, 1, quantile, 0.9))
    qq <- predict(m, D$Q, type = "quantiles", quantiles = c(0.1, 0.9))$predictions; lo_q <- expm1(qq[, 1]); hi_q <- expm1(qq[, 2])
    set.seed(3); idx <- sample(nrow(tr)); ncal <- max(30, nrow(tr) %/% 4); cal <- tr[idx[1:ncal], ]; fit <- tr[idx[-(1:ncal)], ]
    bc <- build_blocks(fit, cal); Dc <- design(fit, cal, bc, FULL); m2 <- rf_fit(Dc$X, log1p(fit$dens))
    qres <- quantile(abs(log1p(cal$dens) - predict(m2, Dc$Q)$predictions), alpha)
    bt <- build_blocks(fit, te); Dt <- design(fit, te, bt, FULL); ctr <- predict(m2, Dt$Q)$predictions; lo_c <- expm1(ctr - qres); hi_c <- expm1(ctr + qres)
    # AOA on the training design (site-years) with impurity weights
    mu <- colMeans(D$X); sd_ <- apply(D$X, 2, sd) + 1e-12; w <- sqrt(pmax(m$variable.importance, 0))
    Z <- sweep(sweep(D$X, 2, mu), 2, sd_, "/") * rep(w, each = nrow(D$X)); Zq <- sweep(sweep(D$Q, 2, mu), 2, sd_, "/") * rep(w, each = nrow(D$Q))
    Dm <- as.matrix(dist(Z)); dbar <- mean(Dm[upper.tri(Dm)]); diag(Dm) <- Inf; thr <- quantile(apply(Dm, 1, min) / dbar, 0.95)
    di <- sapply(seq_len(nrow(Zq)), function(i) min(sqrt(colSums((t(Z) - Zq[i, ])^2)))) / dbar; inside <- di <= thr; y <- te$dens
    for (nm in c("ensemble spread", "QRF", "split conformal")) { lo <- switch(nm, `ensemble spread` = lo_e, QRF = lo_q, `split conformal` = lo_c); hi <- switch(nm, `ensemble spread` = hi_e, QRF = hi_q, `split conformal` = hi_c)
      cv <- y >= lo & y <= hi
      rows[[length(rows) + 1]] <- data.frame(regime = s$lab, interval = nm, coverage = mean(cv), coverage_inside = if (any(inside)) mean(cv[inside]) else NA,
        coverage_outside = if (any(!inside)) mean(cv[!inside]) else NA, width = mean(hi - lo), width_inside = if (any(inside)) mean((hi - lo)[inside]) else NA,
        width_outside = if (any(!inside)) mean((hi - lo)[!inside]) else NA, pct_outside = mean(!inside)) } }
  do.call(rbind, rows)
}

# ----------------------------------------------------------------- leakage check
leakage_check <- function(agg) {
  tr <- agg[agg$region != "CSMC", ]; te <- agg[agg$region == "CSMC", ]; b1 <- build_blocks(tr, te)
  te2 <- te; set.seed(0); te2$dens <- sample(te2$dens) * 3 + 11; b2 <- build_blocks(tr, te2)
  same <- all(sapply(BLOCKS, function(k) isTRUE(all.equal(b1$q[[k]], b2$q[[k]]))))
  j <- which.min(hav(te$longitude[1], te$latitude[1], tr$longitude, tr$latitude))   # a training row that IS a neighbour of the first query
  tr2 <- tr; tr2$dens[j] <- tr2$dens[j] + 50; b3 <- build_blocks(tr2, te)
  list(no_leakage = same, non_vacuous = !isTRUE(all.equal(b1$q$LAG, b3$q$LAG)))
}

# ================================================================= main
if (sys.nframe() == 0) {
  t0 <- Sys.time(); agg <- load_case()
  cat(sprintf("data: %d site-years, %d sites, %d coordinate pairs, %d-%d\n", nrow(agg), length(unique(agg$site_uid)), nrow(unique(agg[c("latitude", "longitude")])), min(agg$year), max(agg$year)))
  write.csv(agg, file.path(OUT, "agg.csv"), row.names = FALSE)
  cv <- cv_regimes(agg); write.csv(cv, file.path(OUT, "cv_regimes.csv"), row.names = FALSE); cat("cv regimes done\n"); print(cv[cv$regime %in% c("random 10-fold", "site-grouped", "LORO"), ], digits = 3)
  imp <- importance_study(agg); write.csv(imp, file.path(OUT, "importance.csv"), row.names = FALSE); cat("importance done\n")
  dc <- dropcol_study(agg); write.csv(dc, file.path(OUT, "dropcol.csv"), row.names = FALSE); cat("drop-column done\n"); print(dc, digits = 3)
  cor_ <- corridor(agg); write.csv(cor_$grid, file.path(OUT, "corridor.csv"), row.names = FALSE); write.csv(cor_$site, file.path(OUT, "sites.csv"), row.names = FALSE)
  cat(sprintf("corridor: %d cells, %.1f%% inside AOA, DI thresh %.2f\n", nrow(cor_$grid), 100 * mean(cor_$grid$inside), cor_$DI_thresh))
  tf <- temporal_field(agg, cor_); write.csv(tf$traj, file.path(OUT, "traj.csv"), row.names = FALSE)
  write.csv(data.frame(band = rep(1:70, length(tf$years)), lat = rep(tf$lat_centres, length(tf$years)), year = rep(tf$years, each = 70), pred = as.vector(tf$M), DI = as.vector(tf$Mdi)), file.path(OUT, "hovmoller.csv"), row.names = FALSE)
  write.csv(tf$depth_prof, file.path(OUT, "depth_profile.csv"), row.names = FALSE)
  cat(sprintf("temporal: OOF Spearman %.3f bias %.2f | temporal AOA thresh %.2f, %.1f%% outside\n", tf$oof_stats[1], tf$oof_stats[2], tf$DI_thresh, 100 * tf$outside_frac))
  ro <- rolling_origin(agg); write.csv(ro, file.path(OUT, "rolling_origin.csv"), row.names = FALSE); cat("rolling origin done\n")
  iv <- interval_study(agg); write.csv(iv, file.path(OUT, "intervals.csv"), row.names = FALSE); cat("intervals done\n")
  lk <- leakage_check(agg)
  meta <- list(n_siteyears = nrow(agg), n_sites = length(unique(agg$site_uid)), n_coords = nrow(unique(agg[c("latitude", "longitude")])),
               corridor_cells = nrow(cor_$grid), static_inside_pct = 100 * mean(cor_$grid$inside), static_DI_thresh = cor_$DI_thresh,
               temporal_DI_thresh = tf$DI_thresh, temporal_outside_pct = 100 * tf$outside_frac, oof_spearman = tf$oof_stats[[1]], oof_bias = tf$oof_stats[[2]],
               leakage = lk, runtime_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  write(toJSON(meta, auto_unbox = TRUE, pretty = TRUE, digits = 6), file.path(OUT, "meta.json"))
  cat("CASE STUDY COMPLETE\n")
}
