#!/usr/bin/env Rscript
# Seagrass latitudinal extrapolation case study in R.
# Expects the Posidonia sinuosa dataset at data/Bayesiandataset_2025_final.csv
# (columns as described in the manuscript). Produces leave-one-region-out skill
# and an area-of-applicability map along the latitudinal transect.
suppressPackageStartupMessages({
  library(marineSim)
})

dat <- read.csv("data/Bayesiandataset_2025_final.csv")
covars <- c("depth", "kd490", "htyperank",
            "meansummermIST", "maxsummermIST", "CVsummermIST")

# aggregate to site-year means
agg <- aggregate(
  cbind(dens = rawcounts, latitude, longitude,
        depth, kd490, htyperank,
        meansummermIST, maxsummermIST, CVsummermIST) ~ Region_5 + Site + YearSeagrassSampling,
  data = dat, FUN = function(z) mean(z, na.rm = TRUE))
names(agg)[names(agg) == "YearSeagrassSampling"] <- "year"
names(agg)[names(agg) == "Region_5"] <- "region"
agg$x <- agg$longitude; agg$y <- agg$latitude; agg$time <- agg$year

regions <- c("JBMP", "MMP", "CSMC", "SIMP", "NCMP")
loro <- data.frame()
for (held in regions) {
  tr <- agg[agg$region != held, ]; te <- agg[agg$region == held, ]
  for (m in c("GLM", "RF", "SpatialRF", "stRF")) {
    # regression variants: reuse the feature layer, fit a ranger regression
    if (m == "stRF") {
      Xtr <- strf_features(tr, tr, covars = covars)
      Xte <- strf_features(tr, te, covars = covars)
    } else if (m == "SpatialRF") {
      Xtr <- as.matrix(tr[, c(covars, "x", "y")]); Xte <- as.matrix(te[, c(covars, "x", "y")])
    } else if (m == "RF") {
      Xtr <- as.matrix(tr[, covars]); Xte <- as.matrix(te[, covars])
    } else {
      Xtr <- as.matrix(tr[, covars]); Xte <- as.matrix(te[, covars])
    }
    ytr <- log1p(tr$dens)
    if (m == "GLM") {
      fit <- lm(ytr ~ ., data = data.frame(Xtr))
      pr <- predict(fit, newdata = data.frame(Xte))
    } else {
      rf <- ranger::ranger(y ~ ., data = data.frame(Xtr, y = ytr),
                           num.trees = 400, min.node.size = 3)
      pr <- predict(rf, data = data.frame(Xte))$predictions
    }
    yhat <- expm1(pr)
    rho <- suppressWarnings(cor(te$dens, yhat, method = "spearman"))
    rmse <- sqrt(mean((te$dens - yhat)^2))
    loro <- rbind(loro, data.frame(held = held, method = m,
                                   Spearman = rho, RMSE = rmse))
  }
}
print(loro)

# area of applicability along a latitudinal transect
aoa_vars <- c(covars, "longitude", "latitude")
site_ref <- aggregate(. ~ site_region, data = transform(
  agg[, aoa_vars], site_region = paste(agg$region, agg$Site)), FUN = mean)
imp <- ranger::ranger(y ~ ., data = data.frame(agg[, aoa_vars], y = log1p(agg$dens)),
                      num.trees = 300, importance = "permutation")$variable.importance
grid <- data.frame(latitude = seq(min(agg$latitude) - 0.15, max(agg$latitude) + 0.15,
                                  length.out = 160))
for (cc in aoa_vars) grid[[cc]] <- approx(agg$latitude, agg[[cc]], grid$latitude,
                                          rule = 2)$y
aoa <- area_of_applicability(agg[, aoa_vars], grid[, aoa_vars], importance = imp)
cat(sprintf("Inside area of applicability: %.0f%%\n", 100 * mean(aoa$inside)))
