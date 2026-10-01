# Case study figures from results/case_study/: Figures 6, 7, 8 and S5, S6, S7
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figstyle.R"))
suppressPackageStartupMessages(library(jsonlite))
O <- function(f) res("case_study", f)
agg <- read.csv(O("agg.csv")); site <- read.csv(O("sites.csv")); grid <- read.csv(O("corridor.csv"))
cv <- read.csv(O("cv_regimes.csv")); imp <- read.csv(O("importance.csv")); dc <- read.csv(O("dropcol.csv"))
coast <- fromJSON(Sys.getenv("CASE_COAST", res("..", "data", "coast_wa.json")), simplifyVector = FALSE)
MKEY <- c(GLM = "GLM", GAM = "GAM", RF = "Standard RF", SpatialRF = "Spatial RF", stRF = "stRF")
FEAT_COL <- c("spatial feature (b, c)" = OI[["blue"]], "temporal (b, c)" = OI[["vermilion"]], "environmental (b, c)" = OI[["grey"]])

# ------------------------------------------------------------------ Figure 6
regimes <- c("random 10-fold" = "random\n10-fold", "site-grouped" = "site\ngrouped", "LORO" = "leave one\nregion out")
cva <- cv[cv$regime %in% names(regimes) & cv$method %in% names(MKEY), ]
cva$value <- ifelse(cva$regime == "LORO" & !is.na(cva$Spearman_within_park), cva$Spearman_within_park, cva$Spearman_pooled)
cva$regime <- factor(regimes[cva$regime], levels = regimes); cva$method <- factor(MKEY[cva$method], levels = MKEY)
pa <- ggplot(cva, aes(regime, value, fill = method)) + geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) + scale_fill_manual(values = METHOD_COL) +
  labs(x = NULL, y = "Spearman rank correlation", fill = NULL) + theme(axis.text.x = element_text(size = 7.5))
grp <- list("longitude" = "COORD1", "latitude" = "COORD2", "spatial lags" = paste0("LAG", 1:3), "distance fields" = paste0("EDF", 1:8),
            "temporal" = "TEMP1", "depth" = "depth", "light (kd490)" = "kd490", "habitat rank" = "htyperank",
            "thermal (3)" = c("meansummermIST", "maxsummermIST", "CVsummermIST"))
sh <- data.frame(group = names(grp), share = sapply(grp, function(v) 100 * sum(imp$share[imp$feature %in% v])))
sh$kind <- ifelse(sh$group %in% c("longitude", "latitude", "spatial lags", "distance fields"), "spatial feature (b, c)",
                  ifelse(sh$group == "temporal", "temporal (b, c)", "environmental (b, c)"))
sh <- sh[order(-sh$share), ]; sh$group <- factor(sh$group, levels = rev(sh$group))
pb <- ggplot(sh, aes(share, group, fill = kind)) + geom_col(width = 0.75) + scale_fill_manual(values = FEAT_COL) +
  labs(x = "Permutation importance (% of total)", y = NULL, fill = NULL) + theme(axis.text.y = element_text(size = 7.5))
lab <- c(longitude = "longitude", latitude = "latitude", both_coords = "both coordinates", spatial_block = "all spatial features", depth = "depth",
         kd490 = "light (kd490)", htyperank = "habitat rank", meansummermIST = "mean summer SST", maxsummermIST = "max summer SST",
         CVsummermIST = "CV summer SST", temporal = "temporal")
d2 <- dc[dc$removed != "full", ]; d2 <- d2[order(-d2$loss), ]
d2$kind <- ifelse(d2$removed %in% c("longitude", "latitude", "both_coords", "spatial_block"), "spatial feature (b, c)",
                  ifelse(d2$removed == "temporal", "temporal (b, c)", "environmental (b, c)"))
d2$label <- factor(lab[d2$removed], levels = rev(lab[d2$removed]))
pc <- ggplot(d2, aes(loss, label, fill = kind)) + geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) + geom_col(width = 0.75) +
  scale_fill_manual(values = FEAT_COL) + labs(x = "Loss in CV rank correlation when removed", y = NULL, fill = NULL) +
  theme(axis.text.y = element_text(size = 7.5))
fig6 <- (pa + pb + pc) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom", legend.text = element_text(size = 7.5))
save_png(fig6, "main/Figure6.png", 11, 4.2)

# ------------------------------------------------------------------ Figure 7
coast_df <- do.call(rbind, lapply(seq_along(coast), function(i) { P <- do.call(rbind, lapply(coast[[i]], unlist))
  if (nrow(P) < 3) return(NULL); data.frame(id = i, lon = P[, 1], lat = P[, 2]) }))
lon0 <- 114.75; lon1 <- 116.05; lat0 <- min(agg$latitude) - 0.12; lat1 <- max(agg$latitude) + 0.12
latbr <- seq(-33.5, -30, 0.5)
base_map <- function() list(
  geom_polygon(data = coast_df, aes(lon, lat, group = id), fill = "#e9e4d8", colour = "grey55", linewidth = 0.2),
  coord_cartesian(xlim = c(lon0, lon1), ylim = c(lat0, lat1), expand = FALSE),
  scale_y_continuous(breaks = latbr, labels = sprintf("%.1f", abs(latbr))),
  labs(x = "Longitude (°E)", y = "Latitude (°S)"),
  theme(axis.text = element_text(size = 7), panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.4), axis.line = element_blank()))
site$region <- factor(site$region, levels = REG)
lab7 <- aggregate(cbind(longitude, latitude) ~ region, site, mean)
p7a <- ggplot() + base_map() +
  geom_point(data = site, aes(longitude, latitude, size = dens, fill = region), shape = 21, stroke = 0.3) +
  geom_text(data = lab7, aes(pmax(lon0 + 0.03, longitude - 0.42), latitude, label = region), size = 2.6, fontface = "bold", hjust = 0) +
  scale_fill_manual(values = REG_COL, guide = "none") + scale_size_area(max_size = 3.2, guide = "none")
vlo <- floor(min(site$dens, grid$pred)); vhi <- ceiling(max(quantile(site$dens, 0.98), max(grid$pred)))
# isobaths: the corridor cells lie on a regular lattice. Cells outside the corridor are filled with the
# depth of the nearest corridor cell so that contours are not drawn along the corridor edge, and the
# resulting lines are clipped back to the corridor (vertices within 0.8 of a lattice step of a corridor cell).
lons <- sort(unique(round(grid$longitude, 6))); lats <- sort(unique(round(grid$latitude, 6)))
sx <- cos(mean(grid$latitude) * pi / 180); step <- min(diff(lats))
Z <- matrix(NA_real_, length(lons), length(lats))
Z[cbind(match(round(grid$longitude, 6), lons), match(round(grid$latitude, 6), lats))] <- grid$depth
na <- which(is.na(Z), arr.ind = TRUE)
for (k in seq_len(nrow(na))) {
  d2 <- ((lons[na[k, 1]] - grid$longitude) * sx)^2 + (lats[na[k, 2]] - grid$latitude)^2
  Z[na[k, 1], na[k, 2]] <- grid$depth[which.min(d2)]
}
cl <- grDevices::contourLines(lons, lats, Z, levels = c(-15, -10, -5))
iso <- do.call(rbind, lapply(seq_along(cl), function(i) data.frame(id = i, level = cl[[i]]$level, lon = cl[[i]]$x, lat = cl[[i]]$y)))
iso$keep <- sapply(seq_len(nrow(iso)), function(i) min(((iso$lon[i] - grid$longitude) * sx)^2 + (iso$lat[i] - grid$latitude)^2) <= (0.8 * step)^2)
iso$id <- cumsum(c(TRUE, diff(iso$id) != 0) | !c(TRUE, iso$keep[-nrow(iso)]))   # break paths where vertices are dropped
iso <- iso[iso$keep, ]
iso_lab <- do.call(rbind, lapply(split(iso, iso$level), function(s) { s <- s[order(s$lat), ]; k <- round(nrow(s) * c(0.25, 0.75)); s[k[k > 0], ] }))
iso_lab$txt <- paste0(abs(iso_lab$level), " m")
sq <- 1.55
p7b <- ggplot() + base_map() +
  geom_point(data = grid, aes(longitude, latitude, colour = pred), shape = 15, size = sq) +
  geom_path(data = iso, aes(lon, lat, group = id), colour = "grey15", linewidth = 0.3, alpha = 0.85) +
  geom_label(data = iso_lab, aes(lon, lat, label = txt), size = 1.9, label.size = 0, label.padding = unit(0.08, "lines"), fill = "white", alpha = 0.8) +
  geom_point(data = site, aes(longitude, latitude, fill = dens, shape = "Observed site"), size = 2.4, stroke = 0.35) +
  scale_colour_distiller(palette = "YlGnBu", direction = 1, limits = c(vlo, vhi), name = "Shoots\n0.04 m²") +
  scale_fill_distiller(palette = "YlGnBu", direction = 1, limits = c(vlo, vhi), guide = "none") +
  scale_shape_manual(values = c("Observed site" = 21), name = NULL) +
  guides(shape = guide_legend(order = 1, override.aes = list(fill = "grey75", size = 2.4)),
         colour = guide_colourbar(order = 2, barheight = unit(45, "pt"), barwidth = unit(5, "pt"), title.position = "top")) +
  theme(legend.position = c(0.84, 0.78), legend.title = element_text(size = 6.5), legend.text = element_text(size = 6),
        legend.background = element_rect(fill = alpha("white", 0.92), colour = NA), legend.spacing.y = unit(2, "pt"),
        legend.key.size = unit(8, "pt"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
grid$aoa <- factor(ifelse(grid$inside, "Inside AOA", "Outside AOA\n(extrapolation)"), levels = c("Inside AOA", "Outside AOA\n(extrapolation)"))
p7c <- ggplot() + base_map() +
  geom_point(data = grid, aes(longitude, latitude, colour = aoa), shape = 15, size = sq) +
  geom_point(data = site, aes(longitude, latitude), size = 0.9, colour = "black") +
  scale_colour_manual(values = c("Inside AOA" = OI[["skyblue"]], "Outside AOA\n(extrapolation)" = OI[["vermilion"]]), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = c(0.97, 0.97), legend.justification = c(1, 1), legend.text = element_text(size = 6), legend.background = element_rect(fill = alpha("white", 0.92), colour = NA),
        legend.key.size = unit(8, "pt"), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
fig7 <- (p7a + p7b + p7c) + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
save_png(fig7, "main/Figure7.png", 7.6, 7.4)

# ------------------------------------------------------------------ Figure 8
hv <- read.csv(O("hovmoller.csv")); traj <- read.csv(O("traj.csv")); dp <- read.csv(O("depth_profile.csv")); meta <- fromJSON(O("meta.json"))
years <- sort(unique(hv$year)); thr <- meta$temporal_DI_thresh; xt <- seq(min(years), max(years), 3); vmax <- quantile(agg$dens, 0.98, na.rm = TRUE)
lat_rng <- range(hv$lat); traj$region <- factor(traj$region, levels = REG); agg$region <- factor(agg$region, levels = REG)
pa8 <- ggplot(traj, aes(year, colour = region)) + geom_point(aes(y = obs), size = 1.3, alpha = 0.8) + geom_line(aes(y = pred), linewidth = 0.8) +
  scale_colour_manual(values = REG_COL, labels = REG_L) + scale_x_continuous(breaks = xt, limits = c(min(years) - 0.5, max(years) + 0.5)) +
  labs(x = NULL, y = "Mean density 0.04 m²", colour = NULL) + theme(legend.position = "bottom", legend.text = element_text(size = 7.5))
reg_lat <- aggregate(latitude ~ region, site, mean)
depth_panel <- function() ggplot(dp) + geom_ribbon(aes(y = lat, xmin = dep_lo, xmax = dep_hi), fill = "#9a9a9a", alpha = 0.3) +
  geom_path(aes(dep_mean, lat), colour = "#2b2b2b", linewidth = 0.5) +
  scale_x_reverse(breaks = c(0, 5, 10, 15), limits = c(ceiling(max(dp$dep_hi)) + 1, 0), expand = c(0, 0)) +
  scale_y_continuous(breaks = reg_lat$latitude, labels = REG_L[as.character(reg_lat$region)], limits = lat_rng, expand = c(0, 0)) +
  labs(x = "Depth (m)", y = NULL) + theme(axis.text = element_text(size = 7), panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3))
right_lat <- function() scale_y_continuous(limits = lat_rng, expand = c(0, 0), breaks = seq(-33, -30, 1), labels = function(v) sprintf("%.0f°S", abs(v)), position = "right")
pb8 <- ggplot() + geom_raster(data = hv, aes(year, lat, fill = pred), interpolate = TRUE) +
  geom_point(data = agg, aes(year, latitude, fill = dens), shape = 21, size = 1.6, stroke = 0.25) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1, limits = c(0, vmax), oob = scales::squish, name = "Shoots 0.04 m²") +
  scale_x_continuous(breaks = xt, expand = c(0, 0)) + right_lat() + labs(x = NULL, y = NULL) +
  theme(legend.position = "right", legend.title = element_text(size = 8), legend.text = element_text(size = 7), axis.text = element_text(size = 7))
pc8 <- ggplot() + geom_raster(data = hv, aes(year, lat, fill = DI), interpolate = TRUE) +
  geom_contour(data = hv, aes(year, lat, z = DI), breaks = thr, colour = "black", linewidth = 0.5) +
  geom_point(data = agg, aes(year, latitude), size = 0.6, colour = "black", alpha = 0.5) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, limits = c(0, quantile(hv$DI, 0.98, na.rm = TRUE)), oob = scales::squish, name = "Dissimilarity\nindex") +
  scale_x_continuous(breaks = xt, expand = c(0, 0)) + right_lat() + labs(x = "Year", y = NULL) +
  theme(legend.position = "right", legend.title = element_text(size = 8), legend.text = element_text(size = 7), axis.text = element_text(size = 7))
design8 <- "AAAA\nBCCC\nDEEE"
fig8 <- wrap_plots(A = pa8, B = depth_panel(), C = pb8, D = depth_panel(), E = pc8, design = design8, heights = c(0.62, 1, 1)) +
  plot_annotation(tag_levels = list(c("(a)", "", "(b)", "", "(c)")))
save_png(fig8, "main/Figure8.png", 9.4, 10.6)

# ------------------------------------------------------------------ Figure S5: latitudinal transect
g2 <- grid[order(grid$latitude), ]; g2$bin <- cut(g2$latitude, 80)
band <- aggregate(cbind(lat = latitude, pred, lo, hi, inside) ~ bin, g2, mean)
out <- band[band$inside < 0.5, ]
figS5 <- ggplot() +
  geom_rect(data = out, aes(xmin = lat - 0.022, xmax = lat + 0.022, ymin = -Inf, ymax = Inf, fill = "outside area of applicability")) +
  geom_ribbon(data = band, aes(lat, ymin = lo, ymax = hi, fill = "80% quantile forest band")) +
  geom_line(data = band, aes(lat, pred, linetype = "stRF prediction (corridor mean)"), colour = OI[["vermilion"]], linewidth = 0.8) +
  geom_point(data = site, aes(latitude, dens, colour = region), size = 1.9, shape = 21, fill = NA, stroke = 0) +
  geom_point(data = site, aes(latitude, dens, fill = region), size = 1.9, shape = 21, stroke = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c("outside area of applicability" = "grey85", "80% quantile forest band" = alpha(OI[["vermilion"]], 0.18), REG_COL),
                    breaks = c("80% quantile forest band", "outside area of applicability"), name = NULL) +
  scale_colour_manual(values = REG_COL, labels = REG_L, name = NULL) + scale_linetype_manual(values = "solid", name = NULL) +
  guides(colour = guide_legend(override.aes = list(fill = REG_COL, stroke = 0.3, size = 2.2))) +
  scale_x_continuous(breaks = latbr, labels = sprintf("%.1f", abs(latbr))) +
  labs(x = "Latitude (°S)", y = "Shoot density 0.04 m²") +
  theme(legend.text = element_text(size = 7.5), legend.box = "vertical", legend.spacing.y = unit(0, "pt"))
save_png(figS5, "si/FigureS5.png", 8.6, 4.6)

# ------------------------------------------------------------------ Figure S6: forecast skill
ro <- read.csv(O("rolling_origin.csv")); FC <- c(stRF = OI[["vermilion"]], persistence = OI[["blue"]], climatology = "#666666")
lv <- aggregate(rho_level ~ horizon + method, ro, mean); lv$method <- factor(lv$method, levels = names(FC))
an <- aggregate(rho_anomaly ~ horizon, ro[ro$method == "stRF", ], mean)
pa <- ggplot(lv, aes(horizon, rho_level, colour = method)) + geom_line(linewidth = 0.8) + geom_point(size = 2) + scale_colour_manual(values = FC) +
  scale_x_continuous(breaks = unique(lv$horizon)) + labs(x = "Forecast horizon (years ahead)", y = "Spearman (predicted vs observed)", colour = NULL)
pb <- ggplot(an, aes(horizon, rho_anomaly)) + geom_hline(yintercept = 0, linetype = "22", colour = "grey60") +
  geom_line(colour = OI[["vermilion"]], linewidth = 0.8) + geom_point(colour = OI[["vermilion"]], size = 2) + ylim(-0.35, 0.35) +
  scale_x_continuous(breaks = unique(an$horizon)) + labs(x = "Forecast horizon (years ahead)", y = "Corr(predicted, observed) anomaly")
figS6 <- (pa + pb) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")
save_png(figS6, "si/FigureS6.png", 9, 4.2)

# ------------------------------------------------------------------ Figure S7: interval coverage and width
iv <- read.csv(O("intervals.csv")); s <- aggregate(cbind(coverage, width) ~ regime + interval, iv, mean)
IC <- c("ensemble spread" = OI[["grey"]], "QRF" = OI[["blue"]], "split conformal" = OI[["green"]])
s$interval <- factor(s$interval, levels = names(IC))
s$regime <- factor(s$regime, levels = c("interpolation (random 5-fold)", "extrapolation (LORO)"), labels = c("interpolation\n(random 5-fold)", "extrapolation\n(leave one region out)"))
pa <- ggplot(s, aes(regime, coverage, fill = interval)) + geom_col(position = position_dodge(width = 0.8), width = 0.78) +
  geom_hline(aes(yintercept = 0.8, linetype = "nominal 80%"), linewidth = 0.5) + scale_fill_manual(values = IC) + scale_linetype_manual(values = "22") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + labs(x = NULL, y = "Empirical coverage of nominal 80% interval", fill = NULL, linetype = NULL)
pb <- ggplot(s, aes(regime, width, fill = interval)) + geom_col(position = position_dodge(width = 0.8), width = 0.78) + scale_fill_manual(values = IC) +
  labs(x = NULL, y = "Mean interval width (shoots 0.04 m²)", fill = NULL)
figS7 <- (pa + pb) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")
save_png(figS7, "si/FigureS7.png", 9, 4.2)
