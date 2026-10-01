# Figures 2, 4 and 5: feature block Shapley decomposition from results/expanded_T*.csv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figstyle.R"))
t1 <- read.csv(res("expanded_T1_decomposition_by_design.csv")); t2 <- read.csv(res("expanded_T2_temporal_by_T.csv"))
t3 <- read.csv(res("expanded_T3_detection.csv")); t4 <- read.csv(res("expanded_T4_abundance.csv")); t5 <- read.csv(res("expanded_T5_NN0_by_design.csv"))
COMP <- c("Algorithmic flexibility", "Coordinates", "Distance fields", "Spatial lags", "Temporal")
short <- function(x) factor(ifelse(x == "Algorithmic flexibility", "Algorithmic", x), levels = names(COMP_COL))
pd <- position_dodge(width = 0.8)

bars <- function(tab, groups, labels, ylab) {
  d <- tab[tab$component %in% COMP & tab$subset %in% groups, ]
  d$group <- factor(d$subset, levels = groups, labels = labels); d$component <- short(d$component)
  ggplot(d, aes(group, value, fill = component)) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
    geom_col(position = pd, width = 0.8) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), position = pd, width = 0.25, linewidth = 0.35) +
    scale_fill_manual(values = COMP_COL, breaks = names(COMP_COL)) +
    labs(x = NULL, y = ylab, fill = NULL)
}

# ---- Figure 2: (a) by design, (b) NN0 control by design
pa <- bars(t1, c("pooled (all designs)", "random", "stratified", "clustered"), c("pooled", "random", "stratified", "clustered"),
           "Shapley contribution (ΔAUC)")
d5 <- t5[t5$quantity %in% c("Temporal marginal, no control", "Temporal marginal, NN0 controlled"), ]
d5$design <- factor(d5$design, levels = c("random", "stratified", "clustered"))
d5$control <- factor(ifelse(grepl("NN0", d5$quantity), "temporal, NN0 controlled (b)", "temporal, no control (b)"),
                     levels = c("temporal, no control (b)", "temporal, NN0 controlled (b)"))
pb <- ggplot(d5, aes(design, value, fill = control)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.68, colour = COMP_COL[["Temporal"]], linewidth = 0.4) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.35) +
  scale_fill_manual(values = c("temporal, no control (b)" = COMP_COL[["Temporal"]], "temporal, NN0 controlled (b)" = "white")) +
  labs(x = NULL, y = "Temporal increment when added last (ΔAUC)", fill = NULL)
fig2 <- (pa + pb) + plot_layout(widths = c(1.35, 1), guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")
save_png(fig2, "main/Figure2.png", 10, 4.2)

# ---- Figure 4: temporal contribution by T
t2$quantity <- factor(t2$quantity, levels = c("Shapley", "Added last", "Added last, NN0 controlled"))
fig4 <- ggplot(t2, aes(n_time, value, shape = quantity, linetype = quantity)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.25, linewidth = 0.4, colour = COMP_COL[["Temporal"]], linetype = "solid") +
  geom_line(colour = COMP_COL[["Temporal"]], linewidth = 0.7) +
  geom_point(colour = COMP_COL[["Temporal"]], fill = COMP_COL[["Temporal"]], size = 2.2) +
  annotate("text", x = 1.2, y = 0.0055, label = "placebo", size = 2.6, colour = "grey35", hjust = 0) +
  scale_shape_manual(values = c(16, 15, 17)) + scale_linetype_manual(values = c("solid", "22", "11")) +
  scale_x_continuous(breaks = c(1, 3, 5, 10)) +
  labs(x = "Number of timesteps", y = "Temporal contribution (ΔAUC)", shape = NULL, linetype = NULL)
save_png(fig4, "main/Figure4.png", 6.4, 3.9)

# ---- Figure 5: (a) detection, (b) occurrence vs abundance shares
gd <- unique(t3$subset); gd <- gd[order(-as.numeric(sub(".* ", "", gd)))]
pa5 <- bars(t3, gd, sprintf("detection %.1f", as.numeric(sub(".* ", "", gd))), "Shapley contribution (ΔAUC)")
bp <- t1[t1$subset == "pooled (all designs)" & t1$component %in% COMP, ]; ab <- t4[t4$subset == "abundance pooled" & t4$component %in% COMP, ]
sh <- rbind(data.frame(component = short(bp$component), share = bp$share_pct, response = "occurrence, AUC (b)"),
            data.frame(component = short(ab$component), share = ab$share_pct, response = "abundance, Spearman ρ (b)"))
sh$response <- factor(sh$response, levels = c("occurrence, AUC (b)", "abundance, Spearman ρ (b)"))
pb5 <- ggplot(sh, aes(component, share, fill = component, alpha = response)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.7, colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = COMP_COL, guide = "none") +
  scale_alpha_manual(values = c(1, 0.4)) +
  guides(alpha = guide_legend(override.aes = list(fill = "grey30"))) +
  labs(x = NULL, y = "Share of total gain (%)", alpha = NULL) +
  theme(axis.text.x = element_text(size = 7.5, angle = 12, hjust = 0.8))
fig5 <- (pa5 + pb5) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")
save_png(fig5, "main/Figure5.png", 10, 4.2)
