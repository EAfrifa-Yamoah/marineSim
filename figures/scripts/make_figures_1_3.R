# Figures 1 and 3: feature block Shapley decomposition from results/expanded_T*.csv
# Figure 1: (a) along the sample size gradient (T9), (b) by sampling design (T1), (c) NN0 control (T5)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figstyle.R"))
t1 <- read.csv(res("expanded_T1_decomposition_by_design.csv")); t2 <- read.csv(res("expanded_T2_temporal_by_T.csv"))
t3 <- read.csv(res("expanded_T3_detection.csv")); t4 <- read.csv(res("expanded_T4_abundance.csv")); t5 <- read.csv(res("expanded_T5_NN0_by_design.csv"))
t9 <- read.csv(res("expanded_T9_decomposition_by_n.csv"))
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

# ---- Figure 1: (a) by sample size, (b) by design, (c) NN0 control by design
d9 <- t9[t9$component %in% COMP, ]; d9$component <- short(d9$component)
pn <- ggplot(d9, aes(n, value, colour = component, fill = component)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8, shape = 21, fill = "white", stroke = 0.8) +
  scale_x_log10(breaks = c(30, 50, 100, 200, 500)) +
  scale_colour_manual(values = COMP_COL, breaks = names(COMP_COL)) +
  scale_fill_manual(values = COMP_COL, breaks = names(COMP_COL), guide = "none") +
  labs(x = "Sites per timestep (n)", y = "Shapley contribution (ΔAUC)", colour = NULL)
pa <- bars(t1, c("pooled (all designs)", "random", "stratified", "clustered"), c("pooled", "random", "stratified", "clustered"),
           "Shapley contribution (ΔAUC)") + scale_fill_manual(values = COMP_COL, breaks = names(COMP_COL), guide = "none")
d5 <- t5[t5$quantity %in% c("Temporal marginal, no control", "Temporal marginal, NN0 controlled"), ]
d5$design <- factor(d5$design, levels = c("random", "stratified", "clustered"))
d5$control <- factor(ifelse(grepl("NN0", d5$quantity), "temporal, NN0 controlled (c)", "temporal, no control (c)"),
                     levels = c("temporal, no control (c)", "temporal, NN0 controlled (c)"))
pb <- ggplot(d5, aes(design, value, fill = control)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.68, colour = COMP_COL[["Temporal"]], linewidth = 0.4) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.35) +
  scale_fill_manual(values = c("temporal, no control (c)" = COMP_COL[["Temporal"]], "temporal, NN0 controlled (c)" = "white")) +
  labs(x = NULL, y = "Temporal increment when added last (ΔAUC)", fill = NULL)
fig1 <- (pn + pa + pb) + plot_layout(widths = c(1.15, 1.25, 0.9), guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")
save_png(fig1, "main/Figure1.png", 12.5, 4.3)

# ---- Figure 3: (a) detection, (b) occurrence vs abundance shares
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
fig3 <- (pa5 + pb5) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")
save_png(fig3, "main/Figure3.png", 10, 4.2)
