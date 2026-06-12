# CFLAR / MHC-I high vs low: NK/T infiltration box plots
# ----------------------------------------------------------------------------
# Reads plot_df_for_R.csv produced by 03_cflar_immune_infiltration.py and draws
# the per-cancer-type and pooled box plots comparing target-cell infiltration
# between MHC-I (CFLAR) high and low tumours.
# ----------------------------------------------------------------------------

library(ggplot2)
library(ggpubr)
library(ggthemes)
library(dplyr)

source("../R/plot_themes.R")   # col_mhc_group

DATA_DIR <- "../data"
df <- read.csv(file.path(DATA_DIR, "cflar_infiltration", "plot_df_for_R.csv"),
               stringsAsFactors = FALSE)

# Order the two groups and drop cancer types with too few samples to plot
df$MHC_group <- factor(df$MHC_group, levels = c("MHC_I_high", "MHC_I_low"))
drop_cancers <- c("BRCA", "MELA", "OV", "NPC", "PC", "NB", "MM",
                  "THCA", "UCEC", "GBM", "CRC", "CTCL", "LYM")
df <- df[!df$CancerType %in% drop_cancers, ]

# --- Faceted: one panel per cancer type -------------------------------------
ggplot(df, aes(x = MHC_group, y = percentage, fill = MHC_group)) +
  geom_boxplot(aes(color = MHC_group), alpha = 0.6, width = 0.4, outlier.shape = NA) +
  geom_point(aes(color = MHC_group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.4),
             size = 0.6, show.legend = FALSE) +
  facet_wrap(~ CancerType, scales = "free_y") +
  scale_color_manual(values = col_mhc_group) +
  scale_fill_manual(values = c(MHC_I_high = "#FFFFFF", MHC_I_low = "#FFFFFF")) +
  labs(x = "MHC-I group", y = "Percentage of target cells (%)") +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  stat_compare_means(aes(label = ..p.format..),
                     comparisons = list(c("MHC_I_high", "MHC_I_low")),
                     method = "t.test", paired = FALSE, size = 3,
                     method.args = list(alternative = "less"))

# --- Pooled: all cancer types together --------------------------------------
df_pooled <- df[df$tumor_cells >= 30, ]
ggplot(df_pooled, aes(x = MHC_group, y = percentage, fill = MHC_group)) +
  geom_boxplot(aes(color = MHC_group), alpha = 0.6, width = 0.5, outlier.shape = NA) +
  geom_point(aes(color = MHC_group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.5),
             size = 0.6, show.legend = FALSE) +
  scale_color_manual(values = col_mhc_group) +
  scale_fill_manual(values = col_mhc_group) +
  labs(x = "MHC-I group", y = "Percentage of target cells (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  stat_compare_means(aes(label = ..p.format..),
                     comparisons = list(c("MHC_I_high", "MHC_I_low")),
                     method = "t.test", paired = FALSE, size = 3,
                     method.args = list(alternative = "less"))
