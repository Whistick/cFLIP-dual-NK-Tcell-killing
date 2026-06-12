# CRC ICB cohort: tumour-cell CFLAR dynamics across treatment (cell level)
# ----------------------------------------------------------------------------
# Single-cell CFLAR expression in tumour (epithelial) cells, compared between
# responders (R) and non-responders (NR) across treatment stages (Pre/On/Post).
#
# Input:  `adata`, a Seurat object already in the session, with metadata columns
#         CellType ("Epi" = tumour), Pathologic_Response (NR/R), Stage (Pre/On/Post).
# ----------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

source("../R/plot_themes.R")   # col_response

# Tumour (epithelial) cells only
tu <- adata[, adata$CellType %in% c("Epi")]
df <- FetchData(tu, vars = c("CFLAR", "Pathologic_Response", "Stage"))
df$Stage <- factor(df$Stage, levels = c("Pre", "On", "Post"))
df$Pathologic_Response <- factor(df$Pathologic_Response, levels = c("NR", "R"))

# Use CFLAR-expressing cells for the distribution plots
df_pos <- df %>% filter(CFLAR > 0)

# --- (1) Violin + box: R vs NR, faceted by stage ----------------------------
ggplot(df_pos, aes(x = Pathologic_Response, y = CFLAR, fill = Pathologic_Response)) +
  geom_violin(position = position_dodge(width = 0.8), alpha = 0.6, trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.8) +
  facet_wrap(~ Stage, scales = "free_y") +
  scale_fill_manual(values = col_response) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.10))) +
  theme_classic() +
  labs(y = "CFLAR expression", x = "Pathologic Response", fill = "Response") +
  stat_compare_means(comparisons = list(c("NR", "R")), method = "wilcox.test",
                     label = "p.format", method.args = list(alternative = "greater"),
                     size = 5) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        strip.text = element_text(size = 13, face = "bold"), legend.position = "none")

# --- (2) Trend of NR-vs-R significance across stages ------------------------
# p-values computed from the Wilcoxon tests above, shown as -log10(p).
trend_df <- data.frame(group = c("Pre", "On", "Post"),
                       pval  = c(0.00034, 1.9e-5, 1.9e-9))
trend_df$group <- factor(trend_df$group, levels = c("Pre", "On", "Post"))
trend_df$neglog10p <- -log10(trend_df$pval)

ggplot(trend_df, aes(x = group, y = neglog10p)) +
  geom_col(fill = "#F4A582", width = 0.6) +
  geom_line(aes(group = 1), color = "#B2182B", linewidth = 1.2) +
  geom_point(size = 3, color = "#B2182B") +
  geom_text(aes(label = round(neglog10p, 2)), vjust = -0.6, size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_classic() +
  labs(x = "", y = expression(-log[10](p)),
       title = "NR vs R: increasing significance over treatment") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(color = "black"), panel.grid.minor = element_blank())

# --- (3) Dot plot: mean expression and percent expressing -------------------
df_summary <- df %>%
  group_by(Stage, Pathologic_Response) %>%
  summarise(mean_expr = mean(CFLAR), pct_expr = mean(CFLAR > 0), .groups = "drop")

ggplot(df_summary, aes(x = Stage, y = Pathologic_Response)) +
  geom_point(aes(size = pct_expr, color = mean_expr)) +
  scale_size_continuous(range = c(2, 10), name = "Pct. expressing") +
  scale_color_viridis_c(name = "Mean expression", option = "C") +
  theme_classic() +
  labs(x = "Treatment Stage", y = "Response Group") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))

# --- (4) Mean +/- SE line plot across stages --------------------------------
df_se <- df %>%
  group_by(Stage, Pathologic_Response) %>%
  summarise(mean_expr = mean(CFLAR, na.rm = TRUE),
            se = sd(CFLAR, na.rm = TRUE) / sqrt(n()), .groups = "drop")

ggplot(df_se, aes(x = Stage, y = mean_expr,
                  group = Pathologic_Response, color = Pathologic_Response)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_expr - se, ymax = mean_expr + se), width = 0.15, size = 0.8) +
  scale_color_manual(values = col_response) +
  theme_classic() +
  labs(y = "Mean CFLAR expression", x = "Treatment Stage", color = "Response") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
