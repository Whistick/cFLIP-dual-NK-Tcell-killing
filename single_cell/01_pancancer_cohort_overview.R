# Pan-cancer single-cell cohort overview
# ----------------------------------------------------------------------------
# Composition of the merged pan-cancer single-cell atlas: number of samples per
# sequencing cohort, per cancer type, and per tissue type, shown as donut charts.
#
# Input:  a metadata table exported from the processed AnnData/Seurat object
#         (one row per cell) with columns: Sample, Cohort, CancerType_final,
#         Tissue_final. Export it from Python with:
#             adata.obs.to_csv("metadata.csv")
#
# Output: donut charts (rendered to the active graphics device / saved by the user).
# ----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggrepel)

source("../R/plot_themes.R")   # col_cohort

# --- config -----------------------------------------------------------------
DATA_DIR <- "../data"
meta <- read.csv(file.path(DATA_DIR, "metadata.csv"), row.names = 1)

# Collapse to one row per sample (cohort composition is sample-level, not cell-level)
sample_stats <- meta %>% distinct(Sample, .keep_all = TRUE)
n_total <- nrow(sample_stats)
cat("Total number of samples:", n_total, "\n")

# Helper: count + percentage for a grouping column
count_perc <- function(df, col) {
  df %>%
    dplyr::count(.data[[col]]) %>%
    dplyr::rename(group = 1) %>%
    dplyr::mutate(perc = n / sum(n) * 100)
}

cohort_stats <- count_perc(sample_stats, "Cohort")
cancer_stats <- count_perc(sample_stats, "CancerType_final")
tissue_stats <- count_perc(sample_stats, "Tissue_final")

# --- Donut chart with repelled labels ---------------------------------------
# pie_donut(): functional wrapper that draws a donut with outside labels.
#   stats : data frame from count_perc() (columns: group, n, perc)
#   pal   : optional named colour vector
pie_donut <- function(stats, title, pal = NULL) {
  stats <- stats %>%
    arrange(desc(group)) %>%
    mutate(
      ymax       = cumsum(perc),
      ymin       = lag(ymax, default = 0),
      label_pos  = (ymax + ymin) / 2,
      label_text = paste0(group, "\n", round(perc, 1), "%")
    )
  p <- ggplot(stats, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = group)) +
    geom_rect(color = "white") +
    coord_polar(theta = "y") +
    theme_void(base_size = 14) +
    labs(title = title, subtitle = paste0("Total samples: ", n_total)) +
    geom_label_repel(
      aes(x = 4.5, y = label_pos, label = label_text),
      nudge_x = 0.5, segment.size = 0.5, segment.color = "grey50",
      box.padding = 0.3, point.padding = 0.5, size = 3.5, show.legend = FALSE
    )
  if (!is.null(pal)) p <- p + scale_fill_manual(values = pal)
  p
}

pie_donut(cohort_stats, "Samples per sequencing cohort", col_cohort)
pie_donut(cancer_stats, "Samples per cancer type")
pie_donut(tissue_stats, "Samples per tissue type")
