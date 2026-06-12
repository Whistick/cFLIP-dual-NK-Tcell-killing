# Pan-cancer UMAP of major cell types
# ----------------------------------------------------------------------------
# Draws a rasterised UMAP coloured by major cell type, using the coordinates
# exported by 02_pancancer_umap_export.py.
# ----------------------------------------------------------------------------

library(ggplot2)
library(ggrastr)

source("../R/plot_themes.R")   # col_celltype_pancancer

DATA_DIR <- "../data"
df <- read.csv(file.path(DATA_DIR, "umap_metadata.csv"))

ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = major_celltypes_tier2)) +
  geom_point_rast(size = 0.1, raster.dpi = 330) +   # rasterise the millions of points
  scale_color_manual(values = col_celltype_pancancer, na.value = "grey80") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 10),
    legend.title    = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
