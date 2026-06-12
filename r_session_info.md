# R environment

Analyses were run on **R 4.3**. Install the packages below from CRAN / Bioconductor.

## CRAN

```r
install.packages(c(
  "Seurat",      # v4.4.0 (single-cell)
  "harmony",     # v1.2.0 (batch integration)
  "ggplot2", "ggpubr", "ggrepel", "ggrastr", "ggthemes",
  "dplyr", "reshape2", "tibble", "viridisLite"
))
```

## Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "DESeq2",            # bulk differential expression
  "limma", "edgeR",    # bulk RNA-seq utilities
  "clusterProfiler",   # KEGG GSEA
  "org.Hs.eg.db",      # human gene annotation
  "pathview"           # KEGG pathway rendering
))
```

## GitHub

```r
# SCP — single-cell plotting (CellDimPlot)
remotes::install_github("zhanghao-njmu/SCP")
```

`ggthemes::theme_few()` is used as the base theme across the figure scripts.
