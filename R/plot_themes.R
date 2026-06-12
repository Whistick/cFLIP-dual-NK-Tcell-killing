# Shared colour palettes used across the R figure scripts.
# Source this file from a script with:  source("../R/plot_themes.R")

# --- Treatment response / stage ---------------------------------------------
col_response   <- c(NR = "#EB5F61", R = "#63A1CB")          # pathologic response
col_response3  <- c(SD = "#EB5F61", PR = "#63A1CB", CR = "#47A364")
col_stage      <- c(Pre = "#FF8820", On = "#3182BA", Post = "#47A364")

# --- MHC-I / CFLAR high vs low ----------------------------------------------
col_mhc_group  <- c(MHC_I_high = "#DE586E", MHC_I_low = "#1B83BB")

# --- Cell types (spatial HCC) -----------------------------------------------
col_celltype_hcc <- c(
  Tumor   = "#A6CEE3",
  Stromal = "#1F78B4",
  T_cell  = "#B2DF8A",
  Plasma  = "#33A02C",
  NK_cell = "#FDBF6F",
  Myeloid = "#FF7F00",
  BECs    = "#FB9A99"
)

# --- Cell types (pan-cancer UMAP) -------------------------------------------
col_celltype_pancancer <- c(
  "T"           = "#FDBF6F",
  "Tumor"       = "#A6CEE3",
  "Macrophage"  = "#B2DF8A",
  "EC"          = "#FB9A99",
  "Mesenchymal" = "#1F78B4",
  "B"           = "#33A02C",
  "Neutrophil"  = "#FF7F00",
  "NK"          = "#FF7F00",
  "Monocyte"    = "#CAB2D6",
  "DC"          = "#6A3D9A",
  "Mast"        = "#B15928",
  "ILC"         = "#E31A1C",
  "N"           = "#FDBF6F",
  "GMP"         = "#B2DF8A",
  "Mesothelial" = "#999999"
)

# --- Sequencing cohorts (single-cell) ---------------------------------------
col_cohort <- c(
  "GSE145370" = "#1B9E77", "GSE152048" = "#D95F02", "GSE191301" = "#7570B3",
  "GSE208653" = "#E7298A", "GSE212966" = "#66A61E", "GSE217845" = "#E6AB02",
  "HRA000051" = "#A6761D", "HRA000321" = "#B3DE69", "HRA000963" = "#8DD3C7",
  "HRA001006" = "#FFFFB3", "HRA004556" = "#BEBADA", "HRA004767" = "#FB8072",
  "HRA004971" = "#80B1D3", "PRJCA003766" = "#FB9A99", "PRJCA007744" = "#FDB462",
  "PRJNA705464" = "#FCCDE5"
)
