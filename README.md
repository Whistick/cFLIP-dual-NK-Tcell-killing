# cFLIP / Extrinsic Apoptosis in Dual NK and T cell-mediated Killing

Code to reproduce the bioinformatic and clinical-genomics analyses for the manuscript:

> **Extrinsic apoptosis activation unleashes vulnerability enhancing dual NK and T cell-mediated killing in cancer**
> Zhen Chen, Gabriele Casagrande Raffi, Ziji Zhao, Chaopin Yang, et al.

This repository covers the **computational analyses only** (public single-cell and
bulk RNA-seq cohorts). Wet-lab assays (CRISPR screen, cytotoxicity, flow cytometry,
in vivo experiments) are described in the manuscript Methods and do not require code.

The central computational finding is that the apoptosis regulator **cFLIP (gene `CFLAR`)**
is negatively associated with intratumoral NK/T-cell infiltration and immunotherapy (ICB)
response, and that immune-driven inflammation acts through the **p38-MAPK** pathway.

## Repository layout

Scripts are organised **by data modality**, not by figure number.

```
.
├── single_cell/        # Pan-cancer single-cell RNA-seq cohort analyses
│   ├── 01_pancancer_cohort_overview.R          # Cohort / cancer-type / tissue composition
│   ├── 02_pancancer_umap_export.py             # Export UMAP coords + labels from AnnData
│   ├── 02_pancancer_umap_plot.R                # Rasterised UMAP of major cell types
│   ├── 03_cflar_immune_infiltration.py         # CFLAR-high/low tumours vs NK/T infiltration
│   ├── 03_cflar_immune_infiltration_plot.R     # Box plots per cancer type
│   ├── 04_crc_treatment_dynamics_cell_level.R  # CRC CFLAR dynamics (Pre/On/Post, R vs NR)
│   └── 05_crc_treatment_dynamics_paired.R      # CRC patient-paired CFLAR+ tumour fraction
│
├── bulk_rnaseq/        # Bulk RNA-seq cohorts
│   ├── 01_breast_prepost_chemo.R               # GSE191127 pre/post chemo, paired DESeq2
│   └── 02_cflip_ko_gsea_p38.R                  # cFLIP-KO A549/SK-Hep1, KEGG GSEA (p38-MAPK)
│
├── R/
│   └── plot_themes.R   # Shared colour palettes used across the R scripts
│
└── data/
    └── README.md       # Public data accessions + download instructions
```

## Input data

Raw and processed data are **not** stored here (single-cell objects are tens of GB).
See [`data/README.md`](data/README.md) for every public accession (GEO / GSA / EGA / ENA)
and how to download or reconstruct the intermediate objects each script expects.

Several scripts assume a pre-processed single-cell object is already loaded in the
session as `adata` (a Seurat object in R, or an AnnData in Python). The pre-processing
(QC, Harmony integration, clustering, annotation) follows the manuscript Methods; the
expected metadata columns are documented at the top of each script.

## Environment

- **Python**: see [`requirements.txt`](requirements.txt) (scanpy / anndata / pandas / scipy).
- **R 4.3+**: see [`r_session_info.md`](r_session_info.md). Key packages: Seurat (v4.4.0),
  Harmony (v1.2.0), SCP, ggplot2, ggpubr, DESeq2, clusterProfiler, ggrastr.

## Citation

Please cite the manuscript above if you use this code.
