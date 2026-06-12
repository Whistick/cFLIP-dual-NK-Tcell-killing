# Data

No raw or processed data are stored in this repository. All datasets are public
(or will be deposited on acceptance). Download the accessions below and place the
files under `data/` (or point the `DATA_DIR` variable at the top of each script to
wherever you keep them).

## Single-cell RNA-seq — pan-cancer cohorts

Used by `single_cell/01`–`single_cell/03`. Raw data were obtained from the
following repositories and processed (QC, Harmony integration, clustering,
annotation) into a single merged object as described in the manuscript Methods.

| Repository | Accessions |
|------------|------------|
| GEO        | GSE145370, GSE152048, GSE191301, GSE208653, GSE212966, GSE217845 |
| GSA-human  | HRA000051, HRA000321, HRA000403, HRA000963, HRA001006, HRA001130, HRA001748, HRA004556, HRA004767, HRA004971 |
| PRJCA      | PRJCA003766, PRJCA007744 |
| EGA        | EGAD00001005054 |
| ENA        | PRJNA705464, PRJNA510251 |

The processed object expected by the scripts is a Seurat object (R) / AnnData
(Python) carrying at least these metadata columns:
`Sample`, `Cohort`, `CancerType_final`, `Tissue_final`,
`major_celltypes_tier2` (fine cell-type labels, with a `Tumor` label).

## Single-cell RNA-seq — CRC, treatment dynamics (Pre / On / Post)

Used by `single_cell/04`–`single_cell/05`. Colorectal cancer ICB cohort with paired
pre-/on-/post-treatment sampling. The object carries:
`Patient`, `Sample`, `Tissue`, `CellType` (`Epi` = tumour), `Stage` (Pre/On/Post),
`Pathologic_Response` (NR/R), `Response` (SD/PR/CR).

> Accession: _to be added (corresponding-author cohort)._

## Bulk RNA-seq

| Dataset | Accession | Used by | Notes |
|---------|-----------|---------|-------|
| Breast cancer pre/post chemotherapy | GSE191127 | `bulk_rnaseq/01` | Paired before/after, raw read counts |
| cFLIP-KO A549 / SK-Hep1 (NK & T effector) | _to be added_ | `bulk_rnaseq/02` | Authors' RNA-seq, FPKM matrices |

File names referenced by the scripts (rename or repath as needed):

- `GSE191127_readcounts_prepost.txt`  — breast pre/post counts
- `NK_RNAseq/gene_expression.xls`      — cFLIP-KO A549, FPKM
- `T_RNAseq/fpkm_merged_count.annot.tsv` — effector-T RNA-seq, FPKM
