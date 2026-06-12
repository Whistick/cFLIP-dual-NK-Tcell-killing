# Breast cancer pre/post chemotherapy bulk RNA-seq (GSE191127)
# ----------------------------------------------------------------------------
# Paired before/after-treatment breast cancer samples. Build DESeq2-normalised
# expression and compare immune marker genes (e.g. NCAM1, CD8A) before vs after,
# with same-patient samples paired.
#
# Input:  data/GSE191127_readcounts_prepost.txt (raw counts + gene annotation)
# ----------------------------------------------------------------------------

library(DESeq2)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(reshape2)

DATA_DIR <- "../data"

# --- 1. Read counts and split annotation from the expression matrix ---------
expr <- read.table(file.path(DATA_DIR, "GSE191127_readcounts_prepost.txt"),
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE, quote = "")

gene_info <- expr[, c("ensembl_gene_id", "external_gene_id",
                      "start_position", "end_position")]
gene_info$length <- gene_info$end_position - gene_info$start_position + 1
rownames(expr) <- expr$ensembl_gene_id

anno_cols <- c("ensembl_gene_id", "gene_biotype", "chromosome_name",
               "start_position", "end_position", "external_gene_id", "description")
expr_mat <- expr[, !(colnames(expr) %in% anno_cols)]

# --- 2. Sample metadata (patient id and Before/After group from sample name) -
sample_names <- colnames(expr_mat)
# Sample name encodes ..._<patient>_<B|S>_...  ; B = Before, S = (after)
patient_id <- sapply(strsplit(sample_names, "_"), `[`, 3)
group <- sub(".*_([BS])_.*", "\\1", sample_names)
sample_info <- data.frame(
  sample  = sample_names,
  patient = patient_id,
  group   = ifelse(group == "B", "Before", "After"),
  stringsAsFactors = FALSE
)

# --- 3. log2(FPKM) (for plotting) -------------------------------------------
gene_length_kb <- gene_info$length / 1000
names(gene_length_kb) <- gene_info$ensembl_gene_id
total_counts <- colSums(expr_mat)
fpkm <- sweep(expr_mat, 1, gene_length_kb[rownames(expr_mat)], FUN = "/")
fpkm <- sweep(fpkm, 2, total_counts / 1e6, FUN = "/")
logfpkm <- log2(fpkm + 1)

# --- 4. DESeq2 (group design) -----------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = expr_mat, colData = sample_info,
                              design = ~ group)
keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep, ]
dds <- DESeq(dds)

# --- 5. Plot selected genes: Before vs After (paired) -----------------------
genes_of_interest <- c("NCAM1", "CD8A")
gene_match <- setNames(gene_info$external_gene_id, gene_info$ensembl_gene_id)
sel_ensembl <- names(gene_match[gene_match %in% genes_of_interest])
if (length(sel_ensembl) == 0) stop("None of the requested genes were found.")

sel_counts <- as.matrix(logfpkm[sel_ensembl, , drop = FALSE])
rownames(sel_counts) <- gene_match[rownames(sel_counts)]

long <- melt(sel_counts, varnames = c("gene", "sample"), value.name = "expression")
long$group   <- factor(sample_info$group[match(long$sample, sample_info$sample)],
                        levels = c("Before", "After"))
long$patient <- sample_info$patient[match(long$sample, sample_info$sample)]

ggplot(long, aes(x = group, y = expression, color = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_line(aes(group = patient), color = "gray70", alpha = 0.6) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.13))) +
  stat_compare_means(aes(label = ..p.format..), method = "t.test",
                     comparisons = list(c("Before", "After")), paired = TRUE, size = 4) +
  theme_few() +
  labs(title = "Treatment Before vs After", x = "Group", y = "log2(FPKM + 1)") +
  scale_color_manual(values = c(Before = "steelblue", After = "firebrick"))
