# cFLIP-KO bulk RNA-seq: KEGG GSEA highlighting the p38-MAPK pathway
# ----------------------------------------------------------------------------
# Effector cells (NK or T) are co-cultured with target cancer cells at increasing
# ratios; for each treated sample, log2 fold-change vs the untreated control is
# ranked and KEGG GSEA is run (clusterProfiler::gseKEGG). The p38-MAPK pathway
# (hsa04010) is extracted and visualised.
#
# Input: FPKM expression tables (see data/README.md):
#   - NK_RNAseq/gene_expression.xls         (A549 vs NK effector ratios)
#   - T_RNAseq/fpkm_merged_count.annot.tsv  (effector-T RNA-seq)
# ----------------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)   # gseaplot2, dotplot
library(pathview)
library(ggplot2)

DATA_DIR <- "../data"
P38_KEGG_ID <- "hsa04010"   # MAPK signalling pathway (contains p38)

# --- Generic GSEA over one treated-vs-control contrast ----------------------
# expr_mat   : matrix with gene symbols as rownames
# control_col: column name of the untreated control
# treated_cols: vector of treated-sample column names
# outdir     : where to write per-sample results
run_gsea_contrasts <- function(expr_mat, control_col, treated_cols, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  control_expr <- expr_mat[, control_col]

  for (sample_id in treated_cols) {
    message("Processing sample: ", sample_id)

    # Ranked gene list by log2 fold-change vs control
    logFC <- log2(expr_mat[, sample_id] + 1) - log2(control_expr + 1)
    names(logFC) <- rownames(expr_mat)
    geneList <- sort(logFC[!is.na(names(logFC))], decreasing = TRUE)

    # SYMBOL -> ENTREZID
    gene_df <- bitr(names(geneList), fromType = "SYMBOL",
                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    geneList_use <- geneList[gene_df$SYMBOL]
    names(geneList_use) <- gene_df$ENTREZID

    gsea_res <- gseKEGG(geneList = geneList_use, organism = "hsa",
                        nPerm = 1000, minGSSize = 10,
                        pvalueCutoff = 0.5, verbose = FALSE)

    tag <- gsub("[^a-zA-Z0-9]", "_", sample_id)
    sample_outdir <- file.path(outdir, tag)
    dir.create(sample_outdir, showWarnings = FALSE)

    write.csv(gsea_res@result, file.path(sample_outdir, "kegg_gsea.csv"), row.names = FALSE)
    write.csv(subset(gsea_res@result, ID == P38_KEGG_ID),
              file.path(sample_outdir, "p38_kegg_gsea.csv"), row.names = FALSE)

    # Dot plot of up/down-regulated pathways
    if (nrow(gsea_res@result) > 0) {
      gsea_res@result$.sign <- ifelse(gsea_res@result$NES > 0, "upregulated", "downregulated")
      p_dot <- dotplot(gsea_res, showCategory = 20, split = ".sign") +
        facet_grid(. ~ .sign) + ggtitle(paste0(sample_id, " - GSEA KEGG"))
      ggsave(file.path(sample_outdir, "gsea_dotplot.pdf"), p_dot, width = 11, height = 20)
    }

    # p38-MAPK enrichment curve + KEGG pathway view
    if (P38_KEGG_ID %in% gsea_res@result$ID) {
      p_p38 <- gseaplot2(gsea_res, geneSetID = P38_KEGG_ID,
                         title = paste0(sample_id, " - p38 MAPK signalling"))
      ggsave(file.path(sample_outdir, "gsea_p38_curve.pdf"), p_p38, width = 6, height = 5)

      pathview(gene.data = geneList_use, pathway.id = P38_KEGG_ID, species = "hsa",
               out.suffix = paste0("p38_", tag), kegg.native = TRUE,
               same.layer = TRUE, out.dir = sample_outdir)
    }
    message("Done: ", sample_id)
  }
}

# --- NK effector RNA-seq (A549 target) --------------------------------------
nk_raw <- read.delim(file.path(DATA_DIR, "NK_RNAseq", "gene_expression.xls"), header = TRUE)
nk_mat <- as.matrix(nk_raw[, c("fpkm_A549_1to1", "fpkm_A549_1to2",
                               "fpkm_A549_1to4", "fpkm_A549UT")])
rownames(nk_mat) <- nk_raw$gene_symbol
run_gsea_contrasts(nk_mat, control_col = "fpkm_A549UT",
                   treated_cols = c("fpkm_A549_1to1", "fpkm_A549_1to2", "fpkm_A549_1to4"),
                   outdir = "NK_p38_results_GSEA_only")

# --- Effector-T RNA-seq -----------------------------------------------------
t_raw <- read.delim(file.path(DATA_DIR, "T_RNAseq", "fpkm_merged_count.annot.tsv"), header = TRUE)
t_mat <- as.matrix(t_raw[, c("gk1_1cd133", "gk1_2CD133", "gk1_4CD133", "UTCD133")])
rownames(t_mat) <- t_raw$gene_name
run_gsea_contrasts(t_mat, control_col = "UTCD133",
                   treated_cols = c("gk1_1cd133", "gk1_2CD133", "gk1_4CD133"),
                   outdir = "T_p38_results_GSEA_only")
