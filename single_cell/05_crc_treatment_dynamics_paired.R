# CRC ICB cohort: patient-paired CFLAR+ tumour-cell fraction across treatment
# ----------------------------------------------------------------------------
# For each sample, compute the fraction of CFLAR-positive tumour cells, then track
# how it changes across treatment stages (Pre/On/Post) within each patient, with
# same-patient time points connected by dashed lines.
#
# Input:  `adata`, a Seurat object in the session, with metadata columns
#         Tissue (exclude "B" = blood), CellType ("Epi" = tumour), Patient, Sample,
#         Stage (Pre/On/Post), Pathologic_Response (NR/R), Response (SD/PR/CR).
# ----------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

source("../R/plot_themes.R")   # col_response, col_response3

# 1. Drop blood samples and fix factor levels
adata <- subset(adata, subset = Tissue != "B")
adata$Stage <- factor(adata$Stage, levels = c("Pre", "On", "Post"))
adata$Pathologic_Response <- factor(adata$Pathologic_Response, levels = c("NR", "R"))
adata$Response <- factor(adata$Response, levels = c("SD", "PR", "CR"))

# 2. Per-sample fraction of CFLAR+ tumour cells (CFLAR > 0 counts as positive)
epi <- subset(adata, subset = CellType == "Epi")
epi$CFLAR_pos <- GetAssayData(epi, slot = "data")["CFLAR", ] > 0

ratio_df <- epi@meta.data %>%
  group_by(Patient, Sample, Stage, Pathologic_Response, Response) %>%
  summarise(total = n(), CFLAR_pos = sum(CFLAR_pos),
            ratio = CFLAR_pos / total, .groups = "drop")

# 3. Keep only patients sampled at more than one stage (for paired comparison)
paired_patients <- ratio_df %>%
  group_by(Patient) %>%
  filter(length(unique(Stage)) > 1) %>%
  pull(Patient) %>% unique()
ratio_df_paired <- ratio_df %>% filter(Patient %in% paired_patients)

# Optional paired t-test across stages
pair_test <- pairwise.t.test(ratio_df_paired$ratio, ratio_df_paired$Stage,
                             paired = TRUE, p.adjust.method = "none")
print(pair_test$p.value)

# 4a. Box plot coloured by pathologic response (NR/R)
ggplot(ratio_df_paired, aes(x = Stage, y = ratio)) +
  geom_boxplot(fill = "grey80", color = "grey50", width = 0.5, outlier.shape = NA) +
  geom_point(aes(color = Pathologic_Response, group = Patient),
             size = 2, position = position_jitter(width = 0)) +
  geom_line(aes(group = Patient), color = "grey70", linetype = "dashed", alpha = 0.6) +
  scale_color_manual(values = col_response) +
  labs(x = "Treatment stage", y = "Fraction of CFLAR+ tumour cells",
       title = "CFLAR+ tumour-cell fraction across treatment", color = "Response") +
  theme_classic() +
  stat_compare_means(aes(label = ..p.format..),
                     comparisons = list(c("Pre", "On"), c("On", "Post")),
                     method = "wilcox.test", paired = FALSE, size = 4) +
  theme(legend.position = "right", plot.title = element_text(face = "bold", size = 16))

# 4b. Same plot coloured by RECIST response (SD/PR/CR)
ggplot(ratio_df_paired, aes(x = Stage, y = ratio)) +
  geom_boxplot(fill = "grey80", color = "grey50", width = 0.5, outlier.shape = NA) +
  geom_point(aes(color = Response, group = Patient),
             size = 2, position = position_jitter(width = 0)) +
  geom_line(aes(group = Patient), color = "grey70", linetype = "dashed", alpha = 0.6) +
  scale_color_manual(values = col_response3) +
  labs(x = "Treatment stage", y = "Fraction of CFLAR+ tumour cells",
       title = "CFLAR+ tumour-cell fraction across treatment", color = "Response") +
  theme_classic() +
  stat_compare_means(aes(label = ..p.format..),
                     comparisons = list(c("Pre", "On"), c("On", "Post")),
                     method = "wilcox.test", paired = FALSE, size = 4) +
  theme(legend.position = "right", plot.title = element_text(face = "bold", size = 16))
