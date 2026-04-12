# ============================================================
# DIABETES MELLITUS GENE EXPRESSION ANALYSIS
# ============================================================
# Author:   Saima Imran
# GitHub:   https://github.com/saima-imran
# Dataset:  GSE15932 (NCBI GEO)
# Disease:  Diabetes Mellitus vs Normal Pancreatic Tissue
# ============================================================
# Description:
# This script analyses gene expression differences between
# diabetic patients and normal controls using public GEO
# microarray data. The analysis focuses on identifying
# key genes dysregulated in diabetes mellitus — directly
# relevant to precision medicine and diagnostics.
#
# Groups:
#   - Cancer + Diabetes (8 samples)
#   - Diabetes Only     (8 samples)
#   - Normal            (8 samples)
#   - Other             (8 samples)
#
# Comparison: Diabetes Only vs Normal
# ============================================================

# ============================================================
# STEP 1: LOAD PACKAGES
# ============================================================

library(GEOquery)
library(limma)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(Biobase)

cat("✓ All packages loaded!\n")

# ============================================================
# STEP 2: CREATE FOLDER STRUCTURE
# ============================================================

dir.create("results",                    showWarnings = FALSE)
dir.create("results/figures",            showWarnings = FALSE)
dir.create("results/figures/diabetes",   showWarnings = FALSE)
dir.create("results/data",               showWarnings = FALSE)

cat("✓ Folders created!\n")

# ============================================================
# STEP 3: DOWNLOAD DATA FROM NCBI GEO
# ============================================================

cat("\n--- Downloading GSE15932 Diabetes Dataset ---\n")
cat("Dataset: Pancreatic islet gene expression in diabetes\n\n")

gse_diabetes <- getGEO("GSE15932", GSEMatrix = TRUE, getGPL = FALSE)

cat("✓ Dataset downloaded!\n")

# ============================================================
# STEP 4: EXTRACT EXPRESSION DATA
# ============================================================

cat("\n--- Extracting Expression Data ---\n")

diabetes_expr <- exprs(gse_diabetes[[1]])
sample_data   <- pData(gse_diabetes[[1]])

cat("Expression matrix:", nrow(diabetes_expr), "genes x",
    ncol(diabetes_expr), "samples\n")

# ============================================================
# STEP 5: EXTRACT PATIENT GROUPS
# ============================================================

cat("\n--- Extracting Patient Groups ---\n")

all_titles <- sample_data$title

groups_diabetes <- case_when(
  grepl("pancreatic cancer and diabetes", all_titles) ~ "Cancer+Diabetes",
  grepl("Patient B", all_titles)                      ~ "Diabetes Only",
  grepl("normal|Normal|control|Control", all_titles)  ~ "Normal",
  TRUE                                                 ~ "Other"
)

cat("Patient groups:\n")
print(table(groups_diabetes))

# ============================================================
# STEP 6: DIFFERENTIAL EXPRESSION — DIABETES vs NORMAL
# ============================================================

cat("\n--- Differential Expression: Diabetes vs Normal ---\n")

# Select Diabetes Only and Normal samples
diabetes_idx <- which(groups_diabetes == "Diabetes Only")
normal_idx   <- which(groups_diabetes == "Normal")

cat("Diabetes samples:", length(diabetes_idx), "\n")
cat("Normal samples:  ", length(normal_idx), "\n")

# Subset expression
diab_subset <- diabetes_expr[, c(diabetes_idx, normal_idx)]

# Create groups
groups <- factor(c(rep("Diabetes", length(diabetes_idx)),
                   rep("Normal",   length(normal_idx))))

# Design matrix
design <- model.matrix(~ 0 + groups)
colnames(design) <- c("Diabetes", "Normal")

# Contrast
contrast <- makeContrasts(Diabetes - Normal, levels = design)

# Fit model with limma
fit <- eBayes(contrasts.fit(lmFit(diab_subset, design), contrast))

# Extract results
diab_results <- topTable(fit, coef = 1, number = Inf, sort.by = "P")

# Apply thresholds
diab_deg <- diab_results %>%
  mutate(
    Gene = rownames(.),
    Significance = case_when(
      P.Value < 0.05 & logFC >  1.5 ~ "Up in Diabetes",
      P.Value < 0.05 & logFC < -1.5 ~ "Down in Diabetes",
      TRUE ~ "Not Significant"
    )
  )

cat("\nDifferential Expression Results:\n")
cat("Total genes tested:  ", nrow(diab_results), "\n")
cat("Up in Diabetes:      ", sum(diab_deg$Significance == "Up in Diabetes"), "\n")
cat("Down in Diabetes:    ", sum(diab_deg$Significance == "Down in Diabetes"), "\n")

# ============================================================
# STEP 7: VOLCANO PLOT
# ============================================================

cat("\n--- Generating Volcano Plot ---\n")

top_label <- diab_deg %>%
  filter(Significance != "Not Significant") %>%
  arrange(P.Value) %>%
  head(15)

png("results/figures/diabetes/01_Diabetes_volcano.png",
    width = 1400, height = 1000, res = 120)

ggplot(diab_deg, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(colour = Significance), alpha = 0.5, size = 1.5) +
  scale_colour_manual(values = c(
    "Up in Diabetes"   = "#E74C3C",
    "Down in Diabetes" = "#3498DB",
    "Not Significant"  = "#AAAAAA")) +
  geom_vline(xintercept = c(-1.5, 1.5),
             linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", linewidth = 0.5) +
  geom_text_repel(
    data = top_label,
    aes(label = Gene),
    size = 3, max.overlaps = 15, colour = "black") +
  labs(
    title    = "Diabetes Mellitus vs Normal — Differential Expression",
    subtitle = "GSE15932 | Threshold: p < 0.05 | |logFC| > 1.5",
    x        = "Log2 Fold Change",
    y        = "-Log10 P-value",
    colour   = "Expression",
    caption  = paste0("Up in Diabetes: ",
                      sum(diab_deg$Significance == "Up in Diabetes"),
                      " | Down in Diabetes: ",
                      sum(diab_deg$Significance == "Down in Diabetes"))
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(colour = "grey40"),
    legend.position = "bottom"
  )

dev.off()
cat("✓ Volcano plot saved!\n")

# ============================================================
# STEP 8: HEATMAP — TOP 40 DEGS
# ============================================================

cat("\n--- Generating Heatmap ---\n")

# Get top 40 significant genes
top40_diab <- diab_deg %>%
  filter(Significance != "Not Significant") %>%
  arrange(P.Value) %>%
  head(40)

# Extract expression for top genes
heatmap_diab <- diab_subset[rownames(top40_diab), ]

# Create annotation
annotation_diab <- data.frame(
  Condition = c(rep("Diabetes", length(diabetes_idx)),
                rep("Normal",   length(normal_idx))),
  row.names = colnames(heatmap_diab))

# Define colours
ann_colors_diab <- list(
  Condition = c(
    "Diabetes" = "#E74C3C",
    "Normal"   = "#2ECC71"))

png("results/figures/diabetes/02_Diabetes_heatmap.png",
    width = 1400, height = 1600, res = 120)

pheatmap(heatmap_diab,
         annotation_col    = annotation_diab,
         annotation_colors = ann_colors_diab,
         scale             = "row",
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         color             = colorRampPalette(
           rev(brewer.pal(11, "RdBu")))(100),
         main              = "Top 40 DEGs — Diabetes vs Normal\nGSE15932",
         fontsize          = 9,
         fontsize_row      = 7,
         show_colnames     = FALSE,
         border_color      = NA)

dev.off()
cat("✓ Heatmap saved!\n")

# ============================================================
# STEP 9: PCA PLOT — ALL PATIENT GROUPS
# ============================================================

cat("\n--- Generating PCA Plot ---\n")

# PCA on all groups
pca_diab <- prcomp(t(diabetes_expr), scale. = TRUE)

# Extract coordinates
pca_diab_df <- data.frame(
  PC1   = pca_diab$x[, 1],
  PC2   = pca_diab$x[, 2],
  Group = groups_diabetes
)

# Variance explained
var_exp_diab <- round(summary(pca_diab)$importance[2, 1:2] * 100, 1)
cat("PC1 variance:", var_exp_diab[1], "%\n")
cat("PC2 variance:", var_exp_diab[2], "%\n")

png("results/figures/diabetes/03_Diabetes_PCA.png",
    width = 1400, height = 1000, res = 120)

ggplot(pca_diab_df, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_colour_manual(values = c(
    "Cancer+Diabetes" = "#E74C3C",
    "Diabetes Only"   = "#F39C12",
    "Normal"          = "#2ECC71",
    "Other"           = "#9B59B6")) +
  stat_ellipse(aes(group = Group),
               linetype = "dashed", linewidth = 0.7) +
  labs(
    title    = "PCA — Diabetes Patient Groups",
    subtitle = "GSE15932 | Each point represents one patient sample",
    x        = paste0("PC1 (", var_exp_diab[1], "% variance)"),
    y        = paste0("PC2 (", var_exp_diab[2], "% variance)"),
    colour   = "Patient Group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = 14),
    plot.subtitle   = element_text(colour = "grey40"),
    legend.position = "right"
  )

dev.off()
cat("✓ PCA plot saved!\n")

# ============================================================
# STEP 10: SAVE RESULTS
# ============================================================

cat("\n--- Saving Results ---\n")

write.csv(diab_deg,
          "results/data/Diabetes_DEGs_Diabetes_vs_Normal.csv",
          row.names = TRUE)

cat("✓ Results saved!\n")

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("DIABETES ANALYSIS COMPLETE!\n")
cat("============================================================\n")
cat("Dataset:          GSE15932\n")
cat("Comparison:       Diabetes Only vs Normal\n")
cat("Up in Diabetes:  ", sum(diab_deg$Significance == "Up in Diabetes"), "\n")
cat("Down in Diabetes:", sum(diab_deg$Significance == "Down in Diabetes"), "\n")
cat("\nFigures saved in: results/figures/diabetes/\n")
cat("  01_Diabetes_volcano.png\n")
cat("  02_Diabetes_heatmap.png\n")
cat("  03_Diabetes_PCA.png\n")
cat("\nData saved in: results/data/\n")
cat("  Diabetes_DEGs_Diabetes_vs_Normal.csv\n")
cat("============================================================\n")
cat("Relevance to Precision Medicine:\n")
cat("Identifying key genes dysregulated in diabetes mellitus\n")
cat("contributes to biomarker discovery and personalised\n")
cat("treatment strategies in precision medicine.\n")
cat("============================================================\n")
