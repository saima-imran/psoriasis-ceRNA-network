# ============================================================
# M.Sc. Bioinformatics Assignment 1
# Bioinformatic Analysis with R — BI731A
# ============================================================
# Author:   Saima Imran
# Course:   Bioinformatics Analysis with R
# University: University of Skövde, Sweden
# Dataset:  E-MTAB-7847 (ArrayExpress)
# ============================================================

# ============================================================
# STEP 1: INSTALL AND LOAD PACKAGES
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ArrayExpress", update = FALSE, ask = FALSE)

library(ArrayExpress)
library(Biobase)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

cat("✓ All packages loaded!\n")

# ============================================================
# STEP 2: SET WORKING DIRECTORY
# ============================================================

setwd("C:/Users/saima/OneDrive/Documents/Bioinformatics_portfolio")
cat("✓ Working directory set!\n")

# ============================================================
# STEP 3: DOWNLOAD E-MTAB-7847 FROM ARRAYEXPRESS
# ============================================================

cat("\n--- Downloading E-MTAB-7847 from ArrayExpress ---\n")
cat("Dataset: RNA-seq data — liver cancer vs normal cell lines\n")
cat("Cell lines: HUH7 (cancer) and THLE5B (normal)\n")
cat("Compounds: AD80, DMSO (control), Sorafenib\n\n")

dataset <- getAE("E-MTAB-7847", type = "processed")
cat("✓ Dataset downloaded!\n")

# ============================================================
# STEP 4: READ EXPRESSION FILES
# ============================================================

cat("\n--- Reading Expression Files ---\n")

thle5b <- read.table("THLE5B_DESeq2_normalized_count_matrix.txt",
                      header = TRUE, row.names = 1, sep = "\t")
huh7   <- read.table("HUH7_DESeq2_normalized_count_matrix.txt",
                      header = TRUE, row.names = 1, sep = "\t")

cat("THLE5B dimensions:", nrow(thle5b), "x", ncol(thle5b), "\n")
cat("HUH7 dimensions:",   nrow(huh7),   "x", ncol(huh7),   "\n")

# ============================================================
# STEP 5: DEFINE EXPERIMENTAL DESIGN
# ============================================================

cat("\n--- Defining Experimental Design ---\n")

# Combine both datasets
all_samples <- cbind(thle5b, huh7)

# Define factor vectors
cell_line <- factor(c(rep("THLE5B", 9), rep("HUH7", 9)))
compound  <- factor(c(rep("AD80", 3), rep("DMSO", 3), rep("Sorafenib", 3),
                      rep("AD80", 3), rep("DMSO", 3), rep("Sorafenib", 3)))

# Create design table
design_table <- data.frame(
  CellLine = cell_line,
  Compound = compound,
  row.names = colnames(all_samples)
)

cat("Experimental design table:\n")
print(design_table)

# ============================================================
# STEP 6: LOG2 TRANSFORMATION
# ============================================================

cat("\n--- Log2 Transforming Data ---\n")

all_log <- log2(all_samples + 1)
cat("✓ Log2 transformation complete!\n")

# ============================================================
# STEP 7: QUALITY CONTROL — BOXPLOT
# ============================================================

cat("\n--- Generating QC Boxplot ---\n")

sample_cols <- c(rep("#3498DB", 9), rep("#E74C3C", 9))

png("Assignment1_boxplot.png", width = 1400, height = 800, res = 120)
par(mar = c(8, 4, 4, 2))
boxplot(all_log,
        main     = "E-MTAB-7847 — RNA-seq Expression Quality Control",
        ylab     = "Log2 Normalised Counts",
        col      = sample_cols,
        las      = 2,
        cex.axis = 0.7,
        outline  = FALSE)
legend("topright",
       legend = c("THLE5B (Normal liver)", "HUH7 (Liver cancer)"),
       fill   = c("#3498DB", "#E74C3C"),
       cex    = 0.8)
dev.off()
cat("✓ Boxplot saved!\n")

# ============================================================
# STEP 8: QUALITY CONTROL — HEATMAP
# ============================================================

cat("\n--- Generating Heatmap ---\n")

# Select top 50 most variable genes
gene_var  <- apply(all_log, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE)[1:50])
heatmap_data <- all_log[top_genes, ]

# Sample annotation
annotation_col <- data.frame(
  CellLine = cell_line,
  Compound = compound,
  row.names = colnames(all_log))

# Colours
ann_colors <- list(
  CellLine = c("THLE5B" = "#3498DB", "HUH7" = "#E74C3C"),
  Compound = c("AD80"      = "#F39C12",
               "DMSO"      = "#9B59B6",
               "Sorafenib" = "#2ECC71"))

png("Assignment1_heatmap.png", width = 1400, height = 1600, res = 120)

pheatmap(heatmap_data,
         annotation_col    = annotation_col,
         annotation_colors = ann_colors,
         scale             = "row",
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         color             = colorRampPalette(
           rev(brewer.pal(11, "RdBu")))(100),
         main              = "Top 50 Variable Genes — E-MTAB-7847\nRNA-seq Quality Control",
         fontsize          = 9,
         fontsize_row      = 6,
         show_colnames     = TRUE,
         border_color      = NA)

dev.off()
cat("✓ Heatmap saved!\n")

# ============================================================
# STEP 9: QUALITY CONTROL — PCA PLOT
# ============================================================

cat("\n--- Generating PCA Plot ---\n")

# Remove zero variance genes
all_log_clean <- all_log[apply(all_log, 1, var) > 0, ]
cat("Genes after cleaning:", nrow(all_log_clean), "\n")

# PCA
pca_assign <- prcomp(t(all_log_clean), scale. = TRUE)

# Extract coordinates
pca_df_assign <- data.frame(
  PC1      = pca_assign$x[, 1],
  PC2      = pca_assign$x[, 2],
  CellLine = cell_line,
  Compound = compound
)

# Variance explained
var_exp_assign <- round(summary(pca_assign)$importance[2, 1:2] * 100, 1)
cat("PC1 variance:", var_exp_assign[1], "%\n")
cat("PC2 variance:", var_exp_assign[2], "%\n")

png("Assignment1_PCA.png", width = 1400, height = 1000, res = 120)

ggplot(pca_df_assign, aes(x = PC1, y = PC2,
                           colour = CellLine,
                           shape  = Compound)) +
  geom_point(size = 5, alpha = 0.9) +
  scale_colour_manual(values = c(
    "THLE5B" = "#3498DB",
    "HUH7"   = "#E74C3C")) +
  scale_shape_manual(values = c(
    "AD80"      = 16,
    "DMSO"      = 17,
    "Sorafenib" = 15)) +
  labs(
    title    = "PCA — E-MTAB-7847 RNA-seq Quality Control",
    subtitle = "Colour = Cell Line | Shape = Compound Treatment",
    x        = paste0("PC1 (", var_exp_assign[1], "% variance)"),
    y        = paste0("PC2 (", var_exp_assign[2], "% variance)"),
    colour   = "Cell Line",
    shape    = "Compound"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(colour = "grey40"),
    legend.position = "right"
  )

dev.off()
cat("✓ PCA plot saved!\n")

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("ASSIGNMENT 1 ANALYSIS COMPLETE!\n")
cat("============================================================\n")
cat("Dataset:      E-MTAB-7847 (ArrayExpress)\n")
cat("Cell lines:   THLE5B (Normal) + HUH7 (Cancer)\n")
cat("Compounds:    AD80, DMSO, Sorafenib\n")
cat("Total samples:", ncol(all_samples), "\n")
cat("Total genes:  ", nrow(all_log), "\n")
cat("\nFigures saved:\n")
cat("  Assignment1_boxplot.png\n")
cat("  Assignment1_heatmap.png\n")
cat("  Assignment1_PCA.png\n")
cat("============================================================\n")
