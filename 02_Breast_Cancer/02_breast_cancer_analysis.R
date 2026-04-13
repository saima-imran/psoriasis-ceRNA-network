# ============================================================
# BREAST CANCER GENE EXPRESSION ANALYSIS
# ============================================================
# Author:   Saima Imran
# GitHub:   https://github.com/saima-imran
# Dataset:  GSE45827 (NCBI GEO)
# Disease:  Breast Cancer Subtypes
# ============================================================
# Description:
# This script analyses gene expression differences across
# breast cancer subtypes using public GEO microarray data.
# Subtypes include: Normal, Basal (TNBC), Luminal A,
# Luminal B, and HER2-enriched.
# Directly relevant to precision medicine and diagnostics.
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

cat("✓ All packages loaded!\n")

# ============================================================
# STEP 2: CREATE FOLDER STRUCTURE
# ============================================================

dir.create("results",                    showWarnings = FALSE)
dir.create("results/figures",            showWarnings = FALSE)
dir.create("results/figures/breast_cancer", showWarnings = FALSE)
dir.create("results/data",               showWarnings = FALSE)

cat("✓ Folders created!\n")

# ============================================================
# STEP 3: DOWNLOAD BREAST CANCER DATA FROM NCBI GEO
# ============================================================

cat("\n--- Downloading GSE45827 Breast Cancer Dataset ---\n")
cat("Dataset: Breast cancer subtypes microarray data\n")
cat("Subtypes: Normal, Basal, Luminal A, Luminal B, HER2\n\n")

gse_bc <- getGEO("GSE45827", GSEMatrix = TRUE, getGPL = FALSE)

cat("✓ Dataset downloaded!\n")
cat("Platforms found:", length(gse_bc), "\n")

# ============================================================
# STEP 4: EXTRACT EXPRESSION DATA
# ============================================================

cat("\n--- Extracting Expression Data ---\n")

# Extract expression matrix
bc_expr <- exprs(gse_bc[[1]])

cat("Expression matrix dimensions:\n")
cat("Genes x Samples:", nrow(bc_expr), "x", ncol(bc_expr), "\n")

# Extract sample information
sample_data <- pData(gse_bc[[1]])

cat("\nSample information columns available:\n")
cat(colnames(sample_data), "\n")

# ============================================================
# STEP 5: EXTRACT CANCER SUBTYPES
# ============================================================

cat("\n--- Extracting Cancer Subtypes ---\n")

# Extract subtype information
# Look for subtype in characteristics columns
subtype_col <- grep("subtype|type|characteristics", 
                    colnames(sample_data), 
                    ignore.case = TRUE, value = TRUE)[1]

cat("Using column:", subtype_col, "\n")

subtypes <- sample_data[[subtype_col]]
cat("Subtypes found:\n")
print(table(subtypes))

# Clean subtype names
subtypes_clean <- gsub("subtype: ", "", subtypes)
subtypes_clean <- gsub("cell type: ", "", subtypes_clean)
subtypes_clean <- trimws(subtypes_clean)

cat("\nCleaned subtypes:\n")
print(table(subtypes_clean))

# ============================================================
# STEP 6: QUALITY CONTROL - BOXPLOT
# ============================================================

cat("\n--- Generating QC Boxplot ---\n")

# Define colours for subtypes
subtype_colours <- c(
  "normal"   = "#2ECC71",
  "basal"    = "#E74C3C", 
  "luminal A" = "#3498DB",
  "luminal B" = "#9B59B6",
  "HER2"     = "#F39C12",
  "Normal"   = "#2ECC71",
  "Basal"    = "#E74C3C",
  "LumA"     = "#3498DB",
  "LumB"     = "#9B59B6",
  "Her2"     = "#F39C12"
)

# Assign colours to samples
sample_cols <- subtype_colours[subtypes_clean]
sample_cols[is.na(sample_cols)] <- "#AAAAAA"

png("results/figures/breast_cancer/01_BC_expression_boxplot.png",
    width = 1600, height = 900, res = 120)

par(mar = c(8, 4, 4, 2))
boxplot(bc_expr,
        main     = "GSE45827 Breast Cancer Gene Expression",
        ylab     = "Log2 Expression",
        col      = sample_cols,
        las      = 2,
        cex.axis = 0.5,
        outline  = FALSE)

# Add legend
unique_subtypes <- unique(subtypes_clean)
legend_cols <- subtype_colours[unique_subtypes]
legend_cols <- legend_cols[!is.na(legend_cols)]

legend("topright",
       legend = names(legend_cols),
       fill   = legend_cols,
       cex    = 0.8,
       title  = "Cancer Subtype")

dev.off()
cat("✓ Boxplot saved!\n")

# ============================================================
# STEP 7: DIFFERENTIAL EXPRESSION — BASAL vs NORMAL
# ============================================================

cat("\n--- Differential Expression: Basal (TNBC) vs Normal ---\n")
cat("Focus: Triple Negative Breast Cancer vs Normal tissue\n\n")

# Select Basal and Normal samples
basal_idx  <- grep("basal|Basal", subtypes_clean, ignore.case = TRUE)
normal_idx <- grep("normal|Normal", subtypes_clean, ignore.case = TRUE)

cat("Basal samples:", length(basal_idx), "\n")
cat("Normal samples:", length(normal_idx), "\n")

# Subset expression data
bc_subset <- bc_expr[, c(basal_idx, normal_idx)]

# Create groups
groups <- factor(c(rep("Basal", length(basal_idx)),
                   rep("Normal", length(normal_idx))))

# Design matrix
design <- model.matrix(~ 0 + groups)
colnames(design) <- c("Basal", "Normal")

# Contrast
contrast <- makeContrasts(Basal - Normal, levels = design)

# Fit model
fit <- eBayes(contrasts.fit(lmFit(bc_subset, design), contrast))

# Extract results
bc_results <- topTable(fit, coef = 1, number = Inf, sort.by = "P")

# Apply thresholds
bc_deg <- bc_results %>%
  mutate(
    Gene = rownames(.),
    Significance = case_when(
      P.Value < 0.05 & logFC >  1.5 ~ "Up in Basal",
      P.Value < 0.05 & logFC < -1.5 ~ "Down in Basal",
      TRUE ~ "Not Significant"
    )
  )

cat("\nDifferential Expression Results:\n")
cat("Total genes tested:", nrow(bc_results), "\n")
cat("Up in Basal (TNBC):", sum(bc_deg$Significance == "Up in Basal"), "\n")
cat("Down in Basal:", sum(bc_deg$Significance == "Down in Basal"), "\n")

# ============================================================
# STEP 8: VOLCANO PLOT — BASAL vs NORMAL
# ============================================================

cat("\n--- Generating Volcano Plot ---\n")

# Get top genes for labelling
top_label <- bc_deg %>%
  filter(Significance != "Not Significant") %>%
  arrange(P.Value) %>%
  head(15)

png("results/figures/breast_cancer/02_BC_volcano_Basal_vs_Normal.png",
    width = 1400, height = 1000, res = 120)

ggplot(bc_deg, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(colour = Significance), alpha = 0.5, size = 1.5) +
  scale_colour_manual(values = c(
    "Up in Basal"     = "#E74C3C",
    "Down in Basal"   = "#3498DB",
    "Not Significant" = "#AAAAAA")) +
  geom_vline(xintercept = c(-1.5, 1.5), 
             linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed", linewidth = 0.5) +
  geom_text_repel(
    data = top_label,
    aes(label = Gene),
    size = 3,
    max.overlaps = 15,
    colour = "black") +
  labs(
    title    = "Breast Cancer — Basal (TNBC) vs Normal",
    subtitle = "GSE45827 | Threshold: p < 0.05 | |logFC| > 1.5",
    x        = "Log2 Fold Change",
    y        = "-Log10 P-value",
    colour   = "Expression",
    caption  = paste0("Up in Basal: ", 
                      sum(bc_deg$Significance == "Up in Basal"),
                      " | Down in Basal: ",
                      sum(bc_deg$Significance == "Down in Basal"))
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
# STEP 9: HEATMAP — TOP GENES ACROSS ALL SUBTYPES
# ============================================================

cat("\n--- Generating Heatmap across all subtypes ---\n")

# Get top 40 most significant genes
top40 <- bc_deg %>%
  filter(Significance != "Not Significant") %>%
  arrange(P.Value) %>%
  head(40)

# Extract expression for top genes across ALL samples
heatmap_data <- bc_expr[rownames(top40), ]

# Create annotation with subtype colours
annotation_col <- data.frame(
  Subtype = subtypes_clean,
  row.names = colnames(heatmap_data)
)

# Define annotation colours
unique_sub <- unique(subtypes_clean)
ann_colors <- list(
  Subtype = setNames(
    colorRampPalette(brewer.pal(8, "Set1"))(length(unique_sub)),
    unique_sub
  )
)

# Generate heatmap
png("results/figures/breast_cancer/03_BC_heatmap_top40_genes.png",
    width = 1600, height = 1600, res = 120)

pheatmap(heatmap_data,
         annotation_col    = annotation_col,
         annotation_colors = ann_colors,
         scale             = "row",
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         color             = colorRampPalette(
           rev(brewer.pal(11, "RdBu")))(100),
         main              = "Top 40 DEGs — Breast Cancer Subtypes\nGSE45827",
         fontsize          = 9,
         fontsize_row      = 7,
         show_colnames     = FALSE,
         border_color      = NA,
         treeheight_row    = 30,
         treeheight_col    = 30)

dev.off()
cat("✓ Heatmap saved!\n")

# ============================================================
# STEP 10: SAVE RESULTS
# ============================================================

cat("\n--- Saving Results ---\n")

write.csv(bc_deg, 
          "results/data/BreastCancer_DEGs_Basal_vs_Normal.csv",
          row.names = TRUE)

cat("✓ Results saved!\n")

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("BREAST CANCER ANALYSIS COMPLETE!\n")
cat("============================================================\n")
cat("Dataset:        GSE45827\n")
cat("Comparison:     Basal (TNBC) vs Normal\n")
cat("Up in Basal:   ", sum(bc_deg$Significance == "Up in Basal"), "\n")
cat("Down in Basal: ", sum(bc_deg$Significance == "Down in Basal"), "\n")
cat("\nFigures saved in: results/figures/breast_cancer/\n")
cat("  01_BC_expression_boxplot.png\n")
cat("  02_BC_volcano_Basal_vs_Normal.png\n")
cat("  03_BC_heatmap_top40_genes.png\n")
cat("\nData saved in: results/data/\n")
cat("  BreastCancer_DEGs_Basal_vs_Normal.csv\n")
cat("============================================================\n")
cat("This analysis is directly relevant to:\n")
cat("KI DDLS PhD — Bioinformatics: Clinical Breast Cancer Proteomics\n")
cat("Reference: AID 2-1349/2026\n")
cat("============================================================\n")
