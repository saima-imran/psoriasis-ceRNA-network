# ============================================================
# Script 01: GEO Data Download and Preprocessing
# ============================================================
# Project: Comprehensive Analysis of lncRNA and circRNA 
#          Mediated ceRNA Network in Psoriasis
# Author:  Saima Imran
# GitHub:  https://github.com/saima-imran
# Date:    2022
# Dataset: GSE145305 (NCBI GEO)
# ============================================================
# Description:
# This script downloads the GSE145305 microarray dataset from 
# NCBI GEO, reads the raw CEL files using the oligo package,
# and performs RMA normalisation for both mRNA and miRNA data.
# Output: normalised expression matrices ready for downstream
# differential expression analysis.
# ============================================================

# ------------------------------------------------------------
# STEP 1: Install and load required packages
# ------------------------------------------------------------

# Install Bioconductor manager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required Bioconductor packages
BiocManager::install(c("GEOquery", "oligo", "Biobase"), 
                     update = FALSE, ask = FALSE)

# Install CRAN packages
install.packages(c("tidyverse", "ggplot2"))

# Load libraries
library(GEOquery)      # Download GEO datasets
library(oligo)         # Read and preprocess CEL files
library(Biobase)       # ExpressionSet handling
library(tidyverse)     # Data manipulation
library(ggplot2)       # Visualisation

cat("✓ All packages loaded successfully\n")

# ------------------------------------------------------------
# STEP 2: Download GSE145305 from NCBI GEO
# ------------------------------------------------------------

cat("\n--- Downloading GSE145305 from NCBI GEO ---\n")

# Download the dataset (this may take a few minutes)
gse <- getGEO("GSE145305", GSEMatrix = TRUE, getGPL = FALSE)

cat("✓ Dataset downloaded successfully\n")
cat("Number of platforms found:", length(gse), "\n")

# The dataset contains two platforms:
# GPL17586 - mRNA (HTA-2_0 Affymetrix Human Transcriptome Array)
# GPL19117 - miRNA (Affymetrix Multispecies miRNA-4 Array)

# Extract platform information
for(i in seq_along(gse)){
  cat("\nPlatform", i, ":", annotation(gse[[i]]), "\n")
  cat("Samples:", ncol(gse[[i]]), "\n")
}

# ------------------------------------------------------------
# STEP 3: Download raw CEL files
# ------------------------------------------------------------

cat("\n--- Downloading raw CEL files ---\n")

# Create directory for raw data
dir.create("raw_data", showWarnings = FALSE)
dir.create("raw_data/mRNA", showWarnings = FALSE)
dir.create("raw_data/miRNA", showWarnings = FALSE)

# Download supplementary files (CEL files)
getGEOSuppFiles("GSE145305", baseDir = "raw_data/")

cat("✓ Raw CEL files downloaded\n")

# ------------------------------------------------------------
# STEP 4: Read mRNA CEL files with oligo
# ------------------------------------------------------------

cat("\n--- Reading mRNA CEL files ---\n")

# List mRNA CEL files (GPL17586 platform - HTA array)
mrna_cel_files <- list.celfiles("raw_data/GSE145305/", 
                                 full.names = TRUE,
                                 listGzipped = TRUE)

# Filter for mRNA files (HTA platform)
# Sample IDs from thesis: GSM4314529, GSM4314530, GSM4314531 (Psoriasis)
#                         GSM4314532, GSM4314533, GSM4314534 (Healthy)
mrna_files <- mrna_cel_files[grep("GSM431", mrna_cel_files)]

cat("mRNA CEL files found:", length(mrna_files), "\n")
cat("Files:", basename(mrna_files), "\n")

# Read CEL files into oligo object
mrna_raw <- read.celfiles(mrna_files)

cat("✓ mRNA CEL files read successfully\n")
cat("Dimensions (probes x samples):", dim(mrna_raw), "\n")

# ------------------------------------------------------------
# STEP 5: Read miRNA CEL files with oligo
# ------------------------------------------------------------

cat("\n--- Reading miRNA CEL files ---\n")

# Filter for miRNA files (miRNA-4 platform)
# Sample IDs: GSM4305469-GSM4305474 (4 psoriasis + 4 healthy)
mirna_files <- mrna_cel_files[grep("GSM430", mrna_cel_files)]

cat("miRNA CEL files found:", length(mirna_files), "\n")

# Read miRNA CEL files
mirna_raw <- read.celfiles(mirna_files)

cat("✓ miRNA CEL files read successfully\n")
cat("Dimensions (probes x samples):", dim(mirna_raw), "\n")

# ------------------------------------------------------------
# STEP 6: Quality Control - Boxplots before normalisation
# ------------------------------------------------------------

cat("\n--- Generating QC boxplots BEFORE normalisation ---\n")

dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)

# mRNA boxplot before normalisation
png("results/figures/01_mRNA_boxplot_before_normalisation.png",
    width = 1200, height = 800, res = 120)
boxplot(mrna_raw, 
        main = "GSE145305 mRNA - Before Normalisation",
        xlab = "Samples",
        ylab = "Log2 Expression",
        col = c(rep("red", 3), rep("blue", 3)),
        las = 2,
        cex.axis = 0.7)
legend("topright", 
       legend = c("Psoriasis", "Healthy Control"),
       fill = c("red", "blue"))
dev.off()

# miRNA boxplot before normalisation
png("results/figures/01_miRNA_boxplot_before_normalisation.png",
    width = 1200, height = 800, res = 120)
boxplot(mirna_raw,
        main = "GSE145305 miRNA - Before Normalisation",
        xlab = "Samples", 
        ylab = "Log2 Expression",
        col = c(rep("red", 4), rep("blue", 4)),
        las = 2,
        cex.axis = 0.7)
legend("topright",
       legend = c("Psoriasis", "Healthy Control"),
       fill = c("red", "blue"))
dev.off()

cat("✓ QC boxplots saved to results/figures/\n")

# ------------------------------------------------------------
# STEP 7: RMA Normalisation
# ------------------------------------------------------------

cat("\n--- Performing RMA Normalisation ---\n")
cat("RMA performs three steps:\n")
cat("  1. Background correction\n")
cat("  2. Quantile normalisation\n")
cat("  3. Summarisation\n\n")

# Normalise mRNA data using RMA
mrna_normalised <- rma(mrna_raw)

cat("✓ mRNA RMA normalisation complete\n")
cat("Dimensions after normalisation:", dim(mrna_normalised), "\n")

# Normalise miRNA data using RMA
mirna_normalised <- rma(mirna_raw)

cat("✓ miRNA RMA normalisation complete\n")
cat("Dimensions after normalisation:", dim(mirna_normalised), "\n")

# ------------------------------------------------------------
# STEP 8: Quality Control - Boxplots after normalisation
# ------------------------------------------------------------

cat("\n--- Generating QC boxplots AFTER normalisation ---\n")

# Extract expression matrices
mrna_expr <- exprs(mrna_normalised)
mirna_expr <- exprs(mirna_normalised)

# mRNA boxplot AFTER normalisation
png("results/figures/01_mRNA_boxplot_after_normalisation.png",
    width = 1200, height = 800, res = 120)
boxplot(mrna_expr,
        main = "GSE145305 mRNA - After RMA Normalisation",
        xlab = "Samples",
        ylab = "Log2 Expression",
        col = c(rep("red", 3), rep("blue", 3)),
        las = 2,
        cex.axis = 0.7)
legend("topright",
       legend = c("Psoriasis", "Healthy Control"),
       fill = c("red", "blue"))
dev.off()

# miRNA boxplot AFTER normalisation
png("results/figures/01_miRNA_boxplot_after_normalisation.png",
    width = 1200, height = 800, res = 120)
boxplot(mirna_expr,
        main = "GSE145305 miRNA - After RMA Normalisation",
        xlab = "Samples",
        ylab = "Log2 Expression",
        col = c(rep("red", 4), rep("blue", 4)),
        las = 2,
        cex.axis = 0.7)
legend("topright",
       legend = c("Psoriasis", "Healthy Control"),
       fill = c("red", "blue"))
dev.off()

cat("✓ Post-normalisation QC boxplots saved\n")

# ------------------------------------------------------------
# STEP 9: Save normalised data for downstream analysis
# ------------------------------------------------------------

cat("\n--- Saving normalised expression data ---\n")

dir.create("results/data", showWarnings = FALSE)

# Save as RDS for use in Script 02
saveRDS(mrna_normalised, "results/data/mrna_normalised.rds")
saveRDS(mirna_normalised, "results/data/mirna_normalised.rds")

# Save as CSV for inspection
write.csv(mrna_expr, "results/data/mrna_expression_matrix.csv")
write.csv(mirna_expr, "results/data/mirna_expression_matrix.csv")

cat("✓ Normalised data saved to results/data/\n")

# ------------------------------------------------------------
# STEP 10: Summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("PREPROCESSING COMPLETE - SUMMARY\n")
cat("============================================================\n")
cat("Dataset:          GSE145305\n")
cat("mRNA samples:     6 (3 Psoriasis + 3 Healthy Control)\n")
cat("miRNA samples:    8 (4 Psoriasis + 4 Healthy Control)\n")
cat("mRNA probes:     ", nrow(mrna_expr), "\n")
cat("miRNA probes:    ", nrow(mirna_expr), "\n")
cat("Normalisation:    RMA (Background correction + Quantile normalisation)\n")
cat("Output files:     results/data/ and results/figures/\n")
cat("Next step:        Run 02_differential_expression.R\n")
cat("============================================================\n")
