# ============================================================
# Script 02: Differential Expression Analysis
# ============================================================
# Project: Comprehensive Analysis of lncRNA and circRNA 
#          Mediated ceRNA Network in Psoriasis
# Author:  Saima Imran
# GitHub:  https://github.com/saima-imran
# Date:    2022
# ============================================================
# Description:
# This script performs differential expression analysis using
# the limma package on the RMA-normalised mRNA and miRNA data
# from GSE145305. Volcano plots are generated to visualise
# significantly up- and down-regulated genes.
#
# Thresholds used (from thesis):
#   mRNA:  p-value < 0.05, |logFC| > 2
#   miRNA: p-value < 0.05, |logFC| > 1.5
# ============================================================

# ------------------------------------------------------------
# STEP 1: Load required packages
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("limma", "EnhancedVolcano"), 
                     update = FALSE, ask = FALSE)

install.packages(c("tidyverse", "ggplot2", "ggrepel"))

library(limma)            # Differential expression analysis
library(EnhancedVolcano)  # Professional volcano plots
library(tidyverse)        # Data manipulation
library(ggplot2)          # Visualisation
library(ggrepel)          # Non-overlapping labels

cat("✓ All packages loaded successfully\n")

# ------------------------------------------------------------
# STEP 2: Load normalised data from Script 01
# ------------------------------------------------------------

cat("\n--- Loading normalised data ---\n")

mrna_normalised <- readRDS("results/data/mrna_normalised.rds")
mirna_normalised <- readRDS("results/data/mirna_normalised.rds")

mrna_expr <- exprs(mrna_normalised)
mirna_expr <- exprs(mirna_normalised)

cat("✓ mRNA data loaded:", nrow(mrna_expr), "probes x", 
    ncol(mrna_expr), "samples\n")
cat("✓ miRNA data loaded:", nrow(mirna_expr), "probes x", 
    ncol(mirna_expr), "samples\n")

# ------------------------------------------------------------
# STEP 3: Define experimental groups
# ------------------------------------------------------------

cat("\n--- Defining experimental groups ---\n")

# mRNA: 3 Psoriasis (PS) + 3 Healthy Skin (HS)
# Sample order: PS, PS, PS, HS, HS, HS
mrna_groups <- factor(c("Psoriasis", "Psoriasis", "Psoriasis",
                         "Healthy",   "Healthy",   "Healthy"))

# miRNA: 4 Psoriasis (PS) + 4 Healthy Skin (HS)
mirna_groups <- factor(c("Psoriasis", "Psoriasis", "Psoriasis", "Psoriasis",
                          "Healthy",   "Healthy",   "Healthy",   "Healthy"))

cat("mRNA groups:", as.character(mrna_groups), "\n")
cat("miRNA groups:", as.character(mirna_groups), "\n")

# ------------------------------------------------------------
# STEP 4: Create design matrices
# ------------------------------------------------------------

cat("\n--- Creating design matrices ---\n")

# mRNA design matrix
mrna_design <- model.matrix(~ 0 + mrna_groups)
colnames(mrna_design) <- c("Healthy", "Psoriasis")

# miRNA design matrix  
mirna_design <- model.matrix(~ 0 + mirna_groups)
colnames(mirna_design) <- c("Healthy", "Psoriasis")

cat("✓ Design matrices created\n")

# Define contrasts: Psoriasis vs Healthy
mrna_contrast  <- makeContrasts(Psoriasis - Healthy, levels = mrna_design)
mirna_contrast <- makeContrasts(Psoriasis - Healthy, levels = mirna_design)

# ------------------------------------------------------------
# STEP 5: Fit linear models with limma
# ------------------------------------------------------------

cat("\n--- Fitting linear models with limma ---\n")
cat("limma uses empirical Bayes to improve variance estimation\n")
cat("especially important for small sample sizes\n\n")

# mRNA: fit linear model
mrna_fit  <- lmFit(mrna_expr, mrna_design)
mrna_fit2 <- contrasts.fit(mrna_fit, mrna_contrast)
mrna_fit2 <- eBayes(mrna_fit2)   # Empirical Bayes moderation

cat("✓ mRNA linear model fitted\n")

# miRNA: fit linear model
mirna_fit  <- lmFit(mirna_expr, mirna_design)
mirna_fit2 <- contrasts.fit(mirna_fit, mirna_contrast)
mirna_fit2 <- eBayes(mirna_fit2)

cat("✓ miRNA linear model fitted\n")

# ------------------------------------------------------------
# STEP 6: Extract differentially expressed genes
# ------------------------------------------------------------

cat("\n--- Extracting differentially expressed genes ---\n")

# mRNA: threshold p < 0.05, |logFC| > 2
mrna_results <- topTable(mrna_fit2, 
                          coef = 1, 
                          number = Inf,
                          sort.by = "P")

# Apply thresholds
mrna_deg <- mrna_results %>%
  filter(P.Value < 0.05 & abs(logFC) > 2) %>%
  mutate(Regulation = ifelse(logFC > 2, "Up-regulated", "Down-regulated"))

cat("\nmRNA Differential Expression Results:\n")
cat("Total probes tested:", nrow(mrna_results), "\n")
cat("Significant DEGs (p<0.05, |logFC|>2):", nrow(mrna_deg), "\n")
cat("Up-regulated:", sum(mrna_deg$Regulation == "Up-regulated"), "\n")
cat("Down-regulated:", sum(mrna_deg$Regulation == "Down-regulated"), "\n")

# miRNA: threshold p < 0.05, |logFC| > 1.5
mirna_results <- topTable(mirna_fit2,
                           coef = 1,
                           number = Inf,
                           sort.by = "P")

mirna_deg <- mirna_results %>%
  filter(P.Value < 0.05 & abs(logFC) > 1.5) %>%
  mutate(Regulation = ifelse(logFC > 1.5, "Up-regulated", "Down-regulated"))

cat("\nmiRNA Differential Expression Results:\n")
cat("Total probes tested:", nrow(mirna_results), "\n")
cat("Significant DEMs (p<0.05, |logFC|>1.5):", nrow(mirna_deg), "\n")
cat("Up-regulated:", sum(mirna_deg$Regulation == "Up-regulated"), "\n")
cat("Down-regulated:", sum(mirna_deg$Regulation == "Down-regulated"), "\n")

# ------------------------------------------------------------
# STEP 7: Generate Volcano Plots
# ------------------------------------------------------------

cat("\n--- Generating Volcano Plots ---\n")

# Add significance column for colouring
mrna_results <- mrna_results %>%
  mutate(Significance = case_when(
    P.Value < 0.05 & logFC >  2 ~ "Up-regulated",
    P.Value < 0.05 & logFC < -2 ~ "Down-regulated",
    TRUE ~ "Not Significant"
  ))

mirna_results <- mirna_results %>%
  mutate(Significance = case_when(
    P.Value < 0.05 & logFC >  1.5 ~ "Up-regulated",
    P.Value < 0.05 & logFC < -1.5 ~ "Down-regulated",
    TRUE ~ "Not Significant"
  ))

# --- mRNA Volcano Plot ---
png("results/figures/02_mRNA_volcano_plot.png",
    width = 1400, height = 1000, res = 120)

ggplot(mrna_results, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(colour = Significance), alpha = 0.6, size = 1.5) +
  scale_colour_manual(values = c(
    "Up-regulated"   = "#E41A1C",
    "Down-regulated" = "#377EB8", 
    "Not Significant" = "#AAAAAA"
  )) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", 
             colour = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "black", linewidth = 0.5) +
  geom_text_repel(
    data = subset(mrna_results, 
                  Significance != "Not Significant" & -log10(P.Value) > 5),
    aes(label = rownames(subset(mrna_results, 
                                Significance != "Not Significant" & 
                                  -log10(P.Value) > 5))),
    size = 3, max.overlaps = 20
  ) +
  labs(
    title = "GSE145305 mRNA — Differential Expression",
    subtitle = "Psoriasis vs Healthy Skin | Threshold: p<0.05, |logFC|>2",
    x = "Log2 Fold Change",
    y = "-Log10 P-value",
    colour = "Expression",
    caption = paste0("Total = ", nrow(mrna_results), " variables | ",
                     "Up: ", sum(mrna_results$Significance == "Up-regulated"),
                     " | Down: ", sum(mrna_results$Significance == "Down-regulated"))
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

dev.off()
cat("✓ mRNA volcano plot saved\n")

# --- miRNA Volcano Plot ---
png("results/figures/02_miRNA_volcano_plot.png",
    width = 1400, height = 1000, res = 120)

ggplot(mirna_results, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(colour = Significance), alpha = 0.6, size = 1.5) +
  scale_colour_manual(values = c(
    "Up-regulated"    = "#E41A1C",
    "Down-regulated"  = "#377EB8",
    "Not Significant" = "#AAAAAA"
  )) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed",
             colour = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "black", linewidth = 0.5) +
  geom_text_repel(
    data = subset(mirna_results,
                  Significance != "Not Significant" & -log10(P.Value) > 4),
    aes(label = rownames(subset(mirna_results,
                                Significance != "Not Significant" &
                                  -log10(P.Value) > 4))),
    size = 3, max.overlaps = 20
  ) +
  labs(
    title = "GSE145305 miRNA — Differential Expression",
    subtitle = "Psoriasis vs Healthy Skin | Threshold: p<0.05, |logFC|>1.5",
    x = "Log2 Fold Change",
    y = "-Log10 P-value",
    colour = "Expression",
    caption = paste0("Total = ", nrow(mirna_results), " variables | ",
                     "Up: ", sum(mirna_results$Significance == "Up-regulated"),
                     " | Down: ", sum(mirna_results$Significance == "Down-regulated"))
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

dev.off()
cat("✓ miRNA volcano plot saved\n")

# ------------------------------------------------------------
# STEP 8: Save results
# ------------------------------------------------------------

cat("\n--- Saving DEG results ---\n")

# Save full results
write.csv(mrna_results, "results/data/mrna_full_results.csv")
write.csv(mirna_results, "results/data/mirna_full_results.csv")

# Save significant DEGs only
write.csv(mrna_deg, "results/data/mrna_significant_DEGs.csv",
          row.names = TRUE)
write.csv(mirna_deg, "results/data/mirna_significant_DEMs.csv",
          row.names = TRUE)

# Save RDS for Script 03
saveRDS(mrna_deg, "results/data/mrna_deg.rds")
saveRDS(mirna_deg, "results/data/mirna_deg.rds")

cat("✓ All results saved to results/data/\n")

# ------------------------------------------------------------
# STEP 9: Summary table of top genes
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TOP 10 UP-REGULATED mRNA GENES\n")
cat("============================================================\n")
print(head(mrna_deg %>% filter(Regulation == "Up-regulated") %>% 
             arrange(desc(logFC)) %>% 
             select(logFC, P.Value, Regulation), 10))

cat("\n============================================================\n")
cat("TOP 10 DOWN-REGULATED mRNA GENES\n")
cat("============================================================\n")
print(head(mrna_deg %>% filter(Regulation == "Down-regulated") %>%
             arrange(logFC) %>%
             select(logFC, P.Value, Regulation), 10))

cat("\n============================================================\n")
cat("DIFFERENTIAL EXPRESSION ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("Next step: Run 03_mirna_target_prediction.R\n")
cat("============================================================\n")
