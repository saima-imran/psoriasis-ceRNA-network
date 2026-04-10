# ============================================================
# Script 03: miRNA-mRNA Target Prediction
# ============================================================
# Project: Comprehensive Analysis of lncRNA and circRNA 
#          Mediated ceRNA Network in Psoriasis
# Author:  Saima Imran
# GitHub:  https://github.com/saima-imran
# Date:    2022
# ============================================================
# Description:
# This script performs miRNA-mRNA target prediction using the
# miRNAtap R package. Only oppositely expressed mRNA-miRNA 
# pairs are considered:
#   - Up-regulated mRNA paired with Down-regulated miRNA
#   - Down-regulated mRNA paired with Up-regulated miRNA
# Results are saved for ceRNA network construction in Cytoscape
# ============================================================

# ------------------------------------------------------------
# STEP 1: Load required packages
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("miRNAtap", "miRNAtap.db", "AnnotationDbi"),
                     update = FALSE, ask = FALSE)

install.packages(c("tidyverse", "ggplot2"))

library(miRNAtap)      # miRNA target prediction
library(miRNAtap.db)   # miRNA target database
library(tidyverse)     # Data manipulation
library(ggplot2)       # Visualisation

cat("✓ All packages loaded successfully\n")

# ------------------------------------------------------------
# STEP 2: Load DEG results from Script 02
# ------------------------------------------------------------

cat("\n--- Loading differential expression results ---\n")

mrna_deg  <- readRDS("results/data/mrna_deg.rds")
mirna_deg <- readRDS("results/data/mirna_deg.rds")

cat("✓ mRNA DEGs loaded:", nrow(mrna_deg), "genes\n")
cat("✓ miRNA DEMs loaded:", nrow(mirna_deg), "miRNAs\n")

# Separate up and down regulated
mrna_up   <- mrna_deg %>% filter(Regulation == "Up-regulated")
mrna_down <- mrna_deg %>% filter(Regulation == "Down-regulated")
mirna_up  <- mirna_deg %>% filter(Regulation == "Up-regulated")
mirna_down <- mirna_deg %>% filter(Regulation == "Down-regulated")

cat("\nmRNA up-regulated:", nrow(mrna_up), "\n")
cat("mRNA down-regulated:", nrow(mrna_down), "\n")
cat("miRNA up-regulated:", nrow(mirna_up), "\n")
cat("miRNA down-regulated:", nrow(mirna_down), "\n")

# ------------------------------------------------------------
# STEP 3: Create opposite expression pairs
# ------------------------------------------------------------

cat("\n--- Creating opposite expression pairs ---\n")
cat("Biological rationale:\n")
cat("  Up-regulated miRNA → suppresses → Down-regulated mRNA\n")
cat("  Down-regulated miRNA → releases → Up-regulated mRNA\n\n")

# Pair 1: Up-regulated mRNA + Down-regulated miRNA
up_mrna_names   <- rownames(mrna_up)
down_mirna_names <- rownames(mirna_down)

# Pair 2: Down-regulated mRNA + Up-regulated miRNA
down_mrna_names <- rownames(mrna_down)
up_mirna_names  <- rownames(mirna_up)

cat("Pair 1 (Up mRNA + Down miRNA):", 
    length(up_mrna_names), "mRNAs x", length(down_mirna_names), "miRNAs\n")
cat("Pair 2 (Down mRNA + Up miRNA):", 
    length(down_mrna_names), "mRNAs x", length(up_mirna_names), "miRNAs\n")

# ------------------------------------------------------------
# STEP 4: miRNA Target Prediction with miRNAtap
# ------------------------------------------------------------

cat("\n--- Running miRNA target prediction ---\n")
cat("Using miRNAtap: aggregates predictions from multiple databases\n")
cat("Sources: DIANA, Miranda, miRDB, miRWalk, Targetscan\n\n")

dir.create("results/target_prediction", showWarnings = FALSE)

# Function to get miRNA targets
get_mirna_targets <- function(mirna_list, min_sources = 2) {
  
  all_targets <- list()
  
  for(mirna in mirna_list) {
    tryCatch({
      # Get predictions from multiple sources
      targets <- getPredictedTargets(
        mirna,
        species = "hsa",        # Homo sapiens
        method = "geom",        # Geometric mean aggregation
        min_src = min_sources   # Minimum sources agreeing
      )
      
      if(!is.null(targets) && nrow(targets) > 0) {
        targets$miRNA <- mirna
        all_targets[[mirna]] <- targets
        cat("  ✓", mirna, ":", nrow(targets), "targets found\n")
      }
      
    }, error = function(e) {
      cat("  ✗", mirna, ": No predictions available\n")
    })
  }
  
  if(length(all_targets) > 0) {
    return(bind_rows(all_targets))
  } else {
    return(NULL)
  }
}

# Run prediction for down-regulated miRNAs
cat("\nPredicting targets for DOWN-regulated miRNAs:\n")
down_mirna_targets <- get_mirna_targets(down_mirna_names, min_sources = 2)

# Run prediction for up-regulated miRNAs  
cat("\nPredicting targets for UP-regulated miRNAs:\n")
up_mirna_targets <- get_mirna_targets(up_mirna_names, min_sources = 2)

cat("\n✓ Target prediction complete\n")

# ------------------------------------------------------------
# STEP 5: Filter for opposite expression pairs
# ------------------------------------------------------------

cat("\n--- Filtering for opposite expression pairs ---\n")

# Pair 1: Down miRNA targets that overlap with Up mRNA
if(!is.null(down_mirna_targets)) {
  pair1 <- down_mirna_targets %>%
    filter(rownames(.) %in% up_mrna_names) %>%
    mutate(
      mRNA = rownames(.),
      Pair_Type = "Up_mRNA_Down_miRNA"
    ) %>%
    select(miRNA, mRNA, Pair_Type, everything())
  
  cat("Pair 1 interactions found:", nrow(pair1), "\n")
} else {
  pair1 <- data.frame()
  cat("No Pair 1 interactions found\n")
}

# Pair 2: Up miRNA targets that overlap with Down mRNA
if(!is.null(up_mirna_targets)) {
  pair2 <- up_mirna_targets %>%
    filter(rownames(.) %in% down_mrna_names) %>%
    mutate(
      mRNA = rownames(.),
      Pair_Type = "Down_mRNA_Up_miRNA"
    ) %>%
    select(miRNA, mRNA, Pair_Type, everything())
  
  cat("Pair 2 interactions found:", nrow(pair2), "\n")
} else {
  pair2 <- data.frame()
  cat("No Pair 2 interactions found\n")
}

# Combine all interactions
all_interactions <- bind_rows(pair1, pair2)
cat("\nTotal mRNA-miRNA interactions:", nrow(all_interactions), "\n")

# ------------------------------------------------------------
# STEP 6: Visualise interaction network summary
# ------------------------------------------------------------

cat("\n--- Generating interaction summary plot ---\n")

if(nrow(all_interactions) > 0) {
  
  # Count interactions per miRNA
  mirna_interaction_counts <- all_interactions %>%
    group_by(miRNA, Pair_Type) %>%
    summarise(n_targets = n(), .groups = "drop") %>%
    arrange(desc(n_targets))
  
  png("results/figures/03_mirna_interaction_counts.png",
      width = 1400, height = 900, res = 120)
  
  ggplot(head(mirna_interaction_counts, 20), 
         aes(x = reorder(miRNA, n_targets), y = n_targets, fill = Pair_Type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c(
      "Up_mRNA_Down_miRNA"   = "#377EB8",
      "Down_mRNA_Up_miRNA"   = "#E41A1C"
    )) +
    coord_flip() +
    labs(
      title = "Top miRNA-mRNA Interactions in Psoriasis ceRNA Network",
      subtitle = "Based on opposite expression pairs",
      x = "miRNA",
      y = "Number of Target mRNAs",
      fill = "Interaction Type"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  dev.off()
  cat("✓ Interaction plot saved\n")
}

# ------------------------------------------------------------
# STEP 7: Save results for Cytoscape network construction
# ------------------------------------------------------------

cat("\n--- Saving results for Cytoscape ---\n")

# Save all interactions
write.csv(all_interactions, 
          "results/target_prediction/mrna_mirna_interactions.csv",
          row.names = FALSE)

# Save Pair 1 (Up mRNA + Down miRNA)
write.csv(pair1,
          "results/target_prediction/pair1_up_mrna_down_mirna.csv",
          row.names = FALSE)

# Save Pair 2 (Down mRNA + Up miRNA)
write.csv(pair2,
          "results/target_prediction/pair2_down_mrna_up_mirna.csv",
          row.names = FALSE)

# Save edge list for Cytoscape import
if(nrow(all_interactions) > 0) {
  cytoscape_edges <- all_interactions %>%
    select(Source = miRNA, Target = mRNA, Interaction = Pair_Type)
  
  write.csv(cytoscape_edges,
            "results/target_prediction/cytoscape_edge_list.csv",
            row.names = FALSE)
  
  cat("✓ Cytoscape edge list saved\n")
}

cat("\n============================================================\n")
cat("miRNA TARGET PREDICTION COMPLETE\n")
cat("============================================================\n")
cat("Results saved to: results/target_prediction/\n")
cat("Next step: Import cytoscape_edge_list.csv into Cytoscape\n")
cat("           for ceRNA network construction\n")
cat("============================================================\n")
