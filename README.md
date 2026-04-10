# Comprehensive Analysis of lncRNA and circRNA Mediated ceRNA Network in Psoriasis

**Author:** Saima Imran  
**Institution:** University of Skövde, Sweden  
**Degree:** M.Sc. in Informatics — Bioinformatics (60 ECTS)  
**Supervisor:** Zelmina Lubovac  
**Dataset:** GSE145305 (NCBI GEO)  

---

## Project Overview

This repository contains the complete bioinformatics analysis pipeline developed for my Master's thesis in Bioinformatics. The project investigates the role of **competing endogenous RNA (ceRNA) networks** in psoriasis by analysing lncRNA and circRNA mediated regulatory mechanisms using multi-omics microarray data.

Psoriasis is a chronic autoimmune skin disease. The ceRNA hypothesis suggests that lncRNAs and circRNAs can act as molecular sponges for miRNAs, thereby regulating mRNA expression. This project computationally identifies key ceRNA interactions that may serve as therapeutic targets for psoriasis.

---

## Dataset

| Property | Details |
|----------|---------|
| **GEO Accession** | GSE145305 |
| **Platform (mRNA)** | GPL17586 — Affymetrix HTA-2_0 Array |
| **Platform (miRNA)** | GPL19117 — Affymetrix miRNA-4 Array |
| **mRNA Samples** | 6 (3 Psoriasis + 3 Healthy Control) |
| **miRNA Samples** | 8 (4 Psoriasis + 4 Healthy Control) |
| **Data Access** | Publicly available — FAIR compliant |

---

## Analysis Pipeline

```
GSE145305 (NCBI GEO)
        │
        ▼
01_data_preprocessing.R
   ├── CEL file reading (oligo)
   ├── RMA normalisation
   └── Quality control boxplots
        │
        ▼
02_differential_expression.R
   ├── Linear modelling (limma)
   ├── Empirical Bayes moderation
   ├── Volcano plots
   └── DEG/DEM identification
        │
        ▼
03_mirna_target_prediction.R
   ├── Opposite expression pairing
   ├── miRNAtap target prediction
   └── Cytoscape edge list export
        │
        ▼
Cytoscape (external)
   ├── ceRNA network construction
   ├── Hub gene identification (cytoHubba)
   └── PPI network (STRING plugin)
        │
        ▼
Biological Interpretation
   ├── Hub lncRNAs: XIST, NEAT1, KCNQ1OT1, MALAT1
   ├── Hub circRNAs: DSC2, GAN, HEPHL1, IFI44L
   └── Pathway analysis (clueGO)
```

---

## Key Results

### Differentially Expressed mRNA
| Regulation | Count | Threshold |
|-----------|-------|-----------|
| Up-regulated | 86 genes | logFC > 2, p < 0.05 |
| Down-regulated | 10 genes | logFC < -2, p < 0.05 |

**Top up-regulated genes:** SERPINB4, SERPINB3, PI3, S100A9, DEFB4B  
**Top down-regulated genes:** XIST, FLG2, FLG, CDR1, DUSP1

### Differentially Expressed miRNA
| Regulation | Count | Threshold |
|-----------|-------|-----------|
| Up-regulated | 75 miRNAs | logFC > 1.5, p < 0.05 |
| Down-regulated | 48 miRNAs | logFC < -1.5, p < 0.05 |

### Hub Genes Identified
| Type | Hub Genes |
|------|----------|
| Hub lncRNA | XIST, NEAT1, KCNQ1OT1, MALAT1, HCG18, NORAD |
| Hub circRNA | DSC2, GAN, HEPHL1, IFI44L, HK1, STAT1 |
| Hub PPI module | CCNB1, DLGAP5, ANLN, KPNA2, STAT1, OAS2, IFI44L, EIF2AK2 |

---

## Tools and Packages

### R Packages
| Package | Purpose | Version |
|---------|---------|---------|
| `GEOquery` | Download GEO datasets | Bioconductor |
| `oligo` | Read CEL files | Bioconductor |
| `limma` | Differential expression | Bioconductor |
| `miRNAtap` | miRNA target prediction | Bioconductor |
| `EnhancedVolcano` | Volcano plots | Bioconductor |
| `clueGO` | Pathway analysis | Cytoscape plugin |

### External Tools
| Tool | Purpose |
|------|---------|
| Cytoscape | Network visualisation |
| cytoHubba | Hub node identification |
| STRING plugin | PPI network |
| MCODE plugin | Module detection |

### Databases Used
- NCBI GEO — raw microarray data
- miRDip v4.1 — miRNA target predictions
- StarBase — miRNA-lncRNA interactions
- CircBase — circRNA database
- LncBase — lncRNA-miRNA interactions
- STRING — protein-protein interactions

---

## Repository Structure

```
psoriasis-ceRNA-network/
│
├── 01_data_preprocessing.R      # GEO download + RMA normalisation
├── 02_differential_expression.R # limma DEG analysis + volcano plots
├── 03_mirna_target_prediction.R # miRNAtap target prediction
│
├── results/
│   ├── figures/                 # All generated plots
│   └── data/                    # Processed data files
│       └── target_prediction/   # miRNA-mRNA interaction tables
│
└── README.md
```

---

## How to Run

### Requirements
- R version ≥ 4.0
- RStudio (recommended)
- Internet connection (for GEO download)

### Installation
```r
# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install all required packages
BiocManager::install(c(
  "GEOquery", "oligo", "limma", 
  "miRNAtap", "miRNAtap.db",
  "EnhancedVolcano", "Biobase"
))

install.packages(c("tidyverse", "ggplot2", "ggrepel"))
```

### Execution Order
```bash
# Run scripts in order
Rscript 01_data_preprocessing.R
Rscript 02_differential_expression.R
Rscript 03_mirna_target_prediction.R
```

---

## FAIR Data Principles

This project follows FAIR (Findable, Accessible, Interoperable, Reusable) data principles:

- **Findable** — Data sourced from NCBI GEO with persistent accession number (GSE145305)
- **Accessible** — All raw data is publicly available with no access restrictions
- **Interoperable** — Standard file formats (CEL, CSV) and documented workflows
- **Reusable** — Fully reproducible pipeline with documented methods and thresholds

---

## Contact

**Saima Imran**  
M.Sc. Bioinformatics | Data Coordinator | Precision Medicine  
📧 saimaimran4822@hotmail.com  
🔗 GitHub: [saima-imran](https://github.com/saima-imran)  
📍 Gothenburg, Sweden  

---

## Citation

If you use this pipeline in your research, please cite:

> Imran, S. (2022). *Comprehensive Analysis of lncRNA and circRNA Mediated ceRNA Network in Psoriasis*. Master's thesis, University of Skövde, Sweden.
