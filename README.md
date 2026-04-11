# Comprehensive Analysis of lncRNA and circRNA Mediated ceRNA Network in Psoriasis

**Author:** Saima Imran  
**Institution:** University of Skövde, Sweden  
**Degree:** M.Sc. in Informatics — Bioinformatics (60 ECTS)  
**Supervisor:** Zelmina Lubovac  
**Dataset:** GSE145305 (NCBI GEO)  
# Bioinformatics Portfolio — Saima Imran

**M.Sc. Bioinformatics | Data Coordinator | Precision Medicine**  
📧 saimaimran4822@hotmail.com | 📍 Gothenburg, Sweden  
🔗 [GitHub Profile](https://github.com/saima-imran)

---

## About This Portfolio

This repository contains two complete bioinformatics analysis projects developed using R and public genomics data from NCBI GEO. Both projects demonstrate end-to-end data science pipelines including data retrieval, preprocessing, statistical analysis, and publication-quality visualisation — all following FAIR data principles.

---

## Projects

### Project 1 — Psoriasis ceRNA Network Analysis
**Dataset:** GSE145305 | **Disease:** Psoriasis vs Healthy Skin

This project is based on my M.Sc. thesis in Bioinformatics at the University of Skövde, Sweden. It investigates the competing endogenous RNA (ceRNA) regulatory network in psoriasis using multi-omics microarray data.

**Analysis pipeline:**
```
NCBI GEO (GSE145305)
        ↓
Data preprocessing — oligo, RMA normalisation
        ↓
Differential expression — limma, empirical Bayes
        ↓
miRNA target prediction — miRNAtap, miRDip v4.1
        ↓
Network construction — Cytoscape, STRING, cytoHubba
        ↓
Pathway analysis — clueGO
```

**Key results:**

| RNA Type | Up-regulated | Down-regulated | Threshold |
|----------|-------------|----------------|-----------|
| mRNA | 86 genes | 10 genes | p<0.05, \|logFC\|>2 |
| miRNA | 75 miRNAs | 48 miRNAs | p<0.05, \|logFC\|>1.5 |

**Hub genes identified:**

| Type | Genes |
|------|-------|
| Hub lncRNA | XIST, NEAT1, KCNQ1OT1, MALAT1, HCG18, NORAD |
| Hub circRNA | DSC2, GAN, HEPHL1, IFI44L, HK1, STAT1 |
| Hub PPI module | CCNB1, DLGAP5, ANLN, KPNA2, STAT1, OAS2, IFI44L, EIF2AK2 |

**Figures generated:**

| Figure | Description |
|--------|-------------|
| `01_mRNA_boxplot.png` | Quality control — expression distribution across samples |
| `02_mRNA_volcano.png` | Volcano plot — significant DEGs in psoriasis |
| `03_heatmap_top30_DEGs.png` | Heatmap — top 30 differentially expressed genes |

---

### Project 2 — Breast Cancer Gene Expression Analysis
**Dataset:** GSE45827 | **Disease:** Breast Cancer Subtypes

This project analyses gene expression differences across breast cancer subtypes using public microarray data. The analysis focuses on Triple Negative Breast Cancer (Basal/TNBC) vs Normal tissue — directly relevant to precision oncology and personalised medicine.

**Cancer subtypes analysed:**

| Subtype | Samples | Colour |
|---------|---------|--------|
| Basal (TNBC) | 41 | 🔴 Red |
| HER2-enriched | 30 | 🟡 Orange |
| Luminal A | 29 | 🔵 Blue |
| Luminal B | 30 | 🟣 Purple |
| Normal | 11 | 🟢 Green |

**Key results — Basal (TNBC) vs Normal:**

| Direction | Count | Threshold |
|-----------|-------|-----------|
| Up in Basal | 4,383 genes | p<0.05, logFC>1.5 |
| Down in Basal | 1,967 genes | p<0.05, logFC<-1.5 |
| Total significant | 6,350 genes | — |

**Figures generated:**

| Figure | Description |
|--------|-------------|
| `01_BC_expression_boxplot.png` | Colourful boxplot — expression by cancer subtype |
| `02_BC_volcano_Basal_vs_Normal.png` | Volcano plot — TNBC vs Normal tissue |
| `03_BC_heatmap_top40_genes.png` | Heatmap — top 40 DEGs across all subtypes |

---

## Tools and Packages

### R Packages

| Package | Purpose |
|---------|---------|
| `GEOquery` | Download datasets from NCBI GEO |
| `oligo` | Read and preprocess CEL files |
| `limma` | Differential expression analysis |
| `pheatmap` | Professional heatmap visualisation |
| `ggplot2` | Data visualisation |
| `ggrepel` | Non-overlapping gene labels |
| `RColorBrewer` | Colour palettes |
| `tidyverse` | Data manipulation |

### External Tools

| Tool | Purpose |
|------|---------|
| Cytoscape | Network visualisation |
| cytoHubba | Hub gene identification |
| STRING plugin | Protein-protein interaction |
| MCODE plugin | Module detection |
| clueGO | Pathway enrichment analysis |

### Databases Used

| Database | Purpose |
|----------|---------|
| NCBI GEO | Raw microarray data |
| miRDip v4.1 | miRNA target predictions |
| StarBase | miRNA-lncRNA interactions |
| CircBase | circRNA database |
| LncBase | lncRNA-miRNA interactions |
| STRING | Protein-protein interactions |

---

## Repository Structure

```
bioinformatics-portfolio/
│
├── 01_psoriasis_analysis.R          # Complete psoriasis pipeline
├── 02_breast_cancer_analysis.R      # Complete breast cancer pipeline
│
├── results/
│   ├── figures/
│   │   ├── 01_mRNA_boxplot.png
│   │   ├── 02_mRNA_volcano.png
│   │   ├── 03_heatmap_top30_DEGs.png
│   │   └── breast_cancer/
│   │       ├── 01_BC_expression_boxplot.png
│   │       ├── 02_BC_volcano_Basal_vs_Normal.png
│   │       └── 03_BC_heatmap_top40_genes.png
│   └── data/
│       ├── mRNA_DEGs.csv
│       ├── miRNA_DEMs.csv
│       └── BreastCancer_DEGs_Basal_vs_Normal.csv
│
└── README.md
```

---

## How to Run

### Requirements
- R version ≥ 4.0
- RStudio (recommended)
- Internet connection

### Installation
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GEOquery", "oligo", "limma"))

install.packages(c("tidyverse", "ggplot2", "ggrepel", 
                   "pheatmap", "RColorBrewer"))
```

### Run Analysis
```r
# Project 1 — Psoriasis
source("01_psoriasis_analysis.R")

# Project 2 — Breast Cancer
source("02_breast_cancer_analysis.R")
```

---

## FAIR Data Principles

Both projects strictly follow FAIR data principles:

| Principle | Implementation |
|-----------|---------------|
| **Findable** | All data sourced from NCBI GEO with persistent accession numbers |
| **Accessible** | Raw data publicly available — no access restrictions |
| **Interoperable** | Standard file formats (CSV, PNG) and documented workflows |
| **Reusable** | Fully reproducible pipelines with documented methods and thresholds |

---

## Education & Background

| Degree | Institution | Year |
|--------|-------------|------|
| M.Sc. Bioinformatics (60 ECTS) | University of Skövde, Sweden | 2022 |
| M.Phil. Computer Science | University of Lahore, Pakistan | 2011 |
| Master of Computer Science | Bahauddin Zakariya University, Pakistan | 2008 |

**Certifications:**
- Microsoft Azure AI Fundamentals (AI-900)
- Microsoft Azure Cloud Fundamentals (AZ-900)
- Power BI Data Analyst (PL-300 — in progress)
- Microsoft Learn: Level 6 | 31,825 XP | 20 Badges

---

## Contact

**Saima Imran**  
M.Sc. Bioinformatics | Data Coordinator | Precision Medicine  
📧 saimaimran4822@hotmail.com  
📍 Gothenburg, Sweden  
🔗 [GitHub](https://github.com/saima-imran)

---

## References

**Project 1 — Psoriasis:**
> Imran, S. (2022). *Comprehensive Analysis of lncRNA and circRNA Mediated ceRNA Network in Psoriasis*. M.Sc. thesis, University of Skövde, Sweden.

**Project 2 — Breast Cancer:**
> Dataset: GSE45827. Breast cancer gene expression profiling. NCBI GEO. Available at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45827
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
