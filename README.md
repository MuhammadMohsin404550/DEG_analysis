# Differential Gene Expression Analysis in NSCLC Subtypes
### Adenocarcinoma (LUAD) vs Squamous Cell Carcinoma (LUSC) using TCGA RNA-seq Data

---

## Overview

This project performs differential gene expression (DGE) analysis between the two major subtypes of Non-Small Cell Lung Cancer (NSCLC) — lung adenocarcinoma (LUAD) and lung squamous cell carcinoma (LUSC) — using publicly available RNA-seq data from The Cancer Genome Atlas (TCGA).

The analysis identifies genes that are significantly upregulated or downregulated between the two subtypes, visualizes the results, and maps significant genes onto biological pathways using Gene Ontology (GO) enrichment analysis.

---

## Background

NSCLC accounts for approximately 85% of all lung cancer cases. Despite sharing the same organ of origin, LUAD and LUSC arise from different cell types, have distinct molecular profiles, and respond differently to clinical therapies. LUAD originates from glandular alveolar cells and is frequently associated with EGFR and KRAS mutations. LUSC originates from squamous bronchial epithelial cells and is strongly linked to smoking and TP63/KRT5 expression.

This analysis quantifies those transcriptional differences at the genome-wide level.

---

## Data Source

- **Database:** The Cancer Genome Atlas (TCGA) via the GDC Data Portal
- **Projects:** TCGA-LUAD and TCGA-LUSC
- **Data type:** RNA-seq Gene Expression Quantification (STAR - Counts, unstranded)
- **Genome reference:** hg38
- **Sample type:** Primary tumor samples only (TCGA barcode position 14-15 = "01")

| Dataset   | Samples (approx) |
|-----------|-----------------|
| TCGA-LUAD | ~500 tumor samples |
| TCGA-LUSC | ~500 tumor samples |

---

## Tools and R Packages

| Package | Version | Purpose |
|---|---|---|
| DESeq2 | Bioconductor | Differential expression statistical testing |
| TCGAbiolinks | Bioconductor | TCGA data download and preparation |
| SummarizedExperiment | Bioconductor | Data container format |
| apeglm | Bioconductor | Log fold change shrinkage |
| clusterProfiler | Bioconductor | GO pathway enrichment analysis |
| org.Hs.eg.db | Bioconductor | Human genome annotation |
| AnnotationDbi | Bioconductor | Gene ID conversion |
| EnhancedVolcano | Bioconductor | Volcano plot visualization |
| pheatmap | CRAN | Heatmap visualization |
| ggplot2 | CRAN | PCA and general plotting |
| dplyr | CRAN | Data manipulation |
| tibble | CRAN | Dataframe utilities |
| RColorBrewer | CRAN | Color palettes |

---

## Installation

Open R or RStudio and run the following to install all required packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "TCGAbiolinks", "SummarizedExperiment",
  "apeglm", "EnhancedVolcano", "clusterProfiler",
  "org.Hs.eg.db", "AnnotationDbi"
))

install.packages(c("ggplot2", "dplyr", "tibble", "pheatmap", "RColorBrewer", "ashr"))
```

---

## Project Structure

```
├── NSCLC_DEG_analysis.R          # Main analysis script (full pipeline)
├── README.md                     # This file
├── outputs/
│   ├── volcano_LUSC_vs_LUAD.png         # Volcano plot
│   ├── heatmap_top50_DEGs.png           # Heatmap of top 50 DEGs
│   ├── PCA_LUAD_vs_LUSC.png             # PCA plot
│   ├── GO_enrichment_LUSC.png           # GO dotplot for LUSC-high genes
│   ├── GO_enrichment_LUAD.png           # GO dotplot for LUAD-high genes
│   ├── NSCLC_significant_DEGs.csv       # All significant DEGs
│   ├── genes_higher_in_LUSC.csv         # Genes upregulated in LUSC
│   ├── genes_higher_in_LUAD.csv         # Genes upregulated in LUAD
│   ├── GO_enrichment_LUSC.csv           # GO enrichment table for LUSC
│   ├── GO_enrichment_LUAD.csv           # GO enrichment table for LUAD
│   └── dds_NSCLC_LUAD_LUSC.rds         # Saved DESeq2 object
```

---

## Pipeline Workflow

### Stage 1 — Data Download
RNA-seq count data is downloaded from TCGA using `TCGAbiolinks`. Only primary tumor samples are retained. Both LUAD and LUSC count matrices are merged on common genes.

### Stage 2 — Pre-filtering
Genes with fewer than 10 counts in fewer than 10 samples are removed to reduce noise and multiple testing burden.

### Stage 3 — DESeq2 Analysis
A DESeq2 object is constructed with the design formula `~ subtype` (LUAD as reference). The pipeline performs:
- Size factor estimation (sequencing depth normalization)
- Dispersion estimation (gene-level biological variability)
- Negative binomial GLM fitting and Wald testing per gene
- Log fold change shrinkage using `apeglm`

### Stage 4 — Filtering Significant Genes
Significance thresholds applied:
- Adjusted p-value (Benjamini-Hochberg FDR) < 0.05
- Absolute log2 fold change > 1 (minimum 2-fold difference)

### Stage 5 — Visualization
Four plots are generated and saved as high-resolution PNG files:
- Volcano plot with top 30 labeled genes
- Heatmap of top 50 most significant DEGs
- PCA plot of all samples colored by subtype
- GO Biological Process enrichment dotplots for each direction

### Stage 6 — Pathway Enrichment
Significant genes are mapped to GO Biological Process terms using `clusterProfiler::enrichGO()` with BH-adjusted p-value cutoff of 0.05.

---

## How to Run

Clone this repository and open `NSCLC_DEG_analysis.R` in RStudio. Run the script from top to bottom. The data download steps (Stages 2 and 3) require an internet connection and will take significant time — approximately 2.5 GB of data per project.

If downloads fail due to connection issues, reduce the chunk size:

```r
GDCdownload(query_luad, method = "api", files.per.chunk = 5)
```

All output files will be saved to your working directory. Check it with:

```r
getwd()
```

---

## Key Results Expected

Genes consistently found higher in LUSC include TP63, KRT5, KRT6A, and SOX2 — markers of squamous epithelial identity. Genes consistently found higher in LUAD include NKX2-1 (TTF1), NAPSA, SFTPB, and EGFR — markers of alveolar glandular identity.

GO enrichment of LUSC-high genes typically shows terms related to cornification, keratinocyte differentiation, and epidermis development. GO enrichment of LUAD-high genes typically shows terms related to surfactant metabolism, glandular epithelium development, and mucus secretion.

Graph:
<img width="4000" height="3200" alt="volcano_LUSC_vs_LUAD" src="https://github.com/user-attachments/assets/6e75378b-c869-416c-8c75-b99e2c2e9196" />

---

## Statistical Notes

- Multiple testing correction: Benjamini-Hochberg FDR
- Fold change estimation: apeglm shrinkage (more reliable for low-count genes)
- Normalization: DESeq2 median-of-ratios method
- Total variables tested: approximately 39,000 genes before filtering

---

## Limitations

- This analysis does not correct for potential batch effects between LUAD and LUSC cohorts
- Normal adjacent tissue samples are excluded — this is a tumor-vs-tumor comparison only
- Results represent population-level trends across hundreds of patients and may not apply to individual tumor biology
- Some ENSEMBL IDs may not map to gene symbols in the current annotation database version

---

## References

- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* (2014)
- Colaprico A et al. TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. *Nucleic Acids Research* (2016)
- Yu G et al. clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters. *OMICS* (2012)
- The Cancer Genome Atlas Research Network. Comprehensive molecular profiling of lung adenocarcinoma. *Nature* (2014)
- The Cancer Genome Atlas Research Network. Comprehensive genomic characterization of squamous cell lung cancers. *Nature* (2012)

---

## License

This project is for educational and research purposes. TCGA data is publicly available under the NIH TCGA data use policy.

---

## One thing to remember
This is can help you to understand basic, and also that the differentail gene expresion analysis pipeline changes based on what you are working. The core statistical engine — whether DESeq2 or edgeR — stays the same, but almost every other decision in the pipeline changes depending on your biological question, your data type, and your experimental context. There is no single universal DGE pipeline. 

---

## Author

Feel free to open an issue or submit a pull request if you find bugs or want to contribute improvements to the pipeline.
