# NSCLC DEG Analysis — Line by Line Code Explanation

This document explains every line of the `NSCLC_DEG_analysis.R` pipeline in plain language — what each command does technically and why it is needed.

---

## Stage 1: Load Libraries

```r
library(DESeq2)
```
Loads the main statistical engine of the project. DESeq2 handles all the normalization and differential expression testing on RNA-seq count data.

---

```r
library(TCGAbiolinks)
```
Loads the tool that connects to the TCGA/GDC database online and lets you search, download, and prepare cancer genomics data directly inside R.

---

```r
library(SummarizedExperiment)
```
Loads the container format that TCGA data comes packaged in. It holds the count matrix, sample information, and gene information all together in one object.

---

```r
library(EnhancedVolcano)
```
Loads the package that draws the volcano plot — a specific type of scatter plot designed for showing differential expression results with labels and color coding.

---

```r
library(clusterProfiler)
```
Loads the pathway enrichment tool. After you find your significant genes, this package tests which biological pathways those genes belong to.

---

```r
library(org.Hs.eg.db)
```
Loads the human genome annotation database. This is a lookup table that converts between different gene ID formats — ENSEMBL IDs, gene symbols, Entrez IDs, and so on.

---

```r
library(AnnotationDbi)
```
Loads the interface that lets you query the annotation database above using functions like `mapIds()`.

---

```r
library(pheatmap)
```
Loads the heatmap drawing package. "pretty heatmap" — it draws clustered heatmaps with annotation bars along the top.

---

```r
library(ggplot2)
```
Loads the general plotting system used for the PCA plot and as the underlying engine for several other visualizations.

---

```r
library(dplyr)
```
Loads data manipulation tools. Functions like `filter()`, `arrange()`, `head()`, and `pull()` all come from here. Used throughout to subset and sort your results tables.

---

```r
library(tibble)
```
Loads tools for working with dataframes. Specifically used here for `rownames_to_column()` which moves rownames into an actual column so gene IDs are not lost when manipulating the results table.

---

```r
library(RColorBrewer)
```
Loads color palette tools used for making plots visually cleaner and publication-ready.

---

```r
library(apeglm)
```
Loads the statistical method used for log fold change shrinkage. Without this package the `lfcShrink()` call in Stage 7 will throw an error.

---

## Stage 2: Download LUAD Data

```r
query_luad <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
```
Sends a search request to the GDC database. You are asking: give me all RNA-seq files from the lung adenocarcinoma project, specifically the gene-level count files that were generated using the STAR aligner. This does not download anything yet — it just builds a manifest of available files and stores it in `query_luad`.

---

```r
GDCdownload(query_luad, method = "api", files.per.chunk = 10)
```
Actually downloads the files to your computer. `method = "api"` uses the GDC REST API to fetch files. `files.per.chunk = 10` means it downloads 10 files at a time instead of hundreds simultaneously, which prevents connection timeouts on slow or unstable internet connections.

---

```r
luad_data <- GDCprepare(query_luad)
```
Reads all the downloaded files from disk and assembles them into a single SummarizedExperiment object. All the individual count files get merged into one count matrix with genes as rows and samples as columns.

---

## Stage 3: Download LUSC Data

```r
query_lusc <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
```
Same search as above but for the squamous cell carcinoma project instead of adenocarcinoma.

---

```r
GDCdownload(query_lusc, method = "api", files.per.chunk = 10)
```
Downloads the LUSC files in chunks of 10.

---

```r
lusc_data <- GDCprepare(query_lusc)
```
Assembles the LUSC files into a SummarizedExperiment object.

---

## Stage 4: Build Combined Count Matrix

```r
luad_counts <- assay(luad_data, "unstranded")
```
Pulls the raw count matrix out of the LUAD SummarizedExperiment object. `"unstranded"` refers to the specific assay slot containing the unstranded read counts — the total reads mapping to each gene regardless of which DNA strand they came from.

---

```r
lusc_counts <- assay(lusc_data, "unstranded")
```
Same extraction for LUSC.

---

```r
luad_counts <- luad_counts[, substr(colnames(luad_counts), 14, 15) == "01"]
```
Filters columns to keep only primary tumor samples. TCGA barcodes are formatted strings like `TCGA-05-4244-01A-01R-1107-07`. Characters at positions 14 and 15 encode the sample type — `01` means primary tumor, `11` means normal tissue, `06` means metastatic. This line keeps only the tumor samples.

---

```r
lusc_counts <- lusc_counts[, substr(colnames(lusc_counts), 14, 15) == "01"]
```
Same tumor filtering for LUSC.

---

```r
common_genes <- intersect(rownames(luad_counts), rownames(lusc_counts))
```
Finds the genes present in both datasets. Although both come from the same TCGA pipeline, there can occasionally be slight differences in which genes are included. `intersect()` returns only the genes found in both, so the two matrices can be safely merged.

---

```r
luad_counts <- luad_counts[common_genes, ]
lusc_counts <- lusc_counts[common_genes, ]
```
Subsets both matrices to only the common genes, and in the same order, so rows align correctly when merged.

---

```r
combined_counts <- cbind(luad_counts, lusc_counts)
```
Merges the two matrices side by side. `cbind` means column-bind — LUAD sample columns come first, then LUSC sample columns, all sharing the same gene rows.

---

```r
metadata <- data.frame(
  sample  = colnames(combined_counts),
  subtype = factor(c(rep("LUAD", ncol(luad_counts)),
                     rep("LUSC", ncol(lusc_counts)))),
  row.names = colnames(combined_counts)
)
```
Creates the sample information table. Each row is one sample. The `sample` column holds the TCGA barcode. The `subtype` column labels each sample as either LUAD or LUSC — as many LUAD labels as there are LUAD columns, followed by as many LUSC labels as there are LUSC columns. `factor()` converts the text labels into a categorical variable that DESeq2 understands. The row names are set to the sample barcodes so the metadata rows match the count matrix columns exactly.

---

```r
cat("LUAD samples:", ncol(luad_counts), "\n")
cat("LUSC samples:", ncol(lusc_counts), "\n")
cat("Total genes:", nrow(combined_counts), "\n")
```
Prints three numbers to the console so you can confirm how many samples and genes you are working with before proceeding.

---

## Stage 5: Filter Low-Count Genes

```r
keep <- rowSums(combined_counts >= 10) >= 10
```
For each gene (row), this counts how many samples have at least 10 reads mapping to that gene. If that count is 10 or more samples, the gene passes the filter. The result `keep` is a TRUE/FALSE vector with one value per gene.

---

```r
combined_counts_filtered <- combined_counts[keep, ]
```
Applies the filter, keeping only rows where `keep` is TRUE. Genes that are barely detected across samples are removed.

---

```r
cat("Genes before filtering:", nrow(combined_counts), "\n")
cat("Genes after filtering:",  nrow(combined_counts_filtered), "\n")
```
Prints gene counts before and after filtering so you can see how many low-count genes were removed.

---

## Stage 6: Build DESeq2 Object

```r
dds <- DESeqDataSetFromMatrix(
  countData = combined_counts_filtered,
  colData   = metadata,
  design    = ~ subtype
)
```
Packages everything into the DESeq2 format. `countData` is your filtered count matrix. `colData` is your metadata table. `design = ~ subtype` is the statistical formula telling DESeq2 what comparison to make — it means model gene expression as a function of subtype (LUAD or LUSC).

---

```r
dds$subtype <- relevel(dds$subtype, ref = "LUAD")
```
Sets LUAD as the reference group. This means all fold changes will be calculated as LUSC relative to LUAD. A positive fold change means higher in LUSC; negative means higher in LUAD.

---

```r
cat("DESeq2 object:", nrow(dds), "genes,", ncol(dds), "samples\n")
```
Confirms the dimensions of your DESeq2 object.

---

## Stage 7: Run DESeq2

```r
dds <- DESeq(dds)
```
Runs the entire DESeq2 statistical pipeline in one call. Internally it does three things sequentially: estimates size factors to normalize for sequencing depth differences between samples, estimates dispersion (gene-by-gene biological variability), and fits a negative binomial generalized linear model and performs Wald tests for every gene. This is the most computationally intensive step and may take 10 to 30 minutes.

---

```r
cat("Coefficient names:\n")
print(resultsNames(dds))
```
Prints the names DESeq2 assigned to the model coefficients. You need to know the exact name of the LUSC vs LUAD coefficient to use in the next line. It will typically print `subtype_LUSC_vs_LUAD`.

---

```r
results_raw <- results(dds,
                       contrast = c("subtype", "LUSC", "LUAD"),
                       alpha = 0.05)
```
Extracts the raw results table from the fitted model. `contrast` explicitly specifies the comparison direction — LUSC versus LUAD. `alpha = 0.05` sets the FDR threshold used internally for flagging significant genes in the summary output.

---

```r
results_shrunk <- lfcShrink(dds,
                             coef = "subtype_LUSC_vs_LUAD",
                             type = "apeglm")
```
Applies log fold change shrinkage using the apeglm method. Genes with low counts have unreliable fold change estimates that can be spuriously large. Shrinkage pulls those estimates toward zero, making the fold change column more trustworthy for ranking and visualization without changing which genes are statistically significant.

---

```r
summary(results_raw)
```
Prints a summary showing how many genes are upregulated, downregulated, and not significant at your chosen threshold.

---

## Stage 8: Filter Significant Genes

```r
results_df <- as.data.frame(results_shrunk) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj)
```
Converts the results object into a regular dataframe. `rownames_to_column("gene_id")` moves the ENSEMBL gene IDs from rownames into an actual column called `gene_id` so they are not lost during subsequent operations. `arrange(padj)` sorts rows by adjusted p-value, most significant first.

---

```r
results_df$ensembl_clean <- gsub("\\..*", "", results_df$gene_id)
```
Removes version numbers from ENSEMBL IDs. TCGA data uses versioned IDs like `ENSG00000123456.5` — the `.5` is a version suffix that the annotation database does not recognize. `gsub("\\..*", "")` removes the dot and everything after it, leaving just `ENSG00000123456`.

---

```r
results_df$gene_symbol <- mapIds(org.Hs.eg.db,
                                  keys      = results_df$ensembl_clean,
                                  column    = "SYMBOL",
                                  keytype   = "ENSEMBL",
                                  multiVals = "first")
```
Looks up the human-readable gene symbol for each ENSEMBL ID in the annotation database. `keys` is your list of ENSEMBL IDs to look up. `column = "SYMBOL"` says return the gene symbol. `keytype = "ENSEMBL"` says the input format is ENSEMBL. `multiVals = "first"` says if a gene maps to multiple symbols, just take the first one.

---

```r
results_df$gene_symbol[is.na(results_df$gene_symbol)] <-
  results_df$ensembl_clean[is.na(results_df$gene_symbol)]
```
For any gene that did not get a symbol from the database lookup, substitutes the ENSEMBL ID as a fallback so the column is never empty.

---

```r
sig_genes <- results_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)
```
Filters to genes that pass both thresholds simultaneously — adjusted p-value below 0.05 (statistically significant after multiple testing correction) and absolute log2 fold change above 1 (meaning at least a 2-fold difference in expression between the subtypes).

---

```r
lusc_up <- sig_genes %>% filter(log2FoldChange >  1) %>% arrange(desc(log2FoldChange))
luad_up <- sig_genes %>% filter(log2FoldChange < -1) %>% arrange(log2FoldChange)
```
Splits the significant genes into two separate lists. `lusc_up` contains genes with higher expression in LUSC (positive fold change). `luad_up` contains genes with higher expression in LUAD (negative fold change, because LUAD is the reference).

---

```r
cat("Total significant DEGs:", nrow(sig_genes), "\n")
cat("Higher in LUSC:", nrow(lusc_up), "\n")
cat("Higher in LUAD:", nrow(luad_up), "\n")
```
Prints the counts of significant genes in each category.

---

```r
head(lusc_up[, c("gene_symbol", "log2FoldChange", "padj")], 20)
head(luad_up[, c("gene_symbol", "log2FoldChange", "padj")], 20)
```
Prints the top 20 genes from each list to the console so you can immediately see which genes are most strongly distinguishing the two subtypes.

---

## Stage 9: Volcano Plot

```r
top_lusc_labels <- results_df %>%
  filter(log2FoldChange > 1, padj < 0.05) %>%
  arrange(padj) %>%
  head(15) %>%
  pull(gene_symbol)
```
Selects the 15 most statistically significant LUSC-high genes to label on the plot. `pull(gene_symbol)` extracts just the gene symbol column as a plain vector.

---

```r
top_luad_labels <- results_df %>%
  filter(log2FoldChange < -1, padj < 0.05) %>%
  arrange(padj) %>%
  head(15) %>%
  pull(gene_symbol)
```
Same but for the 15 most significant LUAD-high genes.

---

```r
genes_to_label <- c(top_lusc_labels, top_luad_labels)
```
Combines both label lists into one vector of 30 gene names that will be labeled on the volcano plot.

---

```r
png("volcano_LUSC_vs_LUAD.png", width = 4000, height = 3200, res = 300)
```
Opens a PNG file for writing. All plotting commands until `dev.off()` will be written into this file instead of the screen. Width and height are in pixels. `res = 300` means 300 dots per inch which gives publication quality resolution.

---

```r
EnhancedVolcano(results_df, ...)
```
Draws the volcano plot. Each dot is one gene. The x-axis is log2 fold change — how much higher or lower a gene is in LUSC versus LUAD. The y-axis is negative log10 of the adjusted p-value — higher means more statistically significant. `selectLab` means only label the 30 genes you specified. `drawConnectors = TRUE` draws lines from labels to their dots when labels are moved to avoid overlap. `xlim` and `ylim` set the axis ranges. `col` sets the colors for the four categories of points.

---

```r
dev.off()
```
Closes the PNG file and finalizes writing it to disk. Always required after `png()`.

---

## Stage 10: Heatmap

```r
normalized_counts <- counts(dds, normalized = TRUE)
```
Extracts the size-factor normalized count matrix from the DESeq2 object. These counts have been corrected for sequencing depth differences between samples, making them comparable across samples.

---

```r
top50_genes <- results_df %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene_id)
```
Selects the 50 genes with the smallest adjusted p-values — the most statistically confident differential expression signals.

---

```r
top50_genes <- top50_genes[top50_genes %in% rownames(normalized_counts)]
```
Removes any gene IDs from the list that are not present as rownames in the normalized count matrix, preventing a subscript error when subsetting.

---

```r
top50_matrix <- normalized_counts[top50_genes, ]
```
Subsets the normalized count matrix to only the 50 selected genes.

---

```r
top50_log <- log2(top50_matrix + 1)
```
Log2 transforms the counts. The `+ 1` is added before taking the log to avoid log(0) which is undefined — this is called a pseudocount. Log transformation compresses the dynamic range so that very highly expressed genes do not visually dominate the heatmap.

---

```r
top50_scaled <- t(scale(t(top50_log)))
```
Z-score scales each gene across all samples. For each gene, it subtracts the mean expression and divides by the standard deviation. The double transpose `t()` is needed because `scale()` works on columns but you want to scale across rows (genes). After scaling, a value of 0 means average expression for that gene, positive means above average, negative means below average. This makes visual comparison across genes meaningful even when their absolute expression levels differ greatly.

---

```r
top50_scaled[is.nan(top50_scaled)] <- 0
```
Replaces any NaN values with 0. NaN arises when scaling a gene that has identical expression in every sample — dividing by a standard deviation of zero produces NaN. Setting these to 0 prevents the heatmap from failing.

---

```r
row_labels <- results_df$gene_symbol[match(top50_genes, results_df$gene_id)]
```
Looks up the gene symbol for each of the 50 genes so the heatmap shows readable names instead of ENSEMBL IDs. `match()` finds the position of each gene ID in `results_df` and retrieves the corresponding symbol.

---

```r
row_labels[is.na(row_labels)] <- top50_genes[is.na(row_labels)]
```
For any gene symbol that is still NA after the lookup, substitutes the ENSEMBL ID as a fallback row label.

---

```r
annotation_col <- data.frame(
  Subtype   = metadata$subtype,
  row.names = rownames(metadata)
)
```
Creates the annotation dataframe that pheatmap uses to draw the colored bar along the top of the heatmap showing which samples are LUAD and which are LUSC.

---

```r
annotation_colors <- list(
  Subtype = c(LUAD = "steelblue", LUSC = "firebrick")
)
```
Assigns specific colors to the annotation bar — steelblue for LUAD and firebrick red for LUSC.

---

```r
png("heatmap_top50_DEGs.png", width = 4000, height = 4500, res = 300)
```
Opens the PNG file for the heatmap. Taller than the volcano plot to accommodate the 50 gene rows.

---

```r
pheatmap(top50_scaled,
         labels_row        = row_labels,
         annotation_col    = annotation_col,
         annotation_colors = annotation_colors,
         show_colnames     = FALSE,
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         color             = colorRampPalette(c("navy", "white", "firebrick"))(100),
         main              = "Top 50 Differentially Expressed Genes",
         fontsize_row      = 9,
         fontsize          = 11,
         border_color      = NA)
```
Draws the heatmap. `labels_row` uses the gene symbols you prepared. `show_colnames = FALSE` hides sample barcodes since there are hundreds of samples and the names would be unreadable. `cluster_rows = TRUE` and `cluster_cols = TRUE` both apply hierarchical clustering — genes with similar patterns group together, and samples with similar profiles group together. `colorRampPalette(c("navy", "white", "firebrick"))(100)` creates a 100-color gradient from navy (low expression) through white (average) to firebrick (high expression). `border_color = NA` removes the grid lines between cells for a cleaner look.

---

```r
dev.off()
```
Closes and saves the heatmap PNG file.

---

## Stage 11: PCA Plot

```r
vsd <- vst(dds, blind = FALSE)
```
Applies variance stabilizing transformation to the count data. Raw counts and even log-transformed counts have a mean-variance relationship where highly expressed genes also have higher variance. VST removes this relationship, making the data appropriate for PCA. `blind = FALSE` means the transformation uses knowledge of the experimental design, which is slightly more accurate than the blind version.

---

```r
pca_data <- plotPCA(vsd, intgroup = "subtype", returnData = TRUE)
```
Performs PCA on the top 500 most variable genes and returns the coordinates. `intgroup = "subtype"` tells it to color points by subtype. `returnData = TRUE` returns the PCA coordinates as a dataframe instead of immediately drawing a plot, so you can use ggplot2 for a more customized figure.

---

```r
pct_var <- round(100 * attr(pca_data, "percentVar"))
```
Extracts how much variance PC1 and PC2 each explain, rounded to whole numbers, for use in the axis labels.

---

```r
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = subtype)) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_color_manual(values = c(LUAD = "steelblue", LUSC = "firebrick")) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  ggtitle("PCA of TCGA LUAD vs LUSC Samples") +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank())
```
Builds the PCA plot using ggplot2. Each `+` adds a layer. `geom_point` draws the sample dots. `alpha = 0.7` makes dots slightly transparent so overlapping points are visible. `scale_color_manual` assigns the blue and red colors. `xlab` and `ylab` add axis labels with the variance percentages. `theme_classic` gives a clean white background with no grid lines. `theme(legend.title = element_blank())` removes the word "subtype" from the legend since it is obvious from context.

---

```r
ggsave("PCA_LUAD_vs_LUSC.png", plot = pca_plot, width = 10, height = 8, dpi = 300)
```
Saves the ggplot object to a PNG file. `ggsave` is the correct way to save ggplot objects — it is more reliable than wrapping with `png()` and `dev.off()` for ggplot2 plots.

---

## Stage 12: GO Pathway Enrichment

```r
lusc_entrez <- mapIds(org.Hs.eg.db,
                      keys      = gsub("\\..*", "", lusc_up$gene_id),
                      column    = "ENTREZID",
                      keytype   = "ENSEMBL",
                      multiVals = "first")
```
Converts your LUSC significant gene ENSEMBL IDs to Entrez IDs. `clusterProfiler` requires Entrez IDs rather than ENSEMBL IDs. `gsub("\\..*", "")` strips version numbers first.

---

```r
luad_entrez <- mapIds(org.Hs.eg.db, ...)
```
Same conversion for LUAD significant genes.

---

```r
lusc_entrez <- na.omit(lusc_entrez)
luad_entrez <- na.omit(luad_entrez)
```
Removes any genes that failed to map to an Entrez ID. Passing NA values into `enrichGO` would cause an error.

---

```r
ego_lusc <- enrichGO(gene          = lusc_entrez,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     readable      = TRUE)
```
Runs GO over-representation analysis on the LUSC-high genes. `ont = "BP"` tests Biological Process terms only (as opposed to Molecular Function or Cellular Component). `pAdjustMethod = "BH"` applies Benjamini-Hochberg FDR correction to the pathway p-values. `pvalueCutoff = 0.05` keeps only significant pathways. `readable = TRUE` converts Entrez IDs back to gene symbols in the output so results show gene names not numbers.

---

```r
ego_luad <- enrichGO(gene = luad_entrez, ...)
```
Same pathway analysis for LUAD-high genes.

---

```r
png("GO_enrichment_LUSC.png", width = 4000, height = 3500, res = 300)
print(dotplot(ego_lusc, showCategory = 20, title = "GO BP: Pathways enriched in LUSC"))
dev.off()
```
Saves the LUSC pathway dotplot. `dotplot` is a visualization where each GO term is a row, dot size represents how many of your significant genes belong to that pathway, and dot color represents the adjusted p-value. `print()` is required here because inside a `png()` block, ggplot-based plots do not render unless explicitly printed. `showCategory = 20` shows the top 20 enriched pathways.

---

```r
png("GO_enrichment_LUAD.png", width = 4000, height = 3500, res = 300)
print(dotplot(ego_luad, showCategory = 20, title = "GO BP: Pathways enriched in LUAD"))
dev.off()
```
Same for the LUAD pathway dotplot.

---

## Stage 13: Save All Result Tables

```r
write.csv(sig_genes, "NSCLC_significant_DEGs.csv", row.names = FALSE)
```
Saves all significant DEGs to a CSV file. `row.names = FALSE` prevents R from adding an extra numbered column at the start of the file.

---

```r
write.csv(lusc_up, "genes_higher_in_LUSC.csv", row.names = FALSE)
write.csv(luad_up, "genes_higher_in_LUAD.csv", row.names = FALSE)
```
Saves the two directional gene lists separately.

---

```r
write.csv(as.data.frame(ego_lusc), "GO_enrichment_LUSC.csv", row.names = FALSE)
write.csv(as.data.frame(ego_luad), "GO_enrichment_LUAD.csv", row.names = FALSE)
```
Saves the GO enrichment result tables. `as.data.frame()` converts the enrichResult object into a plain dataframe that `write.csv` can handle.

---

```r
saveRDS(dds, "dds_NSCLC_LUAD_LUSC.rds")
```
Saves the entire DESeq2 object to disk in R's binary format. This means next time you open R you can load it back with `readRDS("dds_NSCLC_LUAD_LUSC.rds")` and immediately access your normalized counts, fitted model, and results without rerunning the hours-long DESeq2 computation.

---

```r
cat("All files saved to:", getwd(), "\n")
```
Prints the folder path where all your output files have been saved so you know exactly where to find them.
