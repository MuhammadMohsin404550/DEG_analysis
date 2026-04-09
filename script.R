# ============================================================
# STAGE 1: Load Libraries
# ============================================================
library(DESeq2)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tibble)
library(RColorBrewer)
library(apeglm)


# ============================================================
# STAGE 2: Download LUAD Data
# ============================================================
query_luad <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_luad, method = "api", files.per.chunk = 10)
luad_data <- GDCprepare(query_luad)


# ============================================================
# STAGE 3: Download LUSC Data
# ============================================================
query_lusc <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_lusc, method = "api", files.per.chunk = 10)
lusc_data <- GDCprepare(query_lusc)


# ============================================================
# STAGE 4: Build Combined Count Matrix
# ============================================================
luad_counts <- assay(luad_data, "unstranded")
lusc_counts <- assay(lusc_data, "unstranded")

# Keep only primary tumor samples
luad_counts <- luad_counts[, substr(colnames(luad_counts), 14, 15) == "01"]
lusc_counts <- lusc_counts[, substr(colnames(lusc_counts), 14, 15) == "01"]

# Keep only common genes
common_genes <- intersect(rownames(luad_counts), rownames(lusc_counts))
luad_counts  <- luad_counts[common_genes, ]
lusc_counts  <- lusc_counts[common_genes, ]

# Merge
combined_counts <- cbind(luad_counts, lusc_counts)

# Metadata
metadata <- data.frame(
  sample  = colnames(combined_counts),
  subtype = factor(c(rep("LUAD", ncol(luad_counts)),
                     rep("LUSC", ncol(lusc_counts)))),
  row.names = colnames(combined_counts)
)

cat("LUAD samples:", ncol(luad_counts), "\n")
cat("LUSC samples:", ncol(lusc_counts), "\n")
cat("Total genes:", nrow(combined_counts), "\n")


# ============================================================
# STAGE 5: Filter Low-Count Genes
# ============================================================
keep <- rowSums(combined_counts >= 10) >= 10
combined_counts_filtered <- combined_counts[keep, ]

cat("Genes before filtering:", nrow(combined_counts), "\n")
cat("Genes after filtering:",  nrow(combined_counts_filtered), "\n")


# ============================================================
# STAGE 6: Build DESeq2 Object
# ============================================================
dds <- DESeqDataSetFromMatrix(
  countData = combined_counts_filtered,
  colData   = metadata,
  design    = ~ subtype
)

dds$subtype <- relevel(dds$subtype, ref = "LUAD")
cat("DESeq2 object:", nrow(dds), "genes,", ncol(dds), "samples\n")


# ============================================================
# STAGE 7: Run DESeq2
# ============================================================
dds <- DESeq(dds)

# Check coefficient name — must match exactly
cat("Coefficient names:\n")
print(resultsNames(dds))

results_raw <- results(dds,
                       contrast = c("subtype", "LUSC", "LUAD"),
                       alpha = 0.05)

results_shrunk <- lfcShrink(dds,
                            coef = "subtype_LUSC_vs_LUAD",
                            type = "apeglm")

summary(results_raw)


# ============================================================
# STAGE 8: Filter Significant Genes
# ============================================================
results_df <- as.data.frame(results_shrunk) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj)

# Clean ENSEMBL IDs and map to gene symbols
results_df$ensembl_clean <- gsub("\\..*", "", results_df$gene_id)

results_df$gene_symbol <- mapIds(org.Hs.eg.db,
                                 keys      = results_df$ensembl_clean,
                                 column    = "SYMBOL",
                                 keytype   = "ENSEMBL",
                                 multiVals = "first")

# Fall back to ENSEMBL ID if symbol not found
results_df$gene_symbol[is.na(results_df$gene_symbol)] <-
  results_df$ensembl_clean[is.na(results_df$gene_symbol)]

# Significant genes
sig_genes <- results_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

lusc_up <- sig_genes %>% filter(log2FoldChange >  1) %>% arrange(desc(log2FoldChange))
luad_up <- sig_genes %>% filter(log2FoldChange < -1) %>% arrange(log2FoldChange)

cat("Total significant DEGs:", nrow(sig_genes), "\n")
cat("Higher in LUSC:", nrow(lusc_up), "\n")
cat("Higher in LUAD:", nrow(luad_up), "\n")

head(lusc_up[, c("gene_symbol", "log2FoldChange", "padj")], 20)
head(luad_up[, c("gene_symbol", "log2FoldChange", "padj")], 20)


# ============================================================
# STAGE 9: Volcano Plot — saved to PNG
# ============================================================
top_lusc_labels <- results_df %>%
  filter(log2FoldChange > 1, padj < 0.05) %>%
  arrange(padj) %>%
  head(15) %>%
  pull(gene_symbol)

top_luad_labels <- results_df %>%
  filter(log2FoldChange < -1, padj < 0.05) %>%
  arrange(padj) %>%
  head(15) %>%
  pull(gene_symbol)

genes_to_label <- c(top_lusc_labels, top_luad_labels)

png("volcano_LUSC_vs_LUAD.png", width = 4000, height = 3200, res = 300)

EnhancedVolcano(results_df,
                lab             = results_df$gene_symbol,
                selectLab       = genes_to_label,
                x               = "log2FoldChange",
                y               = "padj",
                title           = "LUSC vs LUAD Differential Expression",
                subtitle        = "Positive fold change = higher in LUSC",
                pCutoff         = 0.05,
                FCcutoff        = 1.0,
                pointSize       = 1.5,
                labSize         = 4.0,
                col             = c("grey70", "steelblue", "orange", "red2"),
                legendLabels    = c("NS", "FC only", "p-value only", "Significant"),
                drawConnectors  = TRUE,
                widthConnectors = 0.4,
                max.overlaps    = Inf,
                xlim            = c(-12, 12),
                ylim            = c(0, 350))

dev.off()
cat("Volcano plot saved.\n")


# ============================================================
# STAGE 10: Heatmap — saved to PNG
# ============================================================
normalized_counts <- counts(dds, normalized = TRUE)

# Get top 50 genes and keep only those present in count matrix
top50_genes <- results_df %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene_id)

top50_genes <- top50_genes[top50_genes %in% rownames(normalized_counts)]
cat("Genes in heatmap:", length(top50_genes), "\n")

top50_matrix <- normalized_counts[top50_genes, ]
top50_log    <- log2(top50_matrix + 1)
top50_scaled <- t(scale(t(top50_log)))
top50_scaled[is.nan(top50_scaled)] <- 0

# Use gene symbols as row labels
row_labels <- results_df$gene_symbol[match(top50_genes, results_df$gene_id)]
row_labels[is.na(row_labels)] <- top50_genes[is.na(row_labels)]

annotation_col <- data.frame(
  Subtype   = metadata$subtype,
  row.names = rownames(metadata)
)

annotation_colors <- list(
  Subtype = c(LUAD = "steelblue", LUSC = "firebrick")
)

png("heatmap_top50_DEGs.png", width = 4000, height = 4500, res = 300)

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

dev.off()
cat("Heatmap saved.\n")


# ============================================================
# STAGE 11: PCA Plot — saved to PNG
# ============================================================
vsd      <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "subtype", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = subtype)) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_color_manual(values = c(LUAD = "steelblue", LUSC = "firebrick")) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  ggtitle("PCA of TCGA LUAD vs LUSC Samples") +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank())

ggsave("PCA_LUAD_vs_LUSC.png", 
       plot   = pca_plot, 
       width  = 10, 
       height = 8, 
       dpi    = 300)

cat("PCA plot saved.\n")


# ============================================================
# STAGE 12: GO Pathway Enrichment — saved to PNG
# ============================================================
lusc_entrez <- mapIds(org.Hs.eg.db,
                      keys      = gsub("\\..*", "", lusc_up$gene_id),
                      column    = "ENTREZID",
                      keytype   = "ENSEMBL",
                      multiVals = "first")

luad_entrez <- mapIds(org.Hs.eg.db,
                      keys      = gsub("\\..*", "", luad_up$gene_id),
                      column    = "ENTREZID",
                      keytype   = "ENSEMBL",
                      multiVals = "first")

lusc_entrez <- na.omit(lusc_entrez)
luad_entrez <- na.omit(luad_entrez)

ego_lusc <- enrichGO(gene          = lusc_entrez,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     readable      = TRUE)

ego_luad <- enrichGO(gene          = luad_entrez,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     readable      = TRUE)

png("GO_enrichment_LUSC.png", width = 4000, height = 3500, res = 300)
print(dotplot(ego_lusc, showCategory = 20, 
              title = "GO BP: Pathways enriched in LUSC"))
dev.off()

png("GO_enrichment_LUAD.png", width = 4000, height = 3500, res = 300)
print(dotplot(ego_luad, showCategory = 20, 
              title = "GO BP: Pathways enriched in LUAD"))
dev.off()

cat("GO plots saved.\n")


# ============================================================
# STAGE 13: Save All Result Tables
# ============================================================
write.csv(sig_genes,               "NSCLC_significant_DEGs.csv",  row.names = FALSE)
write.csv(lusc_up,                 "genes_higher_in_LUSC.csv",    row.names = FALSE)
write.csv(luad_up,                 "genes_higher_in_LUAD.csv",    row.names = FALSE)
write.csv(as.data.frame(ego_lusc), "GO_enrichment_LUSC.csv",      row.names = FALSE)
write.csv(as.data.frame(ego_luad), "GO_enrichment_LUAD.csv",      row.names = FALSE)
saveRDS(dds,                       "dds_NSCLC_LUAD_LUSC.rds")

cat("All files saved to:", getwd(), "\n")
