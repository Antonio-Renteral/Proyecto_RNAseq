# -------------------------
# Load libraries
# -------------------------
library(DESeq2)
library(ggplot2)
library(pheatmap)

# -------------------------
# Load model outputs
# -------------------------

dds <- readRDS("processed-data/dds_act_cd4.rds")
res <- readRDS("processed-data/res_act_cd4.rds")

# -------------------------
# Variance stabilizing transformation
# -------------------------

vsd <- vst(dds, blind = FALSE)

# -------------------------
# PCA plot
# -------------------------

pca_data <- plotPCA(vsd, intgroup = "sra_attribute.genotype", returnData = TRUE)

percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = sra_attribute.genotype)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(text = element_text(size = 14))

ggsave("plots/PCA_plot.png", p_pca, width = 6, height = 5)

# -------------------------
# MA plot
# -------------------------

png("plots/MA_plot.png", width = 800, height = 600)
plotMA(res, ylim = c(-5, 5))
dev.off()

# -------------------------
# Volcano plot
# -------------------------

res_df <- as.data.frame(res)

# Remove NA padj values
res_df <- res_df[!is.na(res_df$padj), ]

res_df$significant <- "Not Sig"
res_df$significant[
  res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) > 1
] <- "Significant"

p_volcano <- ggplot(
  res_df,
  aes(x = log2FoldChange, y = -log10(padj), color = significant)
) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 14))

ggsave("plots/Volcano_plot.png", p_volcano, width = 6, height = 5)

# -------------------------
# Heatmap of top 30 DE genes
# -------------------------

# Order by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Select top 30 genes
top_genes <- head(rownames(res_ordered), 30)

# Extract normalized expression values
mat <- assay(vsd)[top_genes, ]

# Scale by gene (row-wise z-score)
mat_scaled <- t(scale(t(mat)))

# Sample annotation
annotation_col <- data.frame(
  Genotype = colData(dds)$sra_attribute.genotype
)
rownames(annotation_col) <- colnames(mat_scaled)

png("plots/Heatmap_top30.png", width = 800, height = 900)

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  fontsize_col = 12
)

dev.off()
