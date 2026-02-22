# -------------------------
# Load libraries
# -------------------------
library(DESeq2)
library(ggplot2)

# -------------------------
# Load model outputs
# -------------------------

dds <- readRDS("processed-data/dds_act_cd4.rds")
res <- readRDS("processed-data/res_act_cd4.rds")

# -------------------------
# PCA plot
# -------------------------

vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "sra_attribute.genotype", returnData = TRUE)

percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = sra_attribute.genotype)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

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
  theme_bw()

ggsave("plots/Volcano_plot.png", p_volcano, width = 6, height = 5)
