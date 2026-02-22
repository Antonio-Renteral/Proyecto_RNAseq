# -------------------------
# Load libraries
# -------------------------
library(recount3)
library(DESeq2)

# -------------------------
# Load processed object
# -------------------------

rse_gene_SRP131263 <-
  readRDS(file = "processed-data/rse_gene_SRP131263")

# -------------------------
# Quality control (QC)
# -------------------------

rse_gene_SRP131263$assigned_gene_prop <-
  rse_gene_SRP131263$recount_qc.gene_fc_count_all.assigned /
  rse_gene_SRP131263$recount_qc.gene_fc_count_all.total

summary(rse_gene_SRP131263$assigned_gene_prop)

# The proportion of reads assigned to genes ranges between ~0.57 and ~0.61.
# No extreme outliers are observed, indicating acceptable gene-level assignment quality.

# -------------------------
# Metadata exploration
# -------------------------

colData(rse_gene_SRP131263)[, grepl(
  "^sra_attribute",
  colnames(colData(rse_gene_SRP131263))
)]

table(
  rse_gene_SRP131263$sra_attribute.cell_type,
  rse_gene_SRP131263$sra_attribute.genotype
)

# Genotypes are balanced within each cell type (3 vs 3 replicates).
# To avoid confounding effects from different Cre systems,
# we restrict the analysis to a single cell type (Act_CD4).

# -------------------------
# Subset to Act_CD4 cells
# -------------------------

rse_act <- rse_gene_SRP131263[,
  rse_gene_SRP131263$sra_attribute.cell_type == "Act_CD4"
]

table(rse_act$sra_attribute.genotype)

# Convert genotype to factor
colData(rse_act)$sra_attribute.genotype <-
  as.factor(colData(rse_act)$sra_attribute.genotype)

# Relevel factor so that wild type is the reference
colData(rse_act)$sra_attribute.genotype <-
  relevel(colData(rse_act)$sra_attribute.genotype, ref = "Bap1wt/wt_CD4Cre")

# -------------------------
# Ensure 'counts' is first assay (required by DESeq2)
# -------------------------

assays(rse_act) <- assays(rse_act)[c("counts", "raw_counts")]
assayNames(rse_act) # Confirm order

# -------------------------
# Filter low-count genes
# -------------------------

keep <- rowSums(assay(rse_act, "counts") >= 10) >= 3
rse_act <- rse_act[keep, ]

# -------------------------
# Construct DESeq2 object
# -------------------------

dds <- DESeqDataSet(rse_act, design = ~sra_attribute.genotype)

# Quick checks
dds
design(dds)
levels(dds$sra_attribute.genotype)

# -------------------------
# Run differential expression analysis
# -------------------------

dds <- DESeq(dds)

res <- results(dds)

summary(res)

# Order results by adjusted p-value
res <- res[order(res$padj), ]

# Inspect top genes
head(res)

# -------------------------
# Save model outputs
# -------------------------

saveRDS(dds, "processed-data/dds_act_cd4.rds")
saveRDS(res, "processed-data/res_act_cd4.rds")
