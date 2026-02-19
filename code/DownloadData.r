# Load library
library(recount3)

# Check mouse available projects and save the result in a variable
mouse_project <- available_projects(organism = "mouse")

# I will use SRP131263 for the analysis because it has 30 samples, providing sufficient replication
# for reliable differential expression while remaining computationally manageable.

# Download the data and save it in a variable
proj_info <- subset(
  mouse_project,
  project == "SRP131263" & project_type == "data_sources"
)

# Create a RangedSummarizedExperiment object with the data
rse_gene_SRP131263 <- create_rse(proj_info)
rse_gene_SRP131263


# Compute read counts and save the result in the assay slot of the RangedSummarizedExperiment object
assay(rse_gene_SRP131263, "counts") <- compute_read_counts(rse_gene_SRP131263)

# Expand the SRA attributes and save the result in the colData slot of the RangedSummarizedExperiment object
rse_gene_SRP131263 <- expand_sra_attributes(rse_gene_SRP131263)


# Save the RangedSummarizedExperiment object in an RDS file
saveRDS(create_rse(proj_info), file = "raw-data/raw_rse_gene_SRP131263") # Raw RSE
saveRDS(rse_gene_SRP131263, file = "processed-data/rse_gene_SRP131263") # Processed RSE
