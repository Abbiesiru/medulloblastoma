library(NetBID2);

tmp <- read.csv("01_expressMatrix.geneLevel_RSEM_TPM_filtered.csv", check.names=FALSE, row.names = "gene_id");
exp_mat <- as.matrix(tmp[,4:110]);
feature_info <- tmp[,2:3];

meta_data <- read.csv("manifest_source-data_RNA-Seq_MB_subset.csv", row.names = "sample_name")
meta_data$X <- NULL

eset <- generate.eset(
  exp_mat = exp_mat,
  phenotype_info = meta_data,
  feature_info = feature_info,
)
