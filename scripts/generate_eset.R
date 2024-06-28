express_matrix <- read.csv("01_expressMatrix.geneLevel_RSEM_TPM_filtered.csv", check.names=FALSE, row.names = "gene_id");
express_matrix$gene_name <- NULL
express_matrix$gene_type <- NULL
meta_data <- read.csv("manifest_source-data_RNA-Seq_MB.csv", row.names = "sample_name")
meta_data$id <- NULL
meta_data$name <- NULL
meta_data$size <- NULL
meta_data$project <- NULL
meta_data$visible <- NULL
meta_data$Kids.First.Participant.ID <- NULL
meta_data$disease_type <- NULL
meta_data$sample_id <- NULL
meta_data$sample_type <- NULL
meta_data$platform <- NULL
meta_data$Kids.First.Biospecimen.ID <- NULL
meta_data$aliquot_id <- NULL
meta_data$library_id <- NULL
meta_data$reference_genome <- NULL
meta_data$case_id <- NULL
meta_data$age_at_diagnosis <- NULL
meta_data$paired_end <- NULL
meta_data$X <- NULL
meta_data$experimental_strategy <- NULL
mat1 <- data.matrix(express_matrix)
mat2 <- data.matrix(meta_data)

generate.eset(
  exp_mat = mat1,
  phenotype_info = NULL,
  feature_info = NULL,
  annotation_info = ""
)
