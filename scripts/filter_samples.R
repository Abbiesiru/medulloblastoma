# load the expression data matrix first into a variable called "X01_expressMatrix_geneLevel_RSEM_TPM" or something simpler
readr::read_delim

express_matrix <- read.csv("dataset/01_expressMatrix.geneLevel_RSEM_TPM.csv", check.names=FALSE);
meta_data <- read_xlsx("dataset/manifest_source-data_RNA-Seq_MB.xlsx", col_names = TRUE);

sample_name2 = as.list(meta_data$sample_name)

'%!in%' <- function(x,y)!('%in%'(x,y))
for (i in colnames(express_matrix)){
  if(i %!in% sample_name2 && i!= "gene_id" && i!= "gene_name" && i!= "gene_type") {
    express_matrix[[i]] <- NULL
  }
}
colnames(express_matrix)

# Save the filered data into a new file that can be loaded back into R as an ExpressionSet class object

                     