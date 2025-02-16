library(dplyr)

setwd("/Users/abbiesiru/Desktop/research/medulloblastoma")

# read expression matrix file into R as a dataframe
tmp <- read.delim("dataset/01_expressMatrix.geneLevel_RSEM_TPM.txt", check.names=FALSE);
express_matrix <- dplyr::filter(tmp, gene_type=="protein_coding")

meta_data <- read.csv("dataset/manifest_source-data_RNA-Seq_MB.csv");

'%!in%' <- function(x,y)!('%in%'(x,y))

# initiate a empty vector that store the index of found sample in the meta data
meta_data_idx <- c();
rm_col_idx <- c();

# for (i in colnames(express_matrix)){
for (i in 4:ncol(express_matrix)) {
  if (colnames(express_matrix)[i] %!in% meta_data$sample_name) {
    rm_col_idx <- c(rm_col_idx, i);
  } 
  else {
    # return the index of the found sample in the meta data
    idx <- match(colnames(express_matrix)[i], meta_data$sample_name)
    # append the index to the end of the vector created above
    meta_data_idx <- c(meta_data_idx, idx )
    
  }
}

express_matrix[,rm_col_idx] <- NULL;
meta_data_subset <- meta_data[meta_data_idx,];

write.csv(express_matrix, "dataset/01_expressMatrix.geneLevel_RSEM_TPM_filtered.csv",row.names=FALSE)
write.csv(meta_data_subset, "dataset/manifest_source-data_RNA-Seq_MB_subset.csv",row.names=FALSE)

                     