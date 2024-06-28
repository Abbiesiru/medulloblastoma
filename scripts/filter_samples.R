setwd("/Users/abbiesiru/Desktop/research/medulloblastoma/dataset")
library(readxl);

# read expression matrix file into R as a dataframe
express_matrix <- read.csv("01_expressMatrix.geneLevel_RSEM_TPM.csv", check.names=FALSE);
meta_data <- read_xlsx("manifest_source-data_RNA-Seq_MB.xlsx", col_names = TRUE);

sample_name2 = as.list(meta_data$sample_name)

'%!in%' <- function(x,y)!('%in%'(x,y))
meta_data_idx <- list()
# initiate a empty list that store the index of found sample in the meta data
# idx <- [];
for (i in colnames(express_matrix)){
  if(i %!in% sample_name2 && i!= "gene_id" && i!= "gene_name" && i!= "gene_type") {
    express_matrix[[i]] <- NULL
  } else if((i == "gene_id" | i== "gene_name" | i== "gene_type"))
    express_matrix[[i]]
  else {
    idx <- match(i, sample_name2)
    meta_data_idx <- c(meta_data_idx, idx )
    # return the index of the found sample in the meta data
    # append the index to the end of the list created above
  }
}

meta_data_subset <- data.frame(matrix(ncol = 24, nrow = 0))
colnames(meta_data_subset) <- c("id", "name", "sample_name", "size", "project", "visible", "Kids First Participant ID", "ethnicity", "gender", "race", "disease_type", "sample_id", "Tumor Descriptor", "sample_type", "paired_end", "platform", "Kids First Biospecimen ID", "primary_site", "age_at_diagnosis", "aliquot_id", "library_id", "reference_genome", "case_id", "experimental_strategy")
tmp <- data.frame(matrix(ncol = 24, nrow = 0))
colnames(tmp) <- c("id", "name", "sample_name", "size", "project", "visible", "Kids First Participant ID", "ethnicity", "gender", "race", "disease_type", "sample_id", "Tumor Descriptor", "sample_type", "paired_end", "platform", "Kids First Biospecimen ID", "primary_site", "age_at_diagnosis", "aliquot_id", "library_id", "reference_genome", "case_id", "experimental_strategy")

for (i in meta_data_idx){
 tmp <-  filter(meta_data, sample_name %in% sample_name2[i])
  meta_data_subset <- rbind(meta_data_subset, tmp)
}



write.csv(express_matrix, "01_expressMatrix.geneLevel_RSEM_TPM_filtered.csv")
write.csv(meta_data_subset, "manifest_source-data_RNA-Seq_MB_subset.csv")

# Output the filtered meta data, where the row ordering is the same as the column ordering of the expression data


                     