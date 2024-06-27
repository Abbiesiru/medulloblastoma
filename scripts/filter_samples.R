# load the expression data matrix first into a variable called "X01_expressMatrix_geneLevel_RSEM_TPM" or something simpler


sample_name2 = as.list(manifest_source_data_RNA_Seq_MB$sample_name)

'%!in%' <- function(x,y)!('%in%'(x,y))
for (i in colnames(X01_expressMatrix_geneLevel_RSEM_TPM)){
  if(i %!in% sample_name2 && i!= "gene_id" && i!= "gene_name" && i!= "gene_type") {
    X01_expressMatrix_geneLevel_RSEM_TPM[[i]] <- NULL
  }
}
colnames(X01_expressMatrix_geneLevel_RSEM_TPM)

# Save the filered data into a new file

                     