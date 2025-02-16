library(NetBID2);
setwd("/Users/abbiesiru/Desktop/research/medulloblastoma/")

screen_map <- read.csv('CRISPR/downloads/CRISPRScreenMap.csv', check.names=FALSE, header=TRUE);
tmp2 <- read.csv('CRISPR/downloads/CRISPRGeneEffect.csv', check.names=FALSE, header=TRUE);
gene_effect <- tmp2
tmp1 <- read.csv('CRISPR/downloads/CRISPRGeneDependency.csv', check.names=FALSE, header=TRUE);
gene_dependency <- tmp1
model <- read.xlsx('CRISPR/downloads/Model.xlsx')
G3_gene_list <- read.xlsx('surfaceome_proteins/G3_protein_list.xlsx')

colnames(gene_dependency) <- gsub(" \\(\\d+\\)$", "", colnames(tmp1)) # Remove "(number)" from row names
rownames(gene_dependency) <- gene_dependency[[1]]
gene_dependency <- gene_dependency[ , -1]
exprs_matrix <- t(gene_dependency)

colnames(gene_effect) <- gsub(" \\(\\d+\\)$", "", colnames(tmp2)) # Remove "(number)" from row names
rownames(gene_effect) <- gene_effect[[1]]
gene_effect <- gene_effect[ , -1]
fdata <- t(gene_effect)

rownames(model) <- model[[1]]
model <- model[ , -1]

common <- intersect(colnames(exprs_matrix), rownames(model))

model_filtered <- model[common, , drop = FALSE] # Filter model to keep only the common cell lines in the rows

fdata_df <- as.data.frame(fdata)

crispr_eset <- generate.eset(
  exp_mat = exprs_matrix,
  phenotype_info = model_filtered,
  feature_info = fdata_df,
)


