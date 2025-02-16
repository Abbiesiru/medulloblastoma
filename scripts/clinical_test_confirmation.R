library(NetBID2);
library(Biobase);
library(GEOquery);
library(readxl);
library(survival);
library(survminer);
library(dplyr);
library(forcats);
setwd("/Users/abbiesiru/Desktop/research/medulloblastoma/")


### load gene expression and network files and calculate driver activity ### 

# 1. load clinical data
clinical <- read_excel("dataset/GSE85217_clinical.xlsx", col_names = TRUE) 

# 2. load eset
load("dataset/GSE85217.raw.eset")
cal.eset <- generate.eset(exp_mat=exprs(GSE85217.raw.eset),phenotype_info=pData(GSE85217.raw.eset),
                          feature_info=fData(GSE85217.raw.eset))

# 3. ensure alignment of data
cal.eset <- cal.eset[, pData(cal.eset)$title %in% clinical$Study_ID]
clinical <- clinical[match(pData(cal.eset)$title, clinical$Study_ID), ]

# 4. read in network file
network <- readRDS("network/GSE85217_allMB_network.rds")

# 5. convert ENSEMBL gene IDs to gene symbols
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # remove version numbers
ensembl_ids <- rownames(exprs(cal.eset))
ensembl_ids <- gsub("_at", "", ensembl_ids)  # remove suffix
conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values = ensembl_ids,
                    mart = mart)
cal_mat <- exprs(cal.eset)
rownames(cal_mat) <- gsub("_at$", "", rownames(cal_mat))
cal_mat <- cal_mat[rownames(cal_mat) %in% conversion$ensembl_gene_id, ]
rownames(cal_mat) <- conversion$hgnc_symbol[match(rownames(cal_mat), conversion$ensembl_gene_id)]
unmatched_ids <- setdiff(ensembl_ids, conversion$ensembl_gene_id) # get unmatched IDs

# 6. get activity matrix
ac_mat <- cal.Activity(target_list=network$target_list,cal_mat=cal_mat,es.method='weightedmean')

# 7. create eset using activity matrix
merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(cal.eset)[colnames(ac_mat),],
                               feature_info=NULL,annotation_info='activity in net-dataset')


### get differential expression (DE) / differential activity (DA) for drivers ###

# 1. create empty list to store comparison result
DE_G3 <- list()
DA_G3 <- list()

# 2. get sample names from each compared group
comp_name <- 'G3.Vs.others'
phe_info <- pData(cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`characteristics_ch1.1`=='subgroup: Group3')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`characteristics_ch1.1`!='subgroup: Group3')] # Combine other groups as the Control group
DE_gene_bid <- getDE.BID.2G(eset=cal.eset,G1=G1,G0=G0,G1_name='G3',G0_name='others')
DA_driver_bid   <- getDE.BID.2G(eset=merge.ac.eset,G1=G1,G0=G0,G1_name='G3',G0_name='others')

# 3. save comparison result and reload data into R workspace
DE_G3[[comp_name]] <- DE_gene_bid
DA_G3[[comp_name]] <- DA_driver_bid
db.preload(use_level='gene',use_spe='human',update=FALSE)

# 4. get all comparison names
all_comp <- names(DE_G3) 

# 5. create  final master table
ms_tab_G3 <- generate.masterTable(use_comp=all_comp,DE=DE_G3,DA=DA_G3,
                                   target_list=network$target_list,
                                   tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                                   main_id_type='external_gene_name')


### survival analysis for DE ### 

# 1. create survival object
clinical <- clinical[!(clinical$'OS (years)' == "NA"), ] # remove NA vals
clinical <- clinical[!(clinical$'Met status (1 Met, 0 M0)' == "NA"), ] # remove NA vals
surv_obj <- Surv(time = as.numeric(as.character(clinical$'OS (years)')), event = as.numeric(as.character(clinical$'Met status (1 Met, 0 M0)')))

# 2. fit survival curves
high_low_group <- ifelse(exprs(cal.eset)['ENSG00000081479_at', ] > median(exprs(cal.eset)['ENSG00000081479_at', ]), "High", "Low")
clinical$group <- high_low_group
fit <- survfit(surv_obj ~ group, data = clinical)

# 3. plot survival curves
ggsurvplot(fit, data = clinical, pval = TRUE, risk.table = TRUE)

# 4. create g3 eset
clinical_G3 <- clinical[(clinical$'Subgroup' == "Group3"), ] # g3 only
titles <- pData(cal.eset)$title
study_ids <- clinical$Study_ID
filtered_pdata <- pData(cal.eset)[titles %in% study_ids, ] # filter pData(cal.eset)
exprs_colnames <- colnames(exprs(cal.eset)) # extract column names from exprs(cal.eset)
pdata_rownames <- rownames(filtered_pdata) # extract row names from pData(net.eset)
matching_columns <- exprs_colnames %in% pdata_rownames # identify matching columns
filtered_exprs <- exprs(cal.eset)[, matching_columns] # subset exprs(cal.eset) to keep only matching columns
G3.eset <- generate.eset(exp_mat=filtered_exprs, phenotype_info=filtered_pdata)

# 5. filter clinical data and eset for common samples
G3.eset <- G3.eset[, pData(G3.eset)$title %in% clinical_G3$Study_ID]
clinical_G3 <- clinical_G3[match(pData(G3.eset)$title, clinical_G3$Study_ID), ]

# 6. group based on median expression
high_low_group_G3 <- ifelse(exprs(G3.eset)['ENSG00000081479_at', ] > median(exprs(G3.eset)['ENSG00000081479_at', ]), "High", "Low")
clinical_G3$group <- high_low_group_G3

# 7. create survival object, then fit & plot survival curves
surv_obj_G3 <- Surv(time = as.numeric(as.character(clinical_G3$'OS (years)')), event = as.numeric(as.character(clinical_G3$'Met status (1 Met, 0 M0)')))
fit_G3 <- survfit(surv_obj_G3 ~ group, data = clinical_G3)
ggsurvplot(fit_G3, data = clinical_G3, pval = TRUE, risk.table = TRUE)


### survival analysis for DA ### 

# 1. extract DA gene & activity data
da_gene <- ms_tab_G3$geneSymbol[377]
da_activity <- exprs(merge.ac.eset)[5400, ]

# 2. filter clinical data and DA eset for common samples
clinical <- clinical[clinical$Study_ID %in% pData(merge.ac.eset)$title, ]
merge.ac.eset <- merge.ac.eset[, pData(merge.ac.eset)$title %in% clinical$Study_ID]
clinical <- clinical[match(pData(merge.ac.eset)$title, clinical$Study_ID), ]

# 3. group based on median activity
high_low_group_DA <- ifelse(da_activity > median(da_activity), "High", "Low")
clinical$group <- high_low_group_DA

# 4. create survival object, then fit & plot survival curves
surv_obj_DA <- Surv(time = as.numeric(as.character(clinical$`OS (years)`)), 
                    event = as.numeric(as.character(clinical$`Met status (1 Met, 0 M0)`)))
fit_DA <- survfit(surv_obj_DA ~ group, data = clinical)
ggsurvplot(fit_DA, data = clinical, pval = TRUE, risk.table = TRUE, 
           title = paste("Survival Analysis for DA Gene:", da_gene))
