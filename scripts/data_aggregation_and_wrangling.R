library(NetBID2);
setwd("/Users/abbiesiru/Desktop/research/medulloblastoma/")

dataset1 <- readRDS("dataset/dataset1.rds")
dataset2 <- readRDS("dataset/northcott_N169_geneLevel_eset.rds")

#check if gene names are unique
tmp1 <- fData(dataset1)$gene_name
tmp2 <- unique(tmp1)
length(tmp1)
length(tmp2)

#data aggregation
mat1 <- exprs(dataset1)

df <- as.data.frame(mat1)
df$gene_names <- tmp1
a <- seq(1,117)
df[, a] <- apply(df[, a], 2, function(x) as.numeric(as.character(x)))

agg_df <- aggregate(df, list(df$gene_names), FUN = mean)
agg_df <- agg_df[, -ncol(agg_df)]

agg_mat <- as.matrix(agg_df[, -1])  # Exclude the first column (rownames)
rownames(agg_mat) <- agg_df[, 1]    # Set row names from the first column of agg_df

#save aggregated CTBN cohort
CTBN_agg <- generate.eset(
  exp_mat = agg_mat,
  phenotype_info = pData(dataset1),
)

#saveRDS(CTBN_agg, 'dataset/CTBN_agg.rds')

#Combine the 2 matrices:
#Identify common row names
mat2 <- exprs(dataset2)
common_rows <- intersect(rownames(agg_mat), rownames(mat2))

#Subset matrices based on common row names & z-transform
aggmat_common <- t(scale(t(agg_mat[common_rows, ])));
mat2_common <- t(scale(t(mat2[common_rows, ])));

exprsdata <- cbind(aggmat_common, mat2_common)


pdata1 <- pData(dataset1)
pdata1$COHORT <- "CTBN"
pdata1$gender <- ifelse(pdata1$gender == "Male", "M", "F")
names(pdata1)[names(pdata1) == "gender"] <- "GENDER"
pdata1_subset <- pdata1[, c("GENDER", "COHORT")]
pdata1_subset <- pdata1_subset[, c("COHORT", setdiff(names(pdata1_subset), "COHORT"))]
pdata1_subset$Subgroup <- NA

pdata2 <- pData(dataset2)
pdata2_subset <- pdata2[, c("COHORT", "GENDER", "Subgroup")]

combined_pdata <- rbind(pdata1_subset, pdata2_subset)

#generate eset
net_eset2 <- generate.eset(
  exp_mat = exprsdata,
  phenotype_info = combined_pdata,
)

project_main_dir <- './' # user defined main directory for the project, one main directory could have multiple project folders, distinguished by project name.
project_name <- 'project'

if (exists("network.par")) {
  rm(network.par)
}
network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)
network.par$net.eset <- net_eset2

#draw qc no.1
#draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='initial_',emb_plot_type='2D.interactive')

#retrieve expression matrix
exprsdata <- exprs(network.par$net.eset)

#normalize quantiles
exprsdata_norm <- normalizeQuantiles(exprsdata)

#generate eset
net_eset2 <- generate.eset(
  exp_mat = exprsdata_norm,
  phenotype_info = combined_pdata,
)
network.par$net.eset <- net_eset2

#draw qc no.2
#draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='norm_',emb_plot_type='2D.interactive')


#Remove batch effect
exprs_rbe <- limma::removeBatchEffect(net_eset2, net_eset2$COHORT)

# Gene filtering:
# filter out genes with very low expression values (bottom 25%) in most samples (more than 60%).
choose1 <- apply(exprs_rbe<=quantile(exprs_rbe, probs = 0.25), 1, sum)<=(ncol(exprs_rbe)*0.60)
print(table(choose1))
exprs_rbe <- exprs_rbe[choose1,]

# Select the most variable genes across samples
choose1 <- IQR.filter(exp_mat=exprs_rbe,use_genes=rownames(exprs_rbe),thre = 0.5)
print(table(choose1))
exprs_rbe_filt <- exprs_rbe[choose1,]

exprs_rbe_filt<- normalizeQuantiles(exprs_rbe_filt)

net_eset2_rbe <- generate.eset(
  exp_mat = exprs_rbe_filt,
  phenotype_info = pData(net_eset2),
)

network.par$net.eset <- net_eset2_rbe

#draw qc no.3
#draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='norm_rbe_',emb_plot_type='2D.interactive')


### K means clustering ###

#pdf(file="project/DATA/batchcorrection.pdf",width=16,height=8)
cluster_label <- draw.emb.kmeans(
  mat = exprs(net_eset2),
  embedding_method = "pca",
  all_k = 4,
  obs_label = get_obs_label(pData(net_eset2_rbe), use_col="COHORT"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "basic",
  choose_k_strategy = "Jaccard",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
#dev.off()

table(cluster_label, pData(net_eset2_rbe)$Subgroup)
cluster_label_backup <- cluster_label

for(i in 1:nrow(pData(net_eset2_rbe))){
  if(cluster_label[i]==1 & is.na(pData(net_eset2_rbe)$Subgroup[i])){
    pData(net_eset2_rbe)$Subgroup[i] <- "G3"
  } else if(cluster_label[i]==2 & is.na(pData(net_eset2_rbe)$Subgroup[i])){
    pData(net_eset2_rbe)$Subgroup[i] <- "G4"
  } else if(cluster_label[i]==3 & is.na(pData(net_eset2_rbe)$Subgroup[i])){
    pData(net_eset2_rbe)$Subgroup[i] <- "SHH"
  } else if(cluster_label[i]==4 & is.na(pData(net_eset2_rbe)$Subgroup[i])){
    pData(net_eset2_rbe)$Subgroup[i] <- "WNT"
  }
}
View(pData(net_eset2_rbe))

network.par$net.eset <- net_eset2_rbe

saveRDS(net_eset2_rbe, 'dataset/dataset_combined.rds')
