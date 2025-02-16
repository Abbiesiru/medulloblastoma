library(NetBID2);

setwd("/Users/abbiesiru/Desktop/research/medulloblastoma/")

tmp <- read.csv("dataset/01_expressMatrix.geneLevel_RSEM_TPM_filtered.csv", check.names=FALSE, row.names="gene_id");
exp_mat <- as.matrix(tmp[,3:ncol(tmp)]);
feature_info <- tmp[,1:2];

meta_data <- read.csv("dataset/manifest_source-data_RNA-Seq_MB_subset.csv", row.names = "sample_name")

# These are the meta data columns to use in the QC
idx <- c("ethnicity","gender","race","Tumor.Descriptor","age_at_diagnosis"); # "primary_site"
meta_data <- meta_data[,idx];

#log2-transformation
net_eset <- generate.eset(
  exp_mat = log2(exp_mat+1),
  phenotype_info = meta_data,
  feature_info = feature_info,
)

# exprs(net_eset)<-log2(exprs(net_eset)+1) 

project_main_dir <- './' # user defined main directory for the project, one main directory could have multiple project folders, distinguished by project name.
project_name <- 'project'

if (exists("network.par")) {
  rm(network.par)
}
network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)
network.par$net.eset <- net_eset

draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='beforeQC_',emb_plot_type='2D.interactive')

###  The following QC steps are optional
## Firstly, the handling of missing data.
# Get the expression matrix from ExpressionSet object
mat <- exprs(network.par$net.eset)
# Count and show number of NAs across samples and genes
sample_na_count <- apply(mat,1,function(x){length(which(is.na(x)==TRUE))})
print(table(sample_na_count))
gene_na_count <- apply(mat,2,function(x){length(which(is.na(x)==TRUE))})
print(table(gene_na_count))
# Perform imputation
if(sum(sample_na_count)+sum(gene_na_count)>0) mat <- impute.knn(mat)$data

## Secondly, the log2 transformation.
# med_val <- median(apply(mat,2,median)); print(med_val)
# if(med_val>16){mat <- log2(mat)}

## Thirdly, the quantile normalization across samples.
# Perform limma quantile normalization
mat_norm <- normalizeQuantiles(mat)

## Fourthly, filter out genes with very low expression values (bottom 5%) in most samples (more than 90%).
# Filter out low-expression genes
choose1 <- apply(mat_norm<=quantile(mat_norm, probs = 0.05), 1, sum)<=(ncol(mat_norm)*0.80)
print(table(choose1))
mat_norm_filt <- mat_norm[choose1,]

# Update eset with normalized expression matrixs
net_eset <- generate.eset(exp_mat=mat_norm_filt, phenotype_info=pData(network.par$net.eset)[colnames(mat_norm_filt),],
                          feature_info=fData(network.par$net.eset)[rownames(mat_norm_filt),])
# Update network.par with new eset
network.par$net.eset <- net_eset

# QC for the normalized eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='afterQC_',emb_plot_type='2D.interactive')

saveRDS(net_eset, 'dataset/dataset1.rds')


pdf(file="project/DATA/pca_basic_ARI.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "pca",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "basic",
  choose_k_strategy = "ARI",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
dev.off()

pdf(file="project/DATA/pca_basic_NMI.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "pca",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "basic",
  choose_k_strategy = "NMI",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
dev.off()

pdf(file="project/DATA/pca_basic_Jaccard.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "pca",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
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
dev.off()

pdf(file="project/DATA/pca_consensus_ARI.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "pca",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "consensus",
  choose_k_strategy = "ARI",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
dev.off()

pdf(file="project/DATA/pca_consensus_NMI.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "pca",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "consensus",
  choose_k_strategy = "NMI",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
dev.off()

pdf(file="project/DATA/pca_consensus_Jaccard.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "pca",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "consensus",
  choose_k_strategy = "Jaccard",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
dev.off()

pdf(file="project/DATA/mds_basic_ARI.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "mds",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "basic",
  choose_k_strategy = "ARI",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
dev.off()

pdf(file="project/DATA/mds_basic_NMI.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "mds",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "basic",
  choose_k_strategy = "NMI",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
dev.off()

pdf(file="project/DATA/mds_basic_Jaccard.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "mds",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
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
dev.off()

pdf(file="project/DATA/mds_consensus_ARI.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "mds",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "consensus",
  choose_k_strategy = "ARI",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
dev.off()

pdf(file="project/DATA/mds_consensus_NMI.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "mds",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "consensus",
  choose_k_strategy = "NMI",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
dev.off()

pdf(file="project/DATA/mds_consensus_Jaccard.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "mds",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "consensus",
  choose_k_strategy = "Jaccard",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
dev.off()

pdf(file="project/DATA/umap_basic_ARI.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "umap",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "basic",
  choose_k_strategy = "ARI",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
dev.off()

pdf(file="project/DATA/umap_basic_NMI.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "umap",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "basic",
  choose_k_strategy = "NMI",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
dev.off()

pdf(file="project/DATA/umap_basic_Jaccard.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "umap",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
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
dev.off()

pdf(file="project/DATA/umap_consensus_ARI.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "umap",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "consensus",
  choose_k_strategy = "ARI",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
dev.off()

pdf(file="project/DATA/umap_consensus_NMI.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "umap",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "consensus",
  choose_k_strategy = "NMI",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
dev.off()

pdf(file="project/DATA/umap_consensus_Jaccard.pdf",width=16,height=8)
draw.emb.kmeans(
  mat = exprs(net_eset),
  embedding_method = "umap",
  all_k = NULL,
  obs_label = get_obs_label(meta_data, use_col="ethnicity"),
  legend_pos = "topleft",
  legend_cex = 0.8,
  plot_type = "2D.ellipse",
  point_cex = 1,
  kmeans_strategy = "consensus",
  choose_k_strategy = "Jaccard",
  return_type = "optimal",
  main = "",
  verbose = TRUE,
  use_color = NULL,
  pre_define = NULL
)
dev.off()



# If use the same expression dataset as in the network construction, just reload it directly
load(sprintf('%s/DATA/network.par.Step.exp-QC.RData',network.dir)) # RData saved after QC in the network construction step
analysis.par$cal.eset <- network.par$net.eset

# Save Step 1 network.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='exp-QC')

NetBID.loadRData(analysis.par=analysis.par,step='exp-QC')

# Get network information
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)

# Creat QC report for the network
draw.network.QC(analysis.par$tf.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='TF_net_',html_info_limit=FALSE)
draw.network.QC(analysis.par$sig.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='SIG_net_',html_info_limit=TRUE)

