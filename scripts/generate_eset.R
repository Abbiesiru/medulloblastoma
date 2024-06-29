library(NetBID2);


setwd("/Users/abbiesiru/Desktop/research/medulloblastoma/")

tmp <- read.csv("dataset/01_expressMatrix.geneLevel_RSEM_TPM_filtered.csv", check.names=FALSE, row.names="gene_id");
exp_mat <- as.matrix(tmp[,3:ncol(tmp)]);
feature_info <- tmp[,1:2];

meta_data <- read.csv("dataset/manifest_source-data_RNA-Seq_MB_subset.csv", row.names = "sample_name")

net_eset <- generate.eset(
  exp_mat = exp_mat,
  phenotype_info = meta_data,
  feature_info = feature_info,
)

project_main_dir <- './' # user defined main directory for the project, one main directory could have multiple project folders, distinguished by project name.
project_name <- 'project'

network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

network.par$net.eset <- net_eset

draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='beforeQC_',emb_plot_type='2D.interactive')
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
