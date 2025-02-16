library(NetBID2);
library(dplyr);
library(forcats);
setwd("/Users/abbiesiru/Desktop/research/medulloblastoma/")

############### Step 1: Load in gene expression dataset for analysis (exp-load, exp-cluster, exp-QC) ###############

CTBN_eset <- readRDS("dataset/CTBN_agg.rds")
cal.eset <- readRDS("dataset/dataset_combined.rds")
phenotype1 <- pData(cal.eset)
phenotype<- phenotype1[-seq(nrow(phenotype1),nrow(phenotype1)-168),]
pData(CTBN_eset)$Subgroup <- phenotype$Subgroup

############### Step 2: Read in network files and calcualte driver activity (act-get) ###############

network <- readRDS("network/GSE85217_allMB_network.rds")

# Get activity matrix
ac_mat <- cal.Activity(target_list=network$target_list,cal_mat=exprs(CTBN_eset),es.method='weightedmean')

# Create eset using activity matrix
merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(CTBN_eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

# QC plot for activity eset
#draw.eset.QC(merge.ac.eset,outdir="network/",intgroup=NULL,do.logtransform=FALSE,prefix='CTBN_AC_',emb_plot_type ='2D.interactive')


############### Step 3: Get differential expression (DE) / differential activity (DA) for drivers (act-DA) ###############


#WNT

# Create empty list to store comparison result
DE_WNT <- list()
DA_WNT <- list()
DE_all <- list()
DA_all <- list()

comp_name <- 'WNT.Vs.others'
# Get sample names from each compared group
phe_info <- pData(CTBN_eset)
G1  <- rownames(phe_info)[which(phe_info$`Subgroup`=='WNT')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`Subgroup`!='WNT')] # Combine other groups as the Control group
DE_gene_bid_WNT <- getDE.BID.2G(eset=CTBN_eset,G1=G1,G0=G0,G1_name='WNT',G0_name='others')
DA_driver_bid_WNT   <- getDE.BID.2G(eset=merge.ac.eset,G1=G1,G0=G0,G1_name='WNT',G0_name='others')
# Save comparison result to list element in analysis.par, with comparison name
DE_WNT[[comp_name]] <- DE_gene_bid_WNT
DE_all[[comp_name]] <- DE_gene_bid_WNT
DA_WNT[[comp_name]] <- DA_driver_bid_WNT
DA_all[[comp_name]] <- DA_driver_bid_WNT

# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE)
# Get all comparison names
all_comp <- names(DE_WNT) # Users can use index or name to get target ones

# Create the final master table
ms_tab_WNT <- generate.masterTable(use_comp=all_comp,DE=DE_WNT,DA=DA_WNT,
                               target_list=network$target_list,
                               tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                               main_id_type='external_gene_name')

# Visualize top genes
sig_gene_WNT <- draw.volcanoPlot(dat=ms_tab_WNT,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
                                 Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=2,Pv_thre=1e-5,
                                 main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=TRUE,
                                 pdf_file='volcano_plot_CTBN_WNT.pdf',label_cex = 0.5)

gs.preload(use_spe='Homo sapiens',update=FALSE)
# Gene Set Enrichment Analysis
gene_list_up <- rownames(sig_gene_WNT)[which(sig_gene_WNT[,2]>0)] # up
gene_list_down <- rownames(sig_gene_WNT)[which(sig_gene_WNT[,2]<0)] # down
res_up <- funcEnrich.Fisher(input_list=ms_tab_WNT[gene_list_up,'geneSymbol'],bg_list=unique(ms_tab_WNT[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
res_down <- funcEnrich.Fisher(input_list=ms_tab_WNT[gene_list_down,'geneSymbol'],bg_list=unique(ms_tab_WNT[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       

# Gene set enrichment analysis Barplot
draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top genes in WNT Subgroup',display_genes = TRUE,gs_cex=0.6,
                    pdf_file='funcEnrich_bar_WNT_CTBN.pdf')

filtered_ms_tab_WNT <- ms_tab_WNT %>%
  filter(P.Value.WNT.Vs.others_DE < 1e-5  & logFC.WNT.Vs.others_DE > 2)

# Visualize top drivers
sig_driver_WNT <- draw.volcanoPlot(dat=ms_tab_WNT,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DA',comp_name),
                                   Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.15,Pv_thre=1e-15,
                                   main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,
                                   pdf_file='volcano_plot_CTBN_WNT_DA.pdf',label_cex = 0.5)

gs.preload(use_spe='Homo sapiens',update=FALSE)
# Gene Set Enrichment Analysis
# driver_list_up <- rownames(sig_driver_WNT)[which(sig_driver_WNT[,2]>0)] # up
# driver_list_down <- rownames(sig_driver_WNT)[which(sig_driver_WNT[,2]<0)] # down
# res_up <- funcEnrich.Fisher(input_list=ms_tab_WNT[driver_list_up,'geneSymbol'],bg_list=unique(ms_tab_WNT[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
# res_down <- funcEnrich.Fisher(input_list=ms_tab_WNT[driver_list_down,'geneSymbol'],bg_list=unique(ms_tab_WNT[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       
# 
# # Gene set enrichment analysis Barplot
# draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers in WNT Subgroup',display_genes = TRUE,gs_cex=0.6,
#                     pdf_file='funcEnrich_bar_WNT_CTBN_DA.pdf')

filtered_ms_tab_WNT_DA <- ms_tab_WNT %>%
  filter(P.Value.WNT.Vs.others_DA < 1e-15  & logFC.WNT.Vs.others_DA > 0.15)



#SHH

# Create empty list to store comparison result
DE_SHH <- list()
DA_SHH <- list()

comp_name <- 'SHH.Vs.others'
# Get sample names from each compared group
phe_info <- pData(CTBN_eset)
G1  <- rownames(phe_info)[which(phe_info$`Subgroup`=='SHH')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`Subgroup`!='SHH')] # Combine other groups as the Control group
DE_gene_bid_SHH <- getDE.BID.2G(eset=CTBN_eset,G1=G1,G0=G0,G1_name='SHH',G0_name='others')
DA_driver_bid_SHH   <- getDE.BID.2G(eset=merge.ac.eset,G1=G1,G0=G0,G1_name='SHH',G0_name='others')
# Save comparison result to list element in analysis.par, with comparison name
DE_SHH[[comp_name]] <- DE_gene_bid_SHH
DE_all[[comp_name]] <- DE_gene_bid_SHH
DA_SHH[[comp_name]] <- DA_driver_bid_SHH
DA_all[[comp_name]] <- DA_driver_bid_SHH

# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE)
# Get all comparison names
all_comp <- names(DE_SHH) # Users can use index or name to get target ones

# Create the final master table
ms_tab_SHH <- generate.masterTable(use_comp=all_comp,DE=DE_SHH,DA=DA_SHH,
                               target_list=network$target_list,
                               tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                               main_id_type='external_gene_name')

# Visualize top genes
sig_gene_SHH <- draw.volcanoPlot(dat=ms_tab_SHH,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
                                 Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=2,Pv_thre=1e-5,
                                 main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=TRUE,
                                 pdf_file='volcano_plot_CTBN_SHH.pdf',label_cex = 0.5)

gs.preload(use_spe='Homo sapiens',update=FALSE)
# Gene Set Enrichment Analysis
# gene_list_up <- rownames(sig_gene_SHH)[which(sig_gene_SHH[,2]>0)] # up
# gene_list_down <- rownames(sig_gene_SHH)[which(sig_gene_SHH[,2]<0)] # down
# res_up <- funcEnrich.Fisher(input_list=ms_tab_SHH[gene_list_up,'geneSymbol'],bg_list=unique(ms_tab_SHH[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
# res_down <- funcEnrich.Fisher(input_list=ms_tab_SHH[gene_list_down,'geneSymbol'],bg_list=unique(ms_tab_SHH[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       
# 
# # Gene set enrichment analysis Barplot
# draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top genes in SHH Subgroup',display_genes = TRUE,gs_cex=0.6,
#                     pdf_file='funcEnrich_bar_SHH_CTBN.pdf')

filtered_ms_tab_SHH <- ms_tab_SHH %>%
  filter(P.Value.SHH.Vs.others_DE < 1e-5  & logFC.SHH.Vs.others_DE > 2)

# Visualize top drivers
sig_driver_SHH <- draw.volcanoPlot(dat=ms_tab_SHH,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DA',comp_name),
                                   Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.15,Pv_thre=1e-15,
                                   main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,
                                   pdf_file=NULL,label_cex = 0.5)

gs.preload(use_spe='Homo sapiens',update=FALSE)
# Gene Set Enrichment Analysis
driver_list_up <- rownames(sig_driver_SHH)[which(sig_driver_SHH[,2]>0)] # up
driver_list_down <- rownames(sig_driver_SHH)[which(sig_driver_SHH[,2]<0)] # down
res_up <- funcEnrich.Fisher(input_list=ms_tab_SHH[driver_list_up,'geneSymbol'],bg_list=unique(ms_tab_SHH[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
res_down <- funcEnrich.Fisher(input_list=ms_tab_SHH[driver_list_down,'geneSymbol'],bg_list=unique(ms_tab_SHH[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       

# Gene set enrichment analysis Barplot
draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers in SHH Subgroup',display_genes = TRUE,gs_cex=0.6,
                    pdf_file='funcEnrich_bar_SHH_CTBN_DA.pdf')

filtered_ms_tab_SHH_DA <- ms_tab_SHH %>%
  filter(P.Value.SHH.Vs.others_DA < 1e-15  & logFC.SHH.Vs.others_DA > 0.15)


#G3

# Create empty list to store comparison result
DE_G3 <- list()
DA_G3 <- list()

comp_name <- 'G3.Vs.others'
# Get sample names from each compared group
phe_info <- pData(CTBN_eset)
G1  <- rownames(phe_info)[which(phe_info$`Subgroup`=='G3')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`Subgroup`!='G3')] # Combine other groups as the Control group
DE_gene_bid_G3 <- getDE.BID.2G(eset=CTBN_eset,G1=G1,G0=G0,G1_name='G3',G0_name='others')
DA_driver_bid_G3   <- getDE.BID.2G(eset=merge.ac.eset,G1=G1,G0=G0,G1_name='G3',G0_name='others')
# Save comparison result to list element in analysis.par, with comparison name
DE_G3[[comp_name]] <- DE_gene_bid_G3
#DE_all[[comp_name]] <- DE_gene_bid_G3
DA_G3[[comp_name]] <- DA_driver_bid_G3
#DA_all[[comp_name]] <- DA_driver_bid_G3

# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE)
# Get all comparison names
all_comp <- names(DE_G3) # Users can use index or name to get target ones

# Create the final master table
ms_tab_G3 <- generate.masterTable(use_comp=all_comp,DE=DE_G3,DA=DA_G3,
                               target_list=network$target_list,
                               tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                               main_id_type='external_gene_name')

# Visualize top genes
sig_gene_G3 <- draw.volcanoPlot(dat=ms_tab_G3,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
                                Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=2,Pv_thre=1e-5,
                                main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=TRUE,
                                pdf_file='volcano_plot_CTBN_G3.pdf',label_cex = 0.5)

gs.preload(use_spe='Homo sapiens',update=FALSE)
# Gene Set Enrichment Analysis
# gene_list_up <- rownames(sig_gene_G3)[which(sig_gene_G3[,2]>0)] # up
# gene_list_down <- rownames(sig_gene_G3)[which(sig_gene_G3[,2]<0)] # down
# res_up <- funcEnrich.Fisher(input_list=ms_tab_G3[gene_list_up,'geneSymbol'],bg_list=unique(ms_tab_G3[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
# res_down <- funcEnrich.Fisher(input_list=ms_tab_G3[gene_list_down,'geneSymbol'],bg_list=unique(ms_tab_G3[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       
# 
# # Gene set enrichment analysis Barplot
# draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top genes in G3 Subgroup',display_genes = TRUE,gs_cex=0.6,
#                     pdf_file='funcEnrich_bar_G3_CTBN.pdf')

filtered_ms_tab_G3 <- ms_tab_G3 %>%
  filter(P.Value.G3.Vs.others_DE < 1e-5  & logFC.G3.Vs.others_DE > 2)

# Visualize top drivers
sig_driver_G3 <- draw.volcanoPlot(dat=ms_tab_G3,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DA',comp_name),
                                  Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.1,Pv_thre=1e-5,
                                  main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,
                                  pdf_file=NULL,label_cex = 0.5)

gs.preload(use_spe='Homo sapiens',update=FALSE)
# Gene Set Enrichment Analysis
driver_list_up <- rownames(sig_driver_G3)[which(sig_driver_G3[,2]>0)] # up
driver_list_down <- rownames(sig_driver_G3)[which(sig_driver_G3[,2]<0)] # down
res_up <- funcEnrich.Fisher(input_list=ms_tab_G3[driver_list_up,'geneSymbol'],bg_list=unique(ms_tab_G3[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
res_down <- funcEnrich.Fisher(input_list=ms_tab_G3[driver_list_down,'geneSymbol'],bg_list=unique(ms_tab_G3[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       

# Gene set enrichment analysis Barplot
draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers in G3 Subgroup',display_genes = TRUE,gs_cex=0.6,
                    pdf_file='funcEnrich_bar_G3_CTBN_DA.pdf')

filtered_ms_tab_G3_DA <- ms_tab_G3 %>%
  filter(P.Value.G3.Vs.others_DA < 1e-5  & logFC.G3.Vs.others_DA > 0.1)

draw.categoryValue(ac_val=ac_mat['LRP2_SIG',],
                   exp_val=exprs(CTBN_eset)[ms_tab_G3['LRP2_SIG','originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G3['LRP2','gene_label'],
                   main_exp=ms_tab_G3['LRP2','geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file=sprintf('surfaceome_proteins/G3_protein_boxplot_LRP2_CTBN.pdf', comp_name))

draw.categoryValue(ac_val=ac_mat['PODXL_SIG',],
                   exp_val=exprs(CTBN_eset)[ms_tab_G3['PODXL_SIG','originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G3['PODXL','gene_label'],
                   main_exp=ms_tab_G3['PODXL','geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file=sprintf('surfaceome_proteins/G3_protein_boxplot_PODXL_CTBN.pdf', comp_name))

draw.categoryValue(ac_val=ac_mat['TGFBR3_SIG',],
                   exp_val=exprs(CTBN_eset)[ms_tab_G3['TGFBR3_SIG','originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G3['TGFBR3','gene_label'],
                   main_exp=ms_tab_G3['TGFBR3','geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file=sprintf('surfaceome_proteins/G3_protein_boxplot_TGFBR3_CTBN.pdf', comp_name))

draw.categoryValue(ac_val=ac_mat['JAG2_SIG',],
                   exp_val=exprs(CTBN_eset)[ms_tab_G3['JAG2_SIG','originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G3['JAG2','gene_label'],
                   main_exp=ms_tab_G3['JAG2','geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file=sprintf('surfaceome_proteins/G3_protein_boxplot_JAG2_CTBN.pdf', comp_name))

draw.categoryValue(ac_val=ac_mat['VTN_SIG',],
                   exp_val=exprs(CTBN_eset)[ms_tab_G3['VTN_SIG','originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G3['VTN','gene_label'],
                   main_exp=ms_tab_G3['VTN','geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file=sprintf('surfaceome_proteins/G3_protein_boxplot_VTN_CTBN.pdf', comp_name))

#G4

# Create empty list to store comparison result
DE_G4 <- list()
DA_G4 <- list()

comp_name <- 'G4.Vs.others'
# Get sample names from each compared group
phe_info <- pData(CTBN_eset)
G1  <- rownames(phe_info)[which(phe_info$`Subgroup`=='G4')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`Subgroup`!='G4')] # Combine other groups as the Control group
DE_gene_bid_G4 <- getDE.BID.2G(eset=CTBN_eset,G1=G1,G0=G0,G1_name='G4',G0_name='others')
DA_driver_bid_G4   <- getDE.BID.2G(eset=merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='others')
# Save comparison result to list element in analysis.par, with comparison name
DE_G4[[comp_name]] <- DE_gene_bid_G4
DE_all[[comp_name]] <- DE_gene_bid_G4
DA_G4[[comp_name]] <- DA_driver_bid_G4
DA_all[[comp_name]] <- DA_driver_bid_G4

# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE)
# Get all comparison names
all_comp <- names(DE_G4) # Users can use index or name to get target ones

# Create the final master table
ms_tab_G4 <- generate.masterTable(use_comp=all_comp,DE=DE_G4,DA=DA_G4,
                               target_list=network$target_list,
                               tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                               main_id_type='external_gene_name')

# Visualize top genes
sig_gene_G4 <- draw.volcanoPlot(dat=ms_tab_G4,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
                                Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=1.5,Pv_thre=1e-5,
                                main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=TRUE,
                                pdf_file='volcano_plot_CTBN_G4.pdf',label_cex = 0.5)

gs.preload(use_spe='Homo sapiens',update=FALSE)
# Gene Set Enrichment Analysis
# gene_list_up <- rownames(sig_gene_G4)[which(sig_gene_G4[,2]>0)] # up
# gene_list_down <- rownames(sig_gene_G4)[which(sig_gene_G4[,2]<0)] # down
# res_up <- funcEnrich.Fisher(input_list=ms_tab_G4[gene_list_up,'geneSymbol'],bg_list=unique(ms_tab_G4[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
# res_down <- funcEnrich.Fisher(input_list=ms_tab_G4[gene_list_down,'geneSymbol'],bg_list=unique(ms_tab_G4[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       
# 
# # Gene set enrichment analysis Barplot
# draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top genes in G4 Subgroup',display_genes = TRUE,gs_cex=0.6,
#                     pdf_file='funcEnrich_bar_G4_CTBN.pdf')

filtered_ms_tab_G4 <- ms_tab_G4 %>%
  filter(P.Value.G4.Vs.others_DE < 1e-5  & logFC.G4.Vs.others_DE > 1.5)

# Visualize top drivers
sig_driver_G4 <- draw.volcanoPlot(dat=ms_tab_G4,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DA',comp_name),
                                  Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.1,Pv_thre=1e-15,
                                  main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,
                                  pdf_file='volcano_plot_CTBN_G4_DA.pdf',label_cex = 0.5)

gs.preload(use_spe='Homo sapiens',update=FALSE)
# Gene Set Enrichment Analysis
driver_list_up <- rownames(sig_driver_G4)[which(sig_driver_G4[,2]>0)] # up
driver_list_down <- rownames(sig_driver_G4)[which(sig_driver_G4[,2]<0)] # down
res_up <- funcEnrich.Fisher(input_list=ms_tab_G4[driver_list_up,'geneSymbol'],bg_list=unique(ms_tab_G4[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
res_down <- funcEnrich.Fisher(input_list=ms_tab_G4[driver_list_down,'geneSymbol'],bg_list=unique(ms_tab_G4[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       

# Gene set enrichment analysis Barplot
draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers in G4 Subgroup',display_genes = TRUE,gs_cex=0.6,
                    pdf_file='funcEnrich_bar_G4_DKFZ_WGS_DA.pdf')

filtered_ms_tab_G4_DA <- ms_tab_G4 %>%
  filter(P.Value.G4.Vs.others_DA < 1e-15  & logFC.G4.Vs.others_DA > 0.1)


#WNT
# filtered_ms_tab_WNT %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_WNT_2.csv', row.names = FALSE)
# filtered_ms_tab_WNT_DA %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_DA_WNT_2.csv', row.names = FALSE)
# #SHH
# filtered_ms_tab_SHH %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_SHH_2.csv', row.names = FALSE)
# filtered_ms_tab_SHH_DA %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_DA_SHH_2.csv', row.names = FALSE)
# #G3
# filtered_ms_tab_G3 %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_G3_2.csv', row.names = FALSE)
# filtered_ms_tab_G3_DA %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_DA_G3_2.csv', row.names = FALSE)
# #G4
# filtered_ms_tab_G4 %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_G4_2.csv', row.names = FALSE)
# filtered_ms_tab_G4_DA %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_DA_G4_2.csv', row.names = FALSE)

#GENERATING HEATMAP

# Combine the comparison results
# DE_all <- list()
# DA_all <- list()
# DE_list <- list(WNT=DE_gene_bid_WNT,SHH=DE_gene_bid_SHH, G3=DE_gene_bid_G3, G4=DE_gene_bid_G4)
# DE_gene_comb <- combineDE(DE_list,transfer_tab=NULL)
# DA_list <- list(WNT=DA_driver_bid_WNT,SHH=DA_driver_bid_SHH, G3=DA_driver_bid_G3, G4=DA_driver_bid_G4)
# DA_driver_comb <- combineDE(DA_list,transfer_tab=NULL)
# DE_all2[[comp_name]] <- DE_gene_comb$combine
# DA_all2[[comp_name]] <- DA_driver_comb$combine


db.preload(use_level='gene',use_spe='human',update=FALSE)
all_comp <- names(DE_all)

# Create the final master table
ms_tab_final <- generate.masterTable(use_comp=all_comp,DE=DE_all,DA=DA_all,
                                     target_list=network$target_list,
                                     tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                                     main_id_type='external_gene_name')

filtered_ms_tab_final <- ms_tab_final  %>%
  filter(((Z.WNT.Vs.others_DE>=2 | Z.WNT.Vs.others_DE<=-2) & (Z.WNT.Vs.others_DA>=2 | Z.WNT.Vs.others_DA<=-2)) | ((Z.SHH.Vs.others_DE>=2 | Z.SHH.Vs.others_DE<=-2) & (Z.SHH.Vs.others_DA>=2 | Z.SHH.Vs.others_DA<=-2)) | ((Z.G3.Vs.others_DE>=2 | Z.G3.Vs.others_DE<=-2) & (Z.G3.Vs.others_DA>=2 | Z.G3.Vs.others_DA<=-2)) | ((Z.G4.Vs.others_DE>=2 | Z.G4.Vs.others_DE<=-2) & (Z.G4.Vs.others_DA>=2 | Z.G4.Vs.others_DA<=-2)))

#reorder samples
subgroup_order <- c('WNT', 'SHH', 'G3', 'G4')
phe_info <- pData(CTBN_eset)
phe_info$Cohort <- 'CTBN'
phe_info$Subgroup <- factor(phe_info$Subgroup, levels = subgroup_order)
phe_info_sorted <- phe_info %>% arrange(Subgroup)

#draw heatmap
exp_mat <- exprs(cal.eset) # expression matrix, the rownames must be the originalID
sorted_sample_ids <- rownames(phe_info_sorted)
exp_mat_sorted <- exp_mat[, sorted_sample_ids]
ac_mat <- exprs(merge.ac.eset) # activity matrix, the rownames must be the originalID_label

#draw.heatmap(mat=exp_mat,use_genes=ms_tab_final$originalID,use_gene_label=ms_tab_final$geneSymbol,
# use_samples=colnames(exp_mat),
# phenotype_info=phe_info,use_phe=c('GENDER','SUBGROUP','COHORT'),main='Expression for Top drivers',scale='row',
# cluster_rows=TRUE,cluster_columns=TRUE,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
# row_names_gp = gpar(fontsize = 12),pdf_file='heatmap_expression.pdf',
# pre_define=c('WNT'='blue','SHH'='red','G3' = 'yellow','G4'='green'))

draw.heatmap(mat=exp_mat_sorted,use_genes=filtered_ms_tab_final$originalID,use_gene_label=ms_tab_final$geneSymbol,
             use_samples=colnames(exp_mat_sorted),
             phenotype_info=phe_info,use_phe=c('gender','Subgroup','Cohort'),main='Expression for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=FALSE, show_row_names = FALSE, clustering_distance_rows='spearman',
             row_names_gp = gpar(fontsize = 12),pdf_file=NULL,
             pre_define=c('WNT'='blue','SHH'='red','G3' = 'yellow','G4'='green'))
#'heatmap_expression_CTBN.pdf'

draw.heatmap(mat=exp_mat_sorted,use_genes=ms_tab_final$originalID,use_gene_label=ms_tab_final$geneSymbol,
             use_samples=colnames(exp_mat_sorted),
             phenotype_info=phe_info,use_phe=c('gender','Subgroup','Cohort'),main='Activity for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=FALSE,show_row_names = FALSE, clustering_distance_rows='spearman',
             row_names_gp = gpar(fontsize = 12),pdf_file=NULL,
             pre_define=c('WNT'='blue','SHH'='red','G3' = 'yellow','G4'='green'))
#'heatmap_activity_CTBN.pdf'
