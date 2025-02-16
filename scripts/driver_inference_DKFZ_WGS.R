library(NetBID2);
library(dplyr);
library(forcats);
setwd("/Users/abbiesiru/Desktop/research/medulloblastoma/")

############### Step 1: Load in gene expression dataset for analysis (exp-load, exp-cluster, exp-QC) ###############

cal.eset <- readRDS("dataset/northcott_N169_geneLevel_eset.rds")

############### Step 2: Read in network files and calcualte driver activity (act-get) ###############

network <- readRDS("network/GSE85217_allMB_network.rds")

# filter lowly expressed genes
# Get activity matrix
ac_mat <- cal.Activity(target_list=network$target_list,cal_mat=exprs(cal.eset),es.method='weightedmean')

# Create eset using activity matrix
merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(cal.eset)[colnames(ac_mat),],
                               feature_info=NULL,annotation_info='activity in net-dataset')

# QC plot for activity eset
#draw.eset.QC(merge.ac.eset,outdir="network/",intgroup=NULL,do.logtransform=FALSE,prefix='DKFZ_WGS_AC_',emb_plot_type ='2D.interactive')


############### Step 3: Get differential expression (DE) / differential activity (DA) for drivers (act-DA) ###############

#WNT

# Create empty list to store comparison result
DE_WNT <- list()
DA_WNT <- list()
DE_all <- list()
DA_all <- list()

comp_name <- 'WNT.Vs.others'
# Get sample names from each compared group
phe_info <- pData(cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`Subgroup`=='WNT')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`Subgroup`!='WNT')] # Combine other groups as the Control group
DE_gene_bid_WNT <- getDE.BID.2G(eset=cal.eset,G1=G1,G0=G0,G1_name='WNT',G0_name='others')
DA_driver_bid_WNT   <- getDE.BID.2G(eset=merge.ac.eset,G1=G1,G0=G0,G1_name='WNT',G0_name='others')
# Save comparison result to list element in analysis.par, with comparison name
DE_WNT[[comp_name]] <- DE_gene_bid_WNT
DE_all[[comp_name]] <- DE_gene_bid_WNT
DA_WNT[[comp_name]] <- DA_driver_bid_WNT
DA_all[[comp_name]] <- DA_driver_bid_WNT


# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE)
# Get all comparison names
all_comp <- names(DE_WNT) 

# Create the final master table
ms_tab_WNT <- generate.masterTable(use_comp=all_comp,DE=DE_WNT,DA=DA_WNT,
                                   target_list=network$target_list,
                                   tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                                   main_id_type='external_gene_name')

# Visualize top genes
sig_gene_WNT <- draw.volcanoPlot(dat=ms_tab_WNT,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
                                 Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=2,Pv_thre=1e-5,
                                 main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=TRUE,
                                 pdf_file=NULL,label_cex = 0.5)

gs.preload(use_spe='Homo sapiens',update=FALSE)
# Gene Set Enrichment Analysis
#gene_list_up <- rownames(sig_gene_WNT)[which(sig_gene_WNT[,2]>0)] # up
#gene_list_down <- rownames(sig_gene_WNT)[which(sig_gene_WNT[,2]<0)] # down
#res_up <- funcEnrich.Fisher(input_list=ms_tab_WNT[gene_list_up,'geneSymbol'],bg_list=unique(ms_tab_WNT[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
#res_down <- funcEnrich.Fisher(input_list=ms_tab_WNT[gene_list_down,'geneSymbol'],bg_list=unique(ms_tab_WNT[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       

# Gene set enrichment analysis Barplot
#draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top genes in WNT Subgroup',display_genes = TRUE,gs_cex=0.6,
                    #pdf_file='funcEnrich_bar_WNT_DKFZ_WGS.pdf')

filtered_ms_tab_WNT <- ms_tab_WNT %>%
  filter(P.Value.WNT.Vs.others_DE < 1e-5  & logFC.WNT.Vs.others_DE > 2)

# Visualize top drivers
#sig_driver_WNT <- draw.volcanoPlot(dat=ms_tab_WNT,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DA',comp_name),
                                 #Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.15,Pv_thre=1e-15,
                                 #main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,
                                 #pdf_file=NULL,label_cex = 0.5)

#gs.preload(use_spe='Homo sapiens',update=FALSE)
# Gene Set Enrichment Analysis
#driver_list_up <- rownames(sig_driver_WNT)[which(sig_driver_WNT[,2]>0)] # up
#driver_list_down <- rownames(sig_driver_WNT)[which(sig_driver_WNT[,2]<0)] # down
#res_up <- funcEnrich.Fisher(input_list=ms_tab_WNT[driver_list_up,'geneSymbol'],bg_list=unique(ms_tab_WNT[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
#res_down <- funcEnrich.Fisher(input_list=ms_tab_WNT[driver_list_down,'geneSymbol'],bg_list=unique(ms_tab_WNT[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       

# Gene set enrichment analysis Barplot
#draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers in WNT Subgroup',display_genes = TRUE,gs_cex=0.6,
                    #pdf_file='funcEnrich_bar_WNT_DKFZ_WGS_DA.pdf')

filtered_ms_tab_WNT_DA <- ms_tab_WNT %>%
  filter(P.Value.WNT.Vs.others_DA < 1e-15  & logFC.WNT.Vs.others_DA > 0.15)

draw.categoryValue(ac_val=ac_mat[7767,],
                   exp_val=exprs(cal.eset)[ms_tab_WNT[395,'originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_WNT[395,'gene_label'],
                   main_exp=ms_tab_WNT[395,'geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'))
#pdf_file=sprintf('surfaceome_proteins/WNT_protein_boxplot.pdf', comp_name)
#SHH

# Create empty list to store comparison result
DE_SHH <- list()
DA_SHH <- list()

comp_name <- 'SHH.Vs.others'
# Get sample names from each compared group
phe_info <- pData(cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`Subgroup`=='SHH')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`Subgroup`!='SHH')] # Combine other groups as the Control group
DE_gene_bid_SHH <- getDE.BID.2G(eset=cal.eset,G1=G1,G0=G0,G1_name='SHH',G0_name='others')
DA_driver_bid_SHH   <- getDE.BID.2G(eset=merge.ac.eset,G1=G1,G0=G0,G1_name='SHH',G0_name='others')
# Save comparison result to list element in analysis.par, with comparison name
DE_SHH[[comp_name]] <- DE_gene_bid_SHH
DE_all[[comp_name]] <- DE_gene_bid_SHH
DA_SHH[[comp_name]] <- DA_driver_bid_SHH
DA_all[[comp_name]] <- DA_driver_bid_SHH

# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE)
# Get all comparison names
all_comp <- names(DE_SHH) 

# Create the final master table
ms_tab_SHH <- generate.masterTable(use_comp=all_comp,DE=DE_SHH,DA=DA_SHH,
                               target_list=network$target_list,
                               tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                               main_id_type='external_gene_name')

# Visualize top genes
sig_gene_SHH <- draw.volcanoPlot(dat=ms_tab_SHH,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
                                Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=2,Pv_thre=1e-5,
                                main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=TRUE,
                                pdf_file=NULL,label_cex = 0.5)

gs.preload(use_spe='Homo sapiens',update=FALSE)
# Gene Set Enrichment Analysis
#gene_list_up <- rownames(sig_gene_SHH)[which(sig_gene_SHH[,2]>0)] # up
#gene_list_down <- rownames(sig_gene_SHH)[which(sig_gene_SHH[,2]<0)] # down
#res_up <- funcEnrich.Fisher(input_list=ms_tab_SHH[gene_list_up,'geneSymbol'],bg_list=unique(ms_tab_SHH[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
#es_down <- funcEnrich.Fisher(input_list=ms_tab_SHH[gene_list_down,'geneSymbol'],bg_list=unique(ms_tab_SHH[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       

# Gene set enrichment analysis Barplot
#draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top genes in SHH Subgroup',display_genes = TRUE,gs_cex=0.6,
                    #pdf_file='funcEnrich_bar_SHH_DKFZ_WGS.pdf')

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
#draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers in SHH Subgroup',display_genes = TRUE,gs_cex=0.6,
                    #pdf_file='funcEnrich_bar_SHH_DKFZ_WGS_DA.pdf')

filtered_ms_tab_SHH_DA <- ms_tab_SHH %>%
  filter(P.Value.SHH.Vs.others_DA < 1e-15  & logFC.SHH.Vs.others_DA > 0.15)

draw.categoryValue(ac_val=ac_mat[7348,],
                   exp_val=exprs(cal.eset)[ms_tab_SHH[91,'originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_SHH[91,'gene_label'],
                   main_exp=ms_tab_SHH[91,'geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'))
#pdf_file=sprintf('surfaceome_proteins/SHH_protein_boxplot.pdf', comp_name)
#G3

# Create empty list to store comparison result
DE_G3 <- list()
DA_G3 <- list()

comp_name <- 'G3.Vs.others'
# Get sample names from each compared group
phe_info <- pData(cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`Subgroup`=='G3')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`Subgroup`!='G3')] # Combine other groups as the Control group
DE_gene_bid_G3 <- getDE.BID.2G(eset=cal.eset,G1=G1,G0=G0,G1_name='G3',G0_name='others')
DA_driver_bid_G3   <- getDE.BID.2G(eset=merge.ac.eset,G1=G1,G0=G0,G1_name='G3',G0_name='others')
# Save comparison result to list element in analysis.par, with comparison name
DE_G3[[comp_name]] <- DE_gene_bid_G3
#DE_all[[comp_name]] <- DE_gene_bid_G3
DA_G3[[comp_name]] <- DA_driver_bid_G3
#DA_all[[comp_name]] <- DA_driver_bid_G3

# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE)
# Get all comparison names
all_comp <- names(DE_G3)

# Create the final master table
ms_tab_G3 <- generate.masterTable(use_comp=all_comp,DE=DE_G3,DA=DA_G3,
                               target_list=network$target_list,
                               tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                               main_id_type='external_gene_name')

# Visualize top genes
sig_gene_G3 <- draw.volcanoPlot(dat=ms_tab_G3,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
                             Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=2,Pv_thre=1e-5,
                             main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=TRUE,
                             pdf_file=NULL,label_cex = 0.5)
#logFC_thre=1; Pv

gs.preload(use_spe='Homo sapiens',update=FALSE)
# Gene Set Enrichment Analysis
#gene_list_up <- rownames(sig_gene_G3)[which(sig_gene_G3[,2]>0)] # up
#gene_list_down <- rownames(sig_gene_G3)[which(sig_gene_G3[,2]<0)] # down
#res_up <- funcEnrich.Fisher(input_list=ms_tab_G3[gene_list_up,'geneSymbol'],bg_list=unique(ms_tab_G3[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
#res_down <- funcEnrich.Fisher(input_list=ms_tab_G3[gene_list_down,'geneSymbol'],bg_list=unique(ms_tab_G3[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       

# Gene set enrichment analysis Barplot
#draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top genes in G3 Subgroup',display_genes = TRUE,gs_cex=0.6,
                    #pdf_file='funcEnrich_bar_G3_DKFZ_WGS.pdf')

filtered_ms_tab_G3 <- ms_tab_G3 %>%
  filter(P.Value.G3.Vs.others_DE < 1e-5  & logFC.G3.Vs.others_DE > 2)


# Visualize top drivers
sig_driver_G3 <- draw.volcanoPlot(dat=ms_tab_G3,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DA',comp_name),
                                   Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.1,Pv_thre=1e-5,
                                   main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,
                                   pdf_file=NULL,label_cex = 0.5)

gs.preload(use_spe='Homo sapiens',update=FALSE)
# Gene Set Enrichment Analysis
#driver_list_up <- rownames(sig_driver_G3)[which(sig_driver_G3[,2]>0)] # up
#driver_list_down <- rownames(sig_driver_G3)[which(sig_driver_G3[,2]<0)] # down
#res_up <- funcEnrich.Fisher(input_list=ms_tab_G3[driver_list_up,'geneSymbol'],bg_list=unique(ms_tab_G3[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
#res_down <- funcEnrich.Fisher(input_list=ms_tab_G3[driver_list_down,'geneSymbol'],bg_list=unique(ms_tab_G3[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       

# Gene set enrichment analysis Barplot
#draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers in G3 Subgroup',display_genes = TRUE,gs_cex=0.6,
                    #pdf_file='funcEnrich_bar_G3_DKFZ_WGS_DA.pdf')

filtered_ms_tab_G3_DA <- ms_tab_G3 %>%
  filter(P.Value.G3.Vs.others_DA < 1e-5  & logFC.G3.Vs.others_DA > 0.1)

draw.categoryValue(ac_val=ac_mat['LRP2_SIG',],
                   exp_val=exprs(cal.eset)[ms_tab_G3['LRP2_SIG','originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G3['LRP2_SIG','gene_label'],
                   main_exp=ms_tab_G3['LRP2_SIG','geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file=sprintf('surfaceome_proteins/G3_protein_boxplot_LRP2_DKFZ_WGS.pdf', comp_name))

draw.categoryValue(ac_val=ac_mat['PODXL_SIG',],
                   exp_val=exprs(cal.eset)[ms_tab_G3['PODXL_SIG','originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G3['PODXL_SIG','gene_label'],
                   main_exp=ms_tab_G3['PODXL_SIG','geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file=sprintf('surfaceome_proteins/G3_protein_boxplot_PODXL_DKFZ_WGS.pdf', comp_name))

draw.categoryValue(ac_val=ac_mat['TGFBR3_SIG',],
                   exp_val=exprs(cal.eset)[ms_tab_G3['TGFBR3_SIG','originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G3['TGFBR3_SIG','gene_label'],
                   main_exp=ms_tab_G3['TGFBR3_SIG','geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file=sprintf('surfaceome_proteins/G3_protein_boxplot_TGFBR3_DKFZ_WGS.pdf', comp_name))

draw.categoryValue(ac_val=ac_mat['JAG2_SIG',],
                   exp_val=exprs(cal.eset)[ms_tab_G3['JAG2_SIG','originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G3['JAG2_SIG','gene_label'],
                   main_exp=ms_tab_G3['JAG2_SIG','geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file=sprintf('surfaceome_proteins/G3_protein_boxplot_JAG2_DKFZ_WGS.pdf', comp_name))

draw.categoryValue(ac_val=ac_mat['VTN_SIG',],
                   exp_val=exprs(cal.eset)[ms_tab_G3['VTN_SIG','originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G3['VTN_SIG','gene_label'],
                   main_exp=ms_tab_G3['VTN_SIG','geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file=sprintf('surfaceome_proteins/G3_protein_boxplot_VTN_DKFZ_WGS.pdf', comp_name))

#G4

# Create empty list to store comparison result
DE_G4 <- list()
DA_G4 <- list()

comp_name <- 'G4.Vs.others'
# Get sample names from each compared group
phe_info <- pData(cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`Subgroup`=='G4')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`Subgroup`!='G4')] # Combine other groups as the Control group
DE_gene_bid_G4 <- getDE.BID.2G(eset=cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='others')
DA_driver_bid_G4   <- getDE.BID.2G(eset=merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='others')
# Save comparison result to list element in analysis.par, with comparison name
DE_G4[[comp_name]] <- DE_gene_bid_G4
DE_all[[comp_name]] <- DE_gene_bid_G4
DA_G4[[comp_name]] <- DA_driver_bid_G4
DA_all[[comp_name]] <- DA_driver_bid_G4

# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE)
# Get all comparison names
all_comp <- names(DE_G4)

# Create the final master table
ms_tab_G4 <- generate.masterTable(use_comp=all_comp,DE=DE_G4,DA=DA_G4,
                               target_list=network$target_list,
                               tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                               main_id_type='external_gene_name')

# Visualize top genes
sig_gene_G4 <- draw.volcanoPlot(dat=ms_tab_G4,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
                                 Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=1.5,Pv_thre=1e-5,
                                 main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=TRUE,
                                 pdf_file=NULL,label_cex = 0.5)

gs.preload(use_spe='Homo sapiens',update=FALSE)
# Gene Set Enrichment Analysis
#gene_list_up <- rownames(sig_gene_G4)[which(sig_gene_G4[,2]>0)] # up
#gene_list_down <- rownames(sig_gene_G4)[which(sig_gene_G4[,2]<0)] # down
#res_up <- funcEnrich.Fisher(input_list=ms_tab_G4[gene_list_up,'geneSymbol'],bg_list=unique(ms_tab_G4[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
#res_down <- funcEnrich.Fisher(input_list=ms_tab_G4[gene_list_down,'geneSymbol'],bg_list=unique(ms_tab_G4[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       

# Gene set enrichment analysis Barplot
#draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top genes in G4 Subgroup',display_genes = TRUE,gs_cex=0.6,
                    #pdf_file='funcEnrich_bar_G4_DKFZ_WGS.pdf')

filtered_ms_tab_G4 <- ms_tab_G4 %>%
  filter(P.Value.G4.Vs.others_DE < 1e-5  & logFC.G4.Vs.others_DE > 1.5)

# Visualize top drivers
sig_driver_G4 <- draw.volcanoPlot(dat=ms_tab_G4,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DA',comp_name),
                                   Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.1,Pv_thre=1e-15,
                                   main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,
                                   pdf_file=NULL,label_cex = 0.5)

gs.preload(use_spe='Homo sapiens',update=FALSE)
# Gene Set Enrichment Analysis
#driver_list_up <- rownames(sig_driver_G4)[which(sig_driver_G4[,2]>0)] # up
#driver_list_down <- rownames(sig_driver_G4)[which(sig_driver_G4[,2]<0)] # down
#res_up <- funcEnrich.Fisher(input_list=ms_tab_G4[driver_list_up,'geneSymbol'],bg_list=unique(ms_tab_G4[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
#res_down <- funcEnrich.Fisher(input_list=ms_tab_G4[driver_list_down,'geneSymbol'],bg_list=unique(ms_tab_G4[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       

# Gene set enrichment analysis Barplot
#draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers in G4 Subgroup',display_genes = TRUE,gs_cex=0.6,
                    #pdf_file='funcEnrich_bar_G4_DKFZ_WGS_DA.pdf')

filtered_ms_tab_G4_DA <- ms_tab_G4 %>%
  filter(P.Value.G4.Vs.others_DA < 1e-15  & logFC.G4.Vs.others_DA > 0.1)

draw.categoryValue(ac_val=ac_mat[3831,],
                   exp_val=exprs(cal.eset)[ms_tab_G4[58,'originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G4[58,'gene_label'],
                   main_exp=ms_tab_G4[58,'geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file=sprintf('surfaceome_proteins/G4_protein_boxplot.pdf', comp_name))


#WNT
# filtered_ms_tab_WNT %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_WNT_1.csv', row.names = FALSE)
# filtered_ms_tab_WNT_DA %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_DA_WNT_1.csv', row.names = FALSE)
# #SHH
# filtered_ms_tab_SHH %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_SHH_1.csv', row.names = FALSE)
# filtered_ms_tab_SHH_DA %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_DA_SHH_1.csv', row.names = FALSE)
# #G3
# filtered_ms_tab_G3 %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_G3_1.csv', row.names = FALSE)
# filtered_ms_tab_G3_DA %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_DA_G3_1.csv', row.names = FALSE)
# #G4
# filtered_ms_tab_G4 %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_G4_1.csv', row.names = FALSE)
# filtered_ms_tab_G4_DA %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_DA_G4_1.csv', row.names = FALSE)

#GENERATING HEATMAP

# Combine the comparison results
# DE_all <- list()
# DA_all <- list()
# DE_list <- list(WNT=DE_gene_bid_WNT,SHH=DE_gene_bid_SHH, G3=DE_gene_bid_G3, G4=DE_gene_bid_G4)
# DE_gene_comb <- combineDE(DE_list,transfer_tab=NULL)
# DA_list <- list(WNT=DA_driver_bid_WNT,SHH=DA_driver_bid_SHH, G3=DA_driver_bid_G3, G4=DA_driver_bid_G4)
# DA_driver_comb <- combineDE(DA_list,transfer_tab=NULL)
# DE_all[[comp_name]] <- DE_gene_comb$combine
# DA_all[[comp_name]] <- DA_driver_comb$combine


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
subgroup_order <- c('WNT', 'SHH', 'Group3', 'Group4')
phe_info$SUBGROUP <- factor(phe_info$SUBGROUP, levels = subgroup_order)
phe_info_sorted <- phe_info %>% arrange(SUBGROUP)

#draw heatmap
exp_mat <- exprs(cal.eset) # expression matrix, the rownames must be the originalID
sorted_sample_ids <- phe_info_sorted$sampleID
exp_mat_sorted <- exp_mat[, sorted_sample_ids]
ac_mat <- exprs(merge.ac.eset) # activity matrix, the rownames must be the originalID_label

#draw.heatmap(mat=exp_mat,use_genes=ms_tab_final$originalID,use_gene_label=ms_tab_final$geneSymbol,
            # use_samples=colnames(exp_mat),
            # phenotype_info=phe_info,use_phe=c('GENDER','SUBGROUP','COHORT'),main='Expression for Top drivers',scale='row',
            # cluster_rows=TRUE,cluster_columns=TRUE,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
            # row_names_gp = gpar(fontsize = 12),pdf_file='heatmap_expression.pdf',
            # pre_define=c('WNT'='blue','SHH'='red','G3' = 'yellow','G4'='green'))

draw.heatmap(mat=exp_mat_sorted,use_genes=ms_tab_final$originalID,use_gene_label=ms_tab_final$geneSymbol,
             use_samples=colnames(exp_mat_sorted),
             phenotype_info=phe_info,use_phe=c('GENDER','SUBGROUP','COHORT'),main='Expression for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=FALSE, show_row_names = FALSE, clustering_distance_rows='spearman',
             row_names_gp = gpar(fontsize = 12),pdf_file=NULL,
             pre_define=c('WNT'='blue','SHH'='red','Group3' = 'yellow','Group4'='green'))
#'heatmap_expression_DKFZ_WGS.pdf'

draw.heatmap(mat=exp_mat_sorted,use_genes=ms_tab_final$originalID,use_gene_label=ms_tab_final$geneSymbol,
             use_samples=colnames(exp_mat_sorted),
             phenotype_info=phe_info,use_phe=c('GENDER','SUBGROUP','COHORT'),main='Activity for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=FALSE,show_row_names = FALSE, clustering_distance_rows='spearman',
             row_names_gp = gpar(fontsize = 12),pdf_file=NULL,
             pre_define=c('WNT'='blue','SHH'='red','Group3' = 'yellow','Group4'='green'))
#'heatmap_activity_DKFZ_WGS.pdf'
