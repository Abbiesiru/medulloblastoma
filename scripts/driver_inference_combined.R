library(NetBID2);
library(dplyr);
library(forcats);
setwd("/Users/abbiesiru/Desktop/research/medulloblastoma/")

### load in gene expression dataset for analysis (exp-load, exp-cluster, exp-QC) ###

cal.eset <- readRDS("dataset/dataset_combined.rds")

### read in network files and calculate driver activity (act-get) ###

network <- readRDS("network/GSE85217_allMB_network.rds")
ac_mat <- cal.Activity(target_list=network$target_list,cal_mat=exprs(cal.eset),es.method='weightedmean')
merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(cal.eset)[colnames(ac_mat),],)


### get differential expression (DE) / differential activity (DA) for drivers (act-DA) ###

# 1. create empty list to store comparison result
DE_G3 <- list()
DA_G3 <- list()

# 2. get sample names from each compared group
comp_name <- 'G3.Vs.others'
phe_info <- pData(cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`Subgroup`=='G3')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`Subgroup`!='G3')] # Combine other groups as the Control group
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

# 6. visualize top genes
sig_gene_G3 <- draw.volcanoPlot(dat=ms_tab_G3,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
                                Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=1,Pv_thre=0.5,
                                main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=TRUE,
                                pdf_file="volcano_plot_combined_G3_DE.pdf",label_cex = 0.5)

gs.preload(use_spe='Homo sapiens',update=FALSE)

# 7. DE gene set enrichment analysis
gene_list_up <- rownames(sig_gene_G3)[which(sig_gene_G3[,2]>0)] # up
gene_list_down <- rownames(sig_gene_G3)[which(sig_gene_G3[,2]<0)] # down
res_up <- funcEnrich.Fisher(input_list=ms_tab_G3[gene_list_up,'geneSymbol'],bg_list=unique(ms_tab_G3[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
res_down <- funcEnrich.Fisher(input_list=ms_tab_G3[gene_list_down,'geneSymbol'],bg_list=unique(ms_tab_G3[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       

draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top genes in G3 Subgroup',display_genes = TRUE,gs_cex=0.6,
                    pdf_file="funcEnrich_bar_G3_combined_DE.pdf") # barplot


# 8. visualize top drivers
sig_driver_G3 <- draw.volcanoPlot(dat=ms_tab_G3,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DA',comp_name),
                                  Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.35,Pv_thre=0.5,
                                  main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,
                                  pdf_file="volcano_plot_combined_G3_DA.pdf",label_cex = 0.5)


# 8. DA gene set enrichment analysis
driver_list_up <- rownames(sig_driver_G3)[which(sig_driver_G3[,2]>0)] # up
driver_list_down <- rownames(sig_driver_G3)[which(sig_driver_G3[,2]<0)] # down
res_up <- funcEnrich.Fisher(input_list=ms_tab_G3[driver_list_up,'geneSymbol'],bg_list=unique(ms_tab_G3[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
res_down <- funcEnrich.Fisher(input_list=ms_tab_G3[driver_list_down,'geneSymbol'],bg_list=unique(ms_tab_G3[,'geneSymbol']),use_gs=NULL, Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)                       
draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers in G3 Subgroup',display_genes = TRUE,gs_cex=0.6,
                    pdf_file="funcEnrich_bar_G3_combined_DA.pdf") # barplot


### bar charts for potential target candidates ###

# 1. LRP2
draw.categoryValue(ac_val=ac_mat['LRP2_SIG',],
                   exp_val=exprs(cal.eset)[ms_tab_G3['LRP2_SIG','originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G3['LRP2_SIG','gene_label'],
                   main_exp=ms_tab_G3['LRP2_SIG','geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file="surfaceome_proteins/G3_protein_boxplot_LRP2_combined.pdf")

# 2. PODXL
draw.categoryValue(ac_val=ac_mat['PODXL_SIG',],
                   exp_val=exprs(cal.eset)[ms_tab_G3['PODXL_SIG','originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G3['PODXL_SIG','gene_label'],
                   main_exp=ms_tab_G3['PODXL_SIG','geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file="surfaceome_proteins/G3_protein_boxplot_PODXL_combined.pdf")

# 3. TGFBR3
draw.categoryValue(ac_val=ac_mat['TGFBR3_SIG',],
                   exp_val=exprs(cal.eset)[ms_tab_G3['TGFBR3_SIG','originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G3['TGFBR3_SIG','gene_label'],
                   main_exp=ms_tab_G3['TGFBR3_SIG','geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file="surfaceome_proteins/G3_protein_boxplot_TGFBR3_combined.pdf")

# 4. JAG2
draw.categoryValue(ac_val=ac_mat['JAG2_SIG',],
                   exp_val=exprs(cal.eset)[ms_tab_G3['JAG2_SIG','originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G3['JAG2_SIG','gene_label'],
                   main_exp=ms_tab_G3['JAG2_SIG','geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file="surfaceome_proteins/G3_protein_boxplot_JAG2_combined.pdf")

# 5. VTN
draw.categoryValue(ac_val=ac_mat['VTN_SIG',],
                   exp_val=exprs(cal.eset)[ms_tab_G3['VTN_SIG','originalID'],],
                   use_obs_class=get_obs_label(phe_info = phe_info,'Subgroup'),
                   class_order=c('WNT','SHH','G3', 'G4'),
                   class_srt=30,
                   main_ac = ms_tab_G3['VTN_SIG','gene_label'],
                   main_exp=ms_tab_G3['VTN_SIG','geneSymbol'],
                   pre_define=c('WNT'='blue','SHH'='red','G3'='yellow', 'G4'='green'), pdf_file="surfaceome_proteins/G3_protein_boxplot_VTN_combined.pdf")

### exporting DE/DA gene list ###
filtered_ms_tab_G3 <- ms_tab_G3 %>%
  filter(P.Value.G3.Vs.others_DE < 0.5  & logFC.G3.Vs.others_DE > 1)
filtered_ms_tab_G3_DA <- ms_tab_G3 %>%
  filter(P.Value.G3.Vs.others_DA < 0.5  & logFC.G3.Vs.others_DA > 1)
# filtered_ms_tab_G3 %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_G3_1.csv', row.names = FALSE)
# filtered_ms_tab_G3_DA %>%
#   select(geneSymbol) %>%
#   write.csv(file = 'geneSymbol_DA_G3_1.csv', row.names = FALSE)
