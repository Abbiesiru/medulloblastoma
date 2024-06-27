install.packages("dplyr")
library(dplyr)
install.packages("stringr")
library(stringr)

manifest_source_data_RNA_Seq <- data.read("");

manifest_source_data_RNA_Seq_MB <- dplyr::filter(manifest_source_data_RNA_Seq, disease_type == "Medulloblastoma")
View(manifest_source_data_RNA-Seq_MB)


manifest_source_data_RNA_Seq_MB$sample_name <- ifelse(grepl(".local.transcript.bam", manifest_source_data_RNA_Seq_MB$name),str_sub(manifest_source_data_RNA_Seq_MB$name, 13, -22), str_sub(manifest_source_data_RNA_Seq_MB$name, 13, -5))
manifest_source_data_RNA_Seq_MB$sample_name <- ifelse(grepl(".fq.gz", manifest_source_data_RNA_Seq_MB$name),str_sub(manifest_source_data_RNA_Seq_MB$name, 13, -7), manifest_source_data_RNA_Seq_MB$sample_name)
View(manifest_source_data_RNA_Seq_MB)

manifest_source_data_RNA_Seq_MB <- manifest_source_data_RNA_Seq_MB %>%
  dplyr::relocate(sample_name,.after = name)
View(manifest_source_data_RNA_Seq_MB)

install.packages("openxlsx")
library(openxlsx)
write.xlsx(manifest_source_data_RNA_Seq_MB, 'manifest_source-data_RNA-Seq_MB.xlsx')
