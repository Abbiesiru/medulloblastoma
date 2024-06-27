install.packages("dplyr")
library(dplyr)
install.packages("stringr")
library(stringr)

library(readxl)

# Read meta data file into R as a dataframe
meta_data <- read_xls("manifest_source-data_RNA-Seq.xls", col_names=TRUE);

meta_data_MB <- dplyr::filter(meta_data, disease_type == "Medulloblastoma")
View(meta_data_MB)


meta_data_MB$sample_name <- ifelse(grepl(".local.transcript.bam", meta_data_MB$name),str_sub(meta_data_MB$name, 13, -22), str_sub(meta_data_MB$name, 13, -5))
meta_data_MB$sample_name <- ifelse(grepl(".fq.gz", meta_data_MB$name),str_sub(meta_data_MB$name, 13, -7), meta_data_MB$sample_name)
View(meta_data_MB)

meta_data_MB <- meta_data_MB %>%
  dplyr::relocate(sample_name,.after = name)
View(meta_data_MB)

install.packages("openxlsx")
library(openxlsx)
write.xlsx(meta_data_MB, 'manifest_source-data_RNA-Seq_MB.xlsx')
