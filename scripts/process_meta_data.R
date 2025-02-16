library(dplyr)
library(stringr)
library(readxl)

setwd("/Users/abbiesiru/Desktop/research/medulloblastoma");

# Read meta data file into R as a data frame
meta_data <- read_xls("dataset/manifest_source-data_RNA-Seq.xls", col_names=TRUE);

meta_data_MB <- dplyr::filter(meta_data, disease_type == "Medulloblastoma")
# View(meta_data_MB)

# go through the sample_name column once to 
# meta_data_MB$sample_name = ifelse(grepl(".local.transcript.bam", meta_data_MB$sample_name),str_sub(meta_data_MB$sample_name, 13, -22), meta_data_MB$sample_name);
# go through the sample_name column a second time to 
# meta_data_MB$sample_name = ifelse(grepl(".bam", meta_data_MB$sample_name),str_sub(meta_data_MB$sample_name, 13, -5), meta_data_MB$sample_name);
# go through the sample_name column a third time to 
# meta_data_MB$sample_name = ifelse(grepl(".fq.gz", meta_data_MB$sample_name),str_sub(meta_data_MB$sample_name, 13, -9), meta_data_MB$sample_name);

# create a new column that starts as a copy of the "name" column
meta_data_MB$sample_name = meta_data_MB$name;
for (i in 1:length(meta_data_MB$sample_name)) {
  # remove all ".local.transcript.bam" postfix
  if (grepl(".local.transcript.bam", meta_data_MB$sample_name[i])) {
    meta_data_MB$sample_name[i] = str_sub(meta_data_MB$sample_name[i], 13, -22);
  }
  # remove all ".bam" postfix
  else if (grepl(".bam", meta_data_MB$sample_name[i])) {
    meta_data_MB$sample_name[i] = str_sub(meta_data_MB$sample_name[i], 13, -5);
  }
  # remove all "_[1|2].fq.gz" postfix
  else if (grepl(".fq.gz", meta_data_MB$sample_name[i])) {
    meta_data_MB$sample_name[i] = str_sub(meta_data_MB$sample_name[i], 13, -9);
  }
  else {
    print("This is impossible");
  }
}

# View(meta_data_MB)

meta_data_MB <- meta_data_MB %>%
  dplyr::relocate(sample_name,.before = id)
# View(meta_data_MB)

# Remove duplicates based on the new sample_name columns
meta_data_MB <- meta_data_MB[!duplicated(meta_data_MB$sample_name), ]

write.csv(meta_data_MB, "dataset/manifest_source-data_RNA-Seq_MB.csv",row.names=FALSE)

# write.xlsx(meta_data_MB, 'dataset/manifest_source-data_RNA-Seq_MB.xlsx')
