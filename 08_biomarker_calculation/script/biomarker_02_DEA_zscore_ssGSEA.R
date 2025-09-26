setwd('~/Documents/projects/project_DGD/immunotherapy_biomarkers/')

library(tidyverse)
library(glue)
library(biomaRt)
library(GSVA)

source('08_biomarker_calculation/script/functions.R')

# Load data
df_clinical <- readRDS('01_metadata/out/00_metadata.rds')
geneset <- readRDS('05_genesets/out/geneset.rds')

# dea_files <- NULL
# 
# # All samples
# samples <- df_clinical %>% dplyr::select(Sample) %>% unlist()
# 
# for (sample in samples){
# 
#   # sample <- samples[[2]]
# 
#   # Load DEA file from 05_DEA.py output
#   dea_file <- suppressMessages(
#     read_csv(glue('04_DEA/out/DEA_{sample}.csv'), show_col_types = FALSE) %>% rename('genes' = '...1')
#   )
# 
#   dea_file <- dea_file %>%
#     mutate(z_score = sign(log2_fold_change) * qnorm(1 - p_value/2))
# 
#   colnames(dea_file) <- c('genes', glue("p_value__{sample}"), glue("q_value__{sample}"),
#                           glue("log2_fold_change__{sample}"), glue("z_score__{sample}"))
# 
#   # dea_file_out <- dea_file %>% dplyr::select(genes, z_score) %>% rename(!!glue("z_score_{sample}") := z_score)
# 
#   # Merge dea files from all samples
#   if (is.null(dea_files)){
#     dea_files <- dea_file
#   } else {
#     dea_files <- dea_files %>% left_join(dea_file, by = 'genes')
#   }
# 
# }
# 
# saveRDS(dea_files, '08_biomarker_calculation/out/DEA_files.rds')
dea_files <- readRDS('08_biomarker_calculation/out/DEA_files.rds')

# Contains score for all samples 
dea_files_mat <- dea_files %>% dplyr::select(genes, starts_with('z_score')) %>% column_to_rownames('genes') %>% as.matrix()

# Change Inf to large value
dea_files_mat[is.infinite(dea_files_mat) & dea_files_mat > 0] <- 10   # replace Inf
dea_files_mat[is.infinite(dea_files_mat) & dea_files_mat < 0] <- -10  # replace -Inf

# Run ssGSEA and save output
ssgsea_on_zscore <- run_ssgsea(dea_files_mat, geneset) 

# Reformat ssGSEA output
ssgsea_on_zscore <- ssgsea_on_zscore %>% t() %>% as.data.frame() %>% rownames_to_column('Sample')
ssgsea_on_zscore$Sample <- str_remove(ssgsea_on_zscore$Sample, 'z_score__')

# Merge with metadata
ssgsea_on_zscore_plot <- ssgsea_on_zscore %>% left_join(df_clinical, by = 'Sample')

# Export 
saveRDS(ssgsea_on_zscore_plot, '08_biomarker_calculation/out/02_DEA_zscore_ssGSEA_df_plot.rds')
