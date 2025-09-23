setwd('~/Documents/projects/project_DGD/immunotherapy_biomarkers/')

library(tidyverse)
library(glue)

df_clinical <- readRDS('01_metadata/out/00_metadata.rds')

df_dgd_list <- readRDS('06_PCA/out/df_dgd_list.rds')
df_tpm_list <- readRDS('06_PCA/out/df_tpm_list.rds')

source('07_prep_biomarker_calculation/script/functions.R')

###############################################################
#################### 01 EDA - Density plots ################### 
###############################################################

# TPM density plots 
lapply(df_tpm_list, plot_density_tpm)

# DGD normalization density plots 
mapply(plot_density_dgd, df = df_dgd_list, study = names(df_dgd_list))

###############################################################
####################### Prep for ssGSEA ####################### 
###############################################################

mat_tpm_list <- list()
mat_dgd_list <- list()

for (study in names(df_tpm_list)){

  # Convert to matrix with genes (rows) x samples (cols)
  mat_tpm <- column_to_rownames(df_tpm_list[[study]], 'genes') %>% dplyr::select(!study) %>% as.matrix()
  mat_tpm_list[[study]] <- mat_tpm

  mat_dgd <- df_dgd_list[[study]] %>% dplyr::select(!c(tissue, study, sample)) %>% as.matrix() %>% t()
  colnames(mat_dgd) <- df_dgd_list[[study]]$sample
  mat_dgd_list[[study]] <- mat_dgd

}

rm(study)

saveRDS(mat_tpm_list, '07_prep_biomarker_calculation/out/mat_dgd_list.rds')
saveRDS(mat_dgd_list, '07_prep_biomarker_calculation/out/mat_tpm_list.rds')




