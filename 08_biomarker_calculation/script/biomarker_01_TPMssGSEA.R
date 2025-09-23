setwd('~/Documents/projects/project_DGD/immunotherapy_biomarkers/')

library(tidyverse)
library(glue)
library(biomaRt)
library(GSVA)

source('08_biomarker_calculation/script/functions.R')

# Load data
mat_tpm_list <- readRDS('07_prep_biomarker_calculation/out/mat_tpm_list.rds')
df_clinical <- readRDS('01_metadata/out/00_metadata.rds')
geneset <- readRDS('05_genesets/out/geneset.rds')

# Prep output
df_mat_tpm <- NULL

for (mat_tpm in mat_tpm_list){
  
  tmp_mat_tpm <- mat_tpm %>% as.data.frame() %>% rownames_to_column('genes')

  # Merge dea files from all samples
  if (is.null(df_mat_tpm)){
    df_mat_tpm <- tmp_mat_tpm
  } else {
    df_mat_tpm <- left_join(df_mat_tpm, tmp_mat_tpm, by = 'genes')
  }
  
}

# Prep for ssGSEA
df_mat_tpm <- df_mat_tpm %>% column_to_rownames('genes') %>% as.matrix()

# Run ssGSEA and save output
ssgsea_on_tpm <- run_ssgsea(df_mat_tpm, geneset) 

# Reformat ssGSEA output
ssgsea_on_tpm <- ssgsea_on_tpm %>% t() %>% as.data.frame() %>% rownames_to_column('Sample')

# Merge with metadata
ssgsea_on_tpm_plot <- ssgsea_on_tpm %>% left_join(df_clinical, by = 'Sample')

# Export 
saveRDS(ssgsea_on_tpm_plot, '08_biomarker_calculation/out/01_TPMssGSEA_df_plot.rds')
