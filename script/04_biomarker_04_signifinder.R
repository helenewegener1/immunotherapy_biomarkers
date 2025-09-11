setwd('~/Documents/projects/project_DGD/cancer_sample_normalization/')

library(tidyverse)
library(glue)

########################### First look at signifinder ########################## 
library(signifinder)

availableSignatures(description = FALSE) %>% view()
availableSignatures() 

# sysdata.rda downlaoded from signifinder github: https://github.com/CaluraLab/signifinder/blob/main/R/sysdata.rda
# load('~/Downloads/sysdata.rda')

Tinflam_Ayers #%>% filter(class == 'TInflam') %>% dplyr::select(SYMBOL) %>% unlist() %in% response_scores$TGEP

################################################################################    

# Load data
df_clinical <- readRDS('rds/00_metadata.rds')
mat_tpm_list <- readRDS('rds/04_mat_tpm_list.rds')

# Prep output
df_tpm_TinflamSign <- NULL

# signifinder
for (mat_tpm in mat_tpm_list){
  
  # Data is log transformed in function so no need to do it before running the function
  # TGEP is calculated with household gene normalization and specific weight on the 18 genes - like it was intended to
  tpm_TinflamSign <- TinflamSign(mat_tpm, nametype = 'ENSEMBL')

  df_tpm_TinflamSign_tpm <- data.frame(Sample = colnames(mat_tpm),
                                       'TGEP' = tpm_TinflamSign$Tinflam_Ayers
  )

  # Merge dea files from all samples
  if (is.null(df_tpm_TinflamSign)){
    df_tpm_TinflamSign <- df_tpm_TinflamSign_tpm
  } else {
    df_tpm_TinflamSign <- rbind(df_tpm_TinflamSign, df_tpm_TinflamSign_tpm)
  }
  
}

# Add metadata 
df_plot <- left_join(df_tpm_TinflamSign, df_clinical, by = 'Sample')

# Export 
saveRDS(df_plot, 'rds/04_04_signifinder_df_plot.rds')


