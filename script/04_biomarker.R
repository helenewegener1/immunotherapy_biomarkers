setwd('~/Documents/projects/project_DGD/cancer_sample_normalization/')

library(tidyverse)
library(glue)
library(biomaRt)
library(GSVA)
library(signifinder)
# Sys.setenv(JAVA_HOME = "/Library/Java/JavaVirtualMachines/temurin-17.jdk/Contents/Home")
library(rJava)
library(xlsx)

source('script/functions.R')
df_clinical <- readRDS('rds/00_metadata.rds')

# df_dgd_list <- readRDS('rds/03_df_dgd_list.rds')
# df_tpm_list <- readRDS('rds/03_df_tpm_list.rds')

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
# 
# mat_tpm_list <- list()
# mat_dgd_list <- list()
# 
# for (study in names(df_tpm_list)){
# 
#   # Convert to matrix with genes (rows) x samples (cols)
#   mat_tpm <- column_to_rownames(df_tpm_list[[study]], 'genes') %>% dplyr::select(!study) %>% as.matrix()
#   mat_tpm_list[[study]] <- mat_tpm
# 
#   mat_dgd <- df_dgd_list[[study]] %>% dplyr::select(!c(tissue, study, sample)) %>% as.matrix() %>% t()
#   colnames(mat_dgd) <- df_dgd_list[[study]]$sample
#   mat_dgd_list[[study]] <- mat_dgd
# 
# }
# 
# rm(study)
# 
# saveRDS(mat_tpm_list, 'rds/04_mat_dgd_list.rds')
# saveRDS(mat_dgd_list, 'rds/04_mat_tpm_list.rds')


###############################################################
########################## 02 ssGSEA ########################## 
###############################################################

source('script/functions.R')

# Load matrix format of DGD and TPM data
mat_dgd_list <- readRDS('rds/04_mat_dgd_list.rds')
mat_tpm_list <- readRDS('rds/04_mat_tpm_list.rds')

# Load metadata
# cohorts <- readRDS('data/wetransfer_bulkrnaseq_2025-08-28_1101/cohorts_all_allNorm.rds')

# Translate gene symbols to ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# List of gene sets 
response_scores <- list(housekeeping = c("ACTB", "GAPDH", "B2M", "PPIA", "RPLP0", "HPRT1", "PGK1"), 
                        TGEP = c("TIGIT", "CD27", "CD8A", "PDCD1LG2", "LAG3", "CD274",
                                 "CXCR6", "CMKLR1", "NKG7", "CCL5", "PSMB10", "IDO1", 
                                 "CXCL9", "HLA-DQA1", "CD276", "STAT1", "HLA-DRB1", "HLA-E"))

# ONLY FOR INITIALIZING
# biomarker_result <- list()
# biomarker_result <- readRDS('rds/04_biomarker_result.rds')

study_list_full <- names(mat_dgd_list)
score_list_full <- names(response_scores)

# One score, one study
# biomarker_result <- biomarker_analysis(study_list = c('Riaz_MEL'), score_list = c('TGEP'))

# Add a new score to all data sets - remember to add the geneset of the score in response_scores
# biomarker_result <- biomarker_analysis(study_list = study_list_full, score_list = c('NEW_SCORE'))

# Add all scores to new data set - remember to preprocess the new data and add it to mat_tpm_list and mat_dgd_list
# biomarker_result <- biomarker_analysis(study_list = c('NEW_STUDY'), score_list = score_list_full)

# Run with all datasets and all scores
biomarker_analysis(study_list = study_list_full, score_list = score_list_full)

# Save in excel 
# openxlsx::write.xlsx(biomarker_result, file = "table/04_biomarker_scores.xlsx")
# saveRDS(biomarker_result, 'rds/04_biomarker_result.rds')

###############################################################
######################### Signifinder ######################### 
###############################################################

library(signifinder)

availableSignatures(description = FALSE) %>% view()
availableSignatures() 

# sysdata.rda downlaoded from signifinder github: https://github.com/CaluraLab/signifinder/blob/main/R/sysdata.rda
# load('~/Downloads/sysdata.rda')

Tinflam_Ayers %>% filter(class == 'TInflam') %>% dplyr::select(SYMBOL) %>% unlist() %in% response_scores$TGEP

ExpandedImmune_Ayers %>% dplyr::select(SYMBOL) %>% unlist() %in% response_scores$TGEP

IFN_Ayers %>% dplyr::select(SYMBOL) %>% unlist() %in% response_scores$TGEP

###############################################################
##################### Merge cancer types ###################### 
###############################################################

# Load matrix format of DGD and TPM data
mat_dgd_list <- readRDS('rds/04_mat_dgd_list.rds')
mat_tpm_list <- readRDS('rds/04_mat_tpm_list.rds')

# Biomarker for cancer type 

cancer_type <- "Melanoma"
cancer_suffix <- "_MEL"
score <- 'TGEP'
# score <- 'housekeeping'

ct_wise_biomarker_analysis(cancer_type = "Melanoma", cancer_suffix = "_MEL", score = "TGEP")






