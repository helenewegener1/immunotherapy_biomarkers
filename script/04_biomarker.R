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
cohorts <- readRDS('data/wetransfer_bulkrnaseq_2025-08-28_1101/cohorts_all_allNorm.rds')

# Translate gene symbols to ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# List of gene sets 
response_scores <- list(housekeeping = c("ACTB", "GAPDH", "B2M", "PPIA", "RPLP0", "HPRT1", "PGK1"), 
                        TGEP = c("TIGIT", "CD27", "CD8A", "PDCD1LG2", "LAG3", "CD274",
                                 "CXCR6", "CMKLR1", "NKG7", "CCL5", "PSMB10", "IDO1", 
                                 "CXCL9", "HLA-DQA1", "CD276", "STAT1", "HLA-DRB1", "HLA-E"))

# ONLY FOR INITIALIZING
biomarker_result <- list()
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
biomarker_result <- biomarker_analysis(study_list = study_list_full, score_list = score_list_full)

# Save in excel 
openxlsx::write.xlsx(biomarker_result, file = "table/04_biomarker_scores.xlsx")
saveRDS(biomarker_result, 'rds/04_biomarker_result.rds')

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

MEL_studies <- names(mat_dgd_list)[grep('_MEL', names(mat_dgd_list))]

MEL_mat_dgd <- mat_dgd_list[[MEL_studies[1]]] %>% as.data.frame() %>% rownames_to_column('genes')
MEL_mat_tpm <- mat_tpm_list[[MEL_studies[1]]] %>% as.data.frame() %>% rownames_to_column('genes')
MEL_df_clinical <- data.frame()
cancer_type <- "Melanoma"

for (MEL_study in MEL_studies[2:length(MEL_studies)]){
  
  tmp <- mat_dgd_list[[MEL_study]] %>% as.data.frame() %>% rownames_to_column('genes')
  MEL_mat_dgd <- left_join(MEL_mat_dgd, tmp, by = 'genes')
  
  tmp <- mat_tpm_list[[MEL_study]] %>% as.data.frame() %>% rownames_to_column('genes')
  MEL_mat_tpm <- left_join(MEL_mat_tpm, tmp, by = 'genes')
  
  tmp <- cohorts[[cancer_type]][[MEL_study]]$Clinical %>% 
    dplyr::select(Sample, Response, Biopsy_time) %>% 
    mutate(Sample = glue('{MEL_study}_{Sample}'))
  
  MEL_df_clinical <- rbind(MEL_df_clinical, tmp)
  
}

rm(tmp)
MEL_mat_dgd <- MEL_mat_dgd %>% column_to_rownames('genes') %>% as.matrix() %>% na.omit()
MEL_mat_tpm <- MEL_mat_tpm %>% column_to_rownames('genes') %>% as.matrix() %>% na.omit()
rownames(MEL_df_clinical) <- MEL_df_clinical$Sample

# Define variables (unction parameters)
mat_tpm <- MEL_mat_tpm
mat_dgd <- MEL_mat_dgd
clinical <- MEL_df_clinical

for (score in names(response_scores))  {
  
  # score <- 'TGEP'
  # score <- 'housekeeping'
  
  genes <- response_scores[[score]]
  
  # Map the ensembl genes in the data to gene symbols
  mapping_ensembl_to_symbol <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = rownames(mat_tpm),
    mart = mart
  )
  
  # Find the ensembl mapping of the gene list 
  geneset_mapping <- mapping_ensembl_to_symbol %>% filter(hgnc_symbol %in% genes)
  
  # Use the ensembl IDs in final geneset since those are the IDs of the data
  geneset <- list(score = geneset_mapping$ensembl_gene_id)
  
  # Run ssGSEA and save output
  ssgsea_tpm <- run_ssgsea(mat_tpm, geneset)  
  ssgsea_dgd <- run_ssgsea(mat_dgd, geneset)
  
  # Prep for comparing output
  df_ssgsea_tpm <- ssgsea_tpm %>% t() %>% as.data.frame() %>% 
    rownames_to_column('Sample') %>% mutate(Data = 'TPM', )
  
  df_ssgsea_dgd <- ssgsea_dgd %>% t() %>% as.data.frame() %>% 
    rownames_to_column('Sample') %>% mutate(Data = 'DGD')
  
  df_plot <- rbind(df_ssgsea_tpm, df_ssgsea_dgd)
  df_plot_meta <- left_join(df_plot, MEL_df_clinical, by = 'Sample')
  
  for (var_split in c('Response', 'Biopsy_time')){
    # var_split <- 'Response'
    df_plot_meta %>% 
      # filter(!is.na(!!sym(var_split))) %>%
      ggplot(aes(y = score, fill = Data, x = interaction(!!sym(var_split), Data))) + 
      geom_boxplot() + 
      geom_jitter(aes(y = score, x = interaction(!!sym(var_split), Data)), alpha = 0.7, shape = 21) + 
      scale_fill_manual(values = c('grey30', 'darkgrey')) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = glue('{cancer_type} - {score} score'),
           subtitle = glue('{score} is based on {length(genes)} genes.'),
           y = glue('{score} ssGSEA score'), 
           x = glue('{var_split} and data type'))
    
    ggsave(glue('plot/04_04_{cancer_type}_{score}_{var_split}.png'), width = 10, height = 8)
    
  }
  
  # Using signifinder
  log_mat_tpm <- log(mat_tpm + 1)
  log_mat_dgd <- log(mat_dgd + 1)
  
  tpm_TinflamSign <- TinflamSign(log_mat_tpm, nametype = 'ENSEMBL')
  dgd_TinflamSign <- TinflamSign(log_mat_dgd, nametype = 'ENSEMBL')
  
  # tpm_TinflamSign <- TinflamSign(mat_tpm, nametype = 'ENSEMBL')
  # dgd_TinflamSign <- TinflamSign(mat_dgd, nametype = 'ENSEMBL')
  
  df_signi <- data.frame(Sample = colnames(mat_tpm),
                         'TPM_TinflamSign' = tpm_TinflamSign$Tinflam_Ayers,
                         'DGD_TinflamSign' = dgd_TinflamSign$Tinflam_Ayers
                        )
  
  df_signi_long <- df_signi %>% pivot_longer(cols = !Sample,
                                             names_to = 'Data',
                                             values_to = 'score')
  # Merge
  df_plot_signi <- rbind(df_plot, df_signi_long)
  
  df_plot_signi_meta <- left_join(df_plot_signi, MEL_df_clinical, by = 'Sample')
  
  var_split <- 'Response'
  df_plot_signi_meta %>% 
    filter(!is.na(!!sym(var_split))) %>%
    ggplot(aes(y = score, fill = Data, x = interaction(!!sym(var_split), Data))) + 
    geom_boxplot() + 
    geom_jitter(aes(y = score, x = interaction(!!sym(var_split), Data)), alpha = 0.7, shape = 21) + 
    # scale_fill_manual(values = c('grey30', 'darkgrey')) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = glue('{cancer_type} - {score} score'),
         subtitle = glue('{score} is based on {length(genes)} genes.'),
         y = glue('{score} ssGSEA score'), 
         x = glue('{var_split} and data type'))
  
  ggsave(glue('plot/04_04_{cancer_type}_signifinder_{score}_{var_split}.png'), width = 10, height = 8)

}






