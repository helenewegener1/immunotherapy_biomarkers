setwd('~/Documents/projects/project_DGD/immunotherapy_biomarkers/')

library(tidyverse)
library(glue)
library(biomaRt)
library(GSVA)

cohorts <- readRDS('00_data/cohorts_all_allNorm.rds')

cohorts$Melanoma$Hugo_MEL$Clinical %>% head()

# Number of data sets per cancer type 
lapply(cohorts, function(x) length(x))

# Reference genome
ref_genomes <- data.frame(study = c("Auslander_MEL", "Cho_NSCLC", "Choueiri_aRCC",
                                    "Clouphesy_GBM", "Du_MEL", "Gide_MEL", "Hugo_MEL", 
                                    "Jung_NSCLC", "Kim_GC", "Liu_MEL", "Mamdani_rEAC", 
                                    "Mariathasan_UC", "McDermott_RCC", "Ravi_NSCLC", 
                                    "Riaz_MEL", "VanAllen_aMEL", "Vandenende_rEAC", 
                                    "Zappasodi_MEL", "Zhao_GBM"),
                          reference_genome = c("GrCh38 - 91", "GrCh38 - 91", "GrCh38 - 103",
                                               "GrCh38 - 103", "GrCh38 - 91", "GrCh38 - 103", "GrCh38 - 103", 
                                               "GrCh38 - 91", "GrCh38 - 103", "GrCh38 - 103", "GrCh38 - 91", 
                                               "GrCh38 - 103", "GrCh38 - 103", "GRCh37",
                                               "GrCh38 - 103", "GrCh38 - 103", "GrCh38 - 107/90", 
                                               "GrCh38 - 107/90", NA))

study_list <- ref_genomes$study

df_clinical <- data.frame()

for (study in study_list){
  
  which_cancer_type <- names(cohorts)[sapply(cohorts, function(x) study %in% names(x))]
  
  df_clinical_tmp <- cohorts[[which_cancer_type]][[study]]$Clinical %>% 
    dplyr::select(Sample, Response, Biopsy_time, Treatment) %>% 
    mutate(Sample = glue('{study}_{Sample}'))
  
  # Sex if present 
  if ('Gender' %in% colnames(cohorts[[which_cancer_type]][[study]]$Clinical)){
    df_clinical_tmp$Sex <- cohorts[[which_cancer_type]][[study]]$Clinical[['Gender']]
  } else if ('Sex' %in% colnames(cohorts[[which_cancer_type]][[study]]$Clinical)){
    df_clinical_tmp$Sex <- cohorts[[which_cancer_type]][[study]]$Clinical[['Sex']]
  } else {
    df_clinical_tmp$Sex <- NA
  }
  
  # Age if present 
  if ('Age' %in% colnames(cohorts[[which_cancer_type]][[study]]$Clinical)){
    df_clinical_tmp$Age <- cohorts[[which_cancer_type]][[study]]$Clinical[['Age']]
  } else {
    df_clinical_tmp$Age <- NA
  }
  
  # RECIST if present 
  if ('RECIST' %in% colnames(cohorts[[which_cancer_type]][[study]]$Clinical)){
    df_clinical_tmp$RECIST <- cohorts[[which_cancer_type]][[study]]$Clinical[['RECIST']]
  } else {
    df_clinical_tmp$RECIST <- NA
  }
  
  df_clinical_tmp$study <- study
  df_clinical_tmp$cancer <- which_cancer_type
  
  df_clinical <- rbind(df_clinical, df_clinical_tmp)
  
}

# Add reference genome 
df_clinical <- df_clinical %>% left_join(ref_genomes, by = 'study')

# Export 
saveRDS(df_clinical, '01_metadata/out/00_metadata.rds')


