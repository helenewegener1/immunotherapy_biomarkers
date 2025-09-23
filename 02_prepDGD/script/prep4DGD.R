setwd('~/Documents/projects/project_DGD/immunotherapy_biomarkers/')

library(dplyr)
library(glue)

cohorts <- readRDS('00_data/cohorts_all_allNorm.rds')

###############################################################
#################### Download all studies ##################### 
###############################################################

cancer_types <- names(cohorts)

for (cancer in cancer_types){
  
  studies <- cohorts[[cancer]] %>% names()
  
  for (study in studies){
    
    if ('raw' %in% names(cohorts[[cancer]][[study]])){

      counts <- cohorts[[cancer]][[study]]$raw
      counts_prep_dgd <- counts %>% t() %>% as.data.frame() %>% rownames_to_column('sample')
      counts_prep_dgd$tissue <- cancer # or other relevant though I don't think that it is used...

      name <- glue('{study}_counts')
      
      write_csv(counts_prep_dgd, glue('03_prepDGD/out/{name}_dgd_input.csv'))

    } else {
      print(glue('{cancer} {study} has no raw counts :( '))
    }
    
  }
  
}

files <- list.files('./03_prepDGD/out/')
csv_files <- files[endsWith(files, '.csv')]
datasets <- csv_files %>% str_remove('_counts_dgd_input.csv')
