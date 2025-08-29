setwd('~/Documents/projects/project_DGD/2025_08_28/')

library(dplyr)
library(glue)

cohorts <- readRDS('data/wetransfer_bulkrnaseq_2025-08-28_1101/cohorts_all_allNorm.rds')

###############################################################
################# Save data for bulkDGD run ################### 
###############################################################

counts <- cohorts$Melanoma$Hugo_MEL$raw %>% t() %>% as.data.frame() %>% rownames_to_column('sample')
counts$tissue <- 'MEL'

name <- 'Hugo_MEL_counts'

write_csv(counts, glue('data/{name}_dgd_input.csv'))

###############################################################
###############################################################
###############################################################

counts <- cohorts$Melanoma$Riaz_MEL$raw %>% t() %>% as.data.frame() %>% rownames_to_column('sample')
counts$tissue <- 'MEL'

name <- 'Riaz_MEL_counts'

write_csv(counts, glue('data/{name}_dgd_input.csv'))

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
      
      write_csv(counts_prep_dgd, glue('data/{name}_dgd_input.csv'))

    } else {
      print(glue('{cancer} {study} has no raw counts :( '))
    }
    
  }
  
}

files <- list.files('./data/')
csv_files <- files[endsWith(files, '.csv')]
datasets <- csv_files %>% str_remove('_counts_dgd_input.csv')
