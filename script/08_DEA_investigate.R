
# 
# # Load data
# dea_files <- readRDS('rds/07_dea_files.rds')
# df_clinical <- readRDS('rds/00_metadata.rds')
# 
# # Longer 
# dea_files_clean <- dea_files %>% 
#   pivot_longer(cols = -genes) %>% 
#   mutate(score = str_split_i(name, '__', 1), 
#          sample = str_split_i(name, '__', 2),
#          .keep = 'unused') %>% 
#   pivot_wider(names_from = score, 
#               values_from = value)
# 
# 
# # Merge with meta data
# dea_files_clean <- dea_files_clean %>% left_join(df_clinical, by = c('sample' = 'Sample'))
# 
# saveRDS(dea_files_clean, 'rds/08_dea_files_clean.rds')
dea_files_clean <- readRDS('rds/08_dea_files_clean.rds')

# From here - test differences in z scores across R/NR for each gene. 