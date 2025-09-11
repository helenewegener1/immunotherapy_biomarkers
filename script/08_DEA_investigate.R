setwd('~/Documents/projects/project_DGD/cancer_sample_normalization/')

library(tidyverse)
library(glue)
library(patchwork)

df_clinical <- readRDS('rds/00_metadata.rds')

# # Load data
# dea_files <- readRDS('rds/04_02_dea_files.rds')
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

cancer_types <- df_clinical$cancer %>% unique()
meta_vars <- c('Response', 'Biopsy_time', 'Treatment', 'Sex', 'RECIST', 'study')

############################### PCA of z scores ################################ 

z_score_pca_prep <- dea_files_clean %>% 
  dplyr::select(genes, sample, z_score) %>% 
  pivot_wider(names_from = genes, 
              values_from = z_score) %>% 
  column_to_rownames('sample') %>% 
  as.matrix()

# Make sure there are no Inf nor -Inf in the matrix
z_score_pca_prep[is.infinite(z_score_pca_prep) & z_score_pca_prep > 0] <- 10   # replace Inf
z_score_pca_prep[is.infinite(z_score_pca_prep) & z_score_pca_prep < 0] <- -10  # replace -Inf

# Remove zero-variance columns
z_score_pca_prep <- z_score_pca_prep[, apply(z_score_pca_prep, 2, function(x) sd(x) > 0)]

# Run PCA
pca_z_score <- prcomp(z_score_pca_prep, scale. = TRUE)

# Save results in object that is ready for gg-plotting 
df_pca_z_score <- pca_z_score$x %>% 
  as.data.frame() %>% 
  dplyr::select(PC1, PC2, PC3) %>% 
  rownames_to_column('sample') %>% 
  left_join(df_clinical, by = c('sample' = 'Sample'))

# Variance explained per PC
pca_z_score_var_explained <- (pca_z_score$sdev^2) / sum(pca_z_score$sdev^2)

# Percentage explained for first 10 PCs
pca_z_score_var_explained_percent <- pca_z_score_var_explained[1:10] * 100

for (this_cancer in cancer_types){
  for (var in meta_vars){
    
    df_pca_z_score %>%
      filter(cancer == this_cancer) %>% 
      ggplot(aes(x = PC1,
                 y = PC3,
                 color = !!sym(var))) + 
      geom_point() + 
      theme_bw() + 
      labs(title = glue('{this_cancer}: Z score --> PCA'), 
           x = glue('PC1 ({round(pca_z_score_var_explained_percent[1], 1)}% variance explained)'),
           y = glue('PC2 ({round(pca_z_score_var_explained_percent[2], 1)}% variance explained)'))
    
    ggsave(glue('plot/08_PCA_zscore_{this_cancer}_{var}.png'), height = 8, width = 12)
    
  }
}


############################ PCA of TGEP z scores ############################## 

geneset <- readRDS('rds/geneset.rds')

TGEP_z_score_pca_prep <- dea_files_clean %>% 
  filter(genes %in% geneset$TGEP) %>% 
  dplyr::select(genes, sample, z_score) %>% 
  pivot_wider(names_from = genes, 
              values_from = z_score) %>% 
  column_to_rownames('sample') %>% 
  as.matrix()

# Remove zero-variance columns
# TGEP_z_score_pca_prep <- TGEP_z_score_pca_prep[, apply(TGEP_z_score_pca_prep, 2, function(x) sd(x) > 0)]

# Make sure there are no Inf nor -Inf in the matrix
TGEP_z_score_pca_prep[is.infinite(TGEP_z_score_pca_prep) & TGEP_z_score_pca_prep > 0] <- 10   # replace Inf
TGEP_z_score_pca_prep[is.infinite(TGEP_z_score_pca_prep) & TGEP_z_score_pca_prep < 0] <- -10  # replace -Inf

# Run PCA
pca_TGEP_z_score <- prcomp(TGEP_z_score_pca_prep, scale. = TRUE)

# Save results in object that is ready for gg-plotting 
df_pca_TGEP_z_score <- pca_TGEP_z_score$x %>% 
  as.data.frame() %>% 
  dplyr::select(PC1, PC2, PC3, PC4, PC5) %>% 
  rownames_to_column('sample') %>% 
  left_join(df_clinical, by = c('sample' = 'Sample'))

# Variance explained per PC
pca_TGEP_z_score_var_explained <- (pca_TGEP_z_score$sdev^2) / sum(pca_TGEP_z_score$sdev^2)

# Percentage explained for first 10 PCs
pca_TGEP_z_score_var_explained_percent <- pca_TGEP_z_score_var_explained[1:10] * 100

for (this_cancer in cancer_types){
  for (var in meta_vars){
    
    df_pca_TGEP_z_score %>%
      filter(cancer == this_cancer) %>% 
      ggplot(aes(x = PC1,
                 y = PC3,
                 color = !!sym(var))) + 
      geom_point() + 
      theme_bw() + 
      labs(title = glue('{this_cancer}: TGEP Z score --> PCA'), 
           x = glue('PC1 ({round(pca_TGEP_z_score_var_explained_percent[1], 1)}% variance explained)'),
           y = glue('PC2 ({round(pca_TGEP_z_score_var_explained_percent[2], 1)}% variance explained)'))
    
    ggsave(glue('plot/08_PCA_TGEP_zscore_{this_cancer}_{var}.png'), height = 8, width = 12)
    
  }
}

