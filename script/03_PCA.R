setwd('~/Documents/projects/project_DGD/2025_08_28/')

library(tidyverse)
library(glue)
library(biomaRt)
library(GSVA)

cohorts <- readRDS('data/wetransfer_bulkrnaseq_2025-08-28_1101/cohorts_all_allNorm.rds')

date <- '2025_08_28'

# studies <- c('Hugo_MEL', 'Riaz_MEL')
# prefix <- 'two'

studies <- c("Auslander_MEL", "Braun_ccRCC", "Cho_NSCLC", "Choueiri_aRCC",
             "Clouphesy_GBM", "Du_MEL", "Gide_MEL", "Hugo_MEL", "Jung_NSCLC",
             "Kim_GC", "Liu_MEL", "Mamdani_rEAC", "Mariathasan_UC",
             "McDermott_RCC", "Miao_ccRCC", "Motzer_aRCC", "Nathanson_MEL",
             "Ravi_NSCLC", "Riaz_MEL", "Rose_mUC", "Snyder_UC",
             "VanAllen_aMEL", "Vandenende_rEAC", "Zappasodi_MEL", "Zhao_GBM")
prefix <- 'all'

###############################################################
##################### Load and merge data ##################### 
###############################################################

df_dgd_latent_space <- data_frame() # DGD Latent space
df_dgd_list <- list() # DGD decoder output (DGD normalization)

for (study in studies){
  
  # DGD Latent space
  tmp_df <- read_csv(glue('out/{date}_{study}_representations.csv'))
  tmp_df$study <- study
  tmp_df$sample <- glue('{study}_{tmp_df$sample}')
  
  df_dgd_latent_space <- rbind(df_dgd_latent_space, tmp_df)
  
  rm(tmp_df)

  # DGD decoder output (DGD normalization)
  tmp_df <- read_csv(glue('out/{date}_{study}_decoder_outputs.csv'))
  tmp_df$study <- study
  tmp_df$sample <- glue('{study}_{tmp_df$sample}')

  df_dgd_list[[study]] <- tmp_df

  rm(tmp_df)
  
}

# TPM

df_tpm_list <- list() # TPM normalization
cancer_types <- names(cohorts)

for (cancer in cancer_types){
  
  tmp_studies <- cohorts[[cancer]] %>% names()
  
  for (study in tmp_studies){
    
    tmp_df <- cohorts[[cancer]][[study]]$tpm %>% rownames_to_column('genes')
    tmp_df$study <- study
    colnames(tmp_df) <- c('genes', glue('{study}_{colnames(tmp_df)[2:(length(colnames(tmp_df))-1)]}'), 'study')
    
    df_tpm_list[[study]] <- tmp_df
    
    rm(tmp_df)
    
  }
  
  rm(tmp_studies)
  
}

# Subset
df_tpm_list <- df_tpm_list[studies]

# Export 
# saveRDS(df_dgd_list, 'rds/03_df_dgd_list.rds')
# saveRDS(df_tpm_list, 'rds/03_df_tpm_list.rds')
df_dgd_list <- readRDS('rds/03_df_dgd_list.rds')
df_tpm_list <- readRDS('rds/03_df_tpm_list.rds')


###############################################################
################## TPM: Gene subset for PCA ################### 
###############################################################

# Extract gene vectors for each dataset in list 
gene_lists <- lapply(df_tpm_list, function(x) x$genes)

# Find intersecting genes between all datasets 
common_genes <- Reduce(intersect, gene_lists)

# Subset each dataset to the intersecting genes
datasets_common <- lapply(df_tpm_list, function(x) {
  dplyr::filter(x, genes %in% common_genes)
})

# Ensure same order of genes across datasets
datasets_common <- lapply(datasets_common, function(x) {
  x[match(common_genes, x$genes), ]
})

# Merge datasets by "genes"
merged <- Reduce(function(x, y) full_join(x, y, by = "genes"), datasets_common)

# Save meta data (study information) separately ...
metadata <- bind_rows(lapply(names(datasets_common), function(study_name) {
  tibble(
    study = study_name,
    sample = colnames(datasets_common[[study_name]])[2:(ncol(datasets_common[[study_name]]) - 1)]
  )
}))

# ... and remove study column from expression 
expr <- merged[, !colnames(merged) %>% startsWith('study')]

# Pre for PCA
expr <- expr %>% column_to_rownames('genes')
expr_t <- t(expr) # Transpose to samples (rows) × genes (columns)

# PCA with scale. = TRUE tries to standardize each column (each gene) to unit variance. 
# If a gene has zero variance across samples, I get an error.
# So, I remove zero-variance columns
expr_t_nzv <- expr_t[, apply(expr_t, 2, function(x) sd(x) > 0)]
# ncol(expr_t) - ncol(expr_t_nzv) # N genes we loose 

# run PCA
pca_tpm <- prcomp(expr_t_nzv, scale. = TRUE)

# Save results in object that is ready for gg-plotting 
df_pca_tpm <- pca_tpm$x %>% as.data.frame() %>% dplyr::select(PC1, PC2, PC3) %>% rownames_to_column('sample')

# Merge with meta data bc we love colors <3 
df_pca_tpm <- df_pca_tpm %>% left_join(metadata, by = 'sample')

###############################################################
########### DGD normalization: Gene subset for PCA ############ 
###############################################################

# Extract gene vectors (ignore non-gene columns)
non_gene_cols <- c("sample", "tissue", "study")
gene_lists <- lapply(df_dgd_list, function(x) {
  setdiff(colnames(x), non_gene_cols)
})

# Find intersecting genes across all datasets
common_genes <- Reduce(intersect, gene_lists)

# Subset each dataset to only the common genes + sample/tissue/study
datasets_common <- lapply(df_dgd_list, function(x) {
  x %>%
    dplyr::select(all_of(c("sample", "tissue", "study", common_genes)))
})

# datasets_common <- lapply(datasets_common, function(x) {
#   x %>%
#     dplyr::mutate(
#       sample = as.character(sample),
#       tissue = as.character(tissue),
#       study  = as.character(study)
#     )
# })

# Merge all datasets by sample, tissue, study (keeping all genes)
merged <- dplyr::bind_rows(datasets_common)

# Save metadata separately
metadata <- merged %>% dplyr::select(sample, tissue, study)

# Prepare expression matrix for PCA (only numeric gene columns)
expr <- merged %>% dplyr::select(all_of(common_genes)) %>% as.data.frame()
rownames(expr) <- metadata$sample   # keep sample names as rownames

# Transpose to samples × genes
expr_t <- as.matrix(expr)   # samples are rows, genes are columns

# Remove zero-variance genes
# expr_t_nzv <- expr_t[, apply(expr_t, 2, sd) > 0]

#  Run PCA
pca_dgd <- prcomp(expr_t, scale. = TRUE)

# Prepare PCA results for plotting
df_pca_dgd <- as.data.frame(pca_dgd$x) %>%
  dplyr::select(PC1, PC2, PC3) %>%
  rownames_to_column("sample") %>%
  left_join(metadata, by = "sample")

###############################################################
##################### Plot TPM in PC space ####################
###############################################################

df_pca_tpm %>% 
  ggplot(aes(x = PC1, 
             y = PC2, 
             color = study)) + 
  geom_point() + 
  theme_bw() + 
  labs(title = 'TPM normalization --> PCA')

ggsave(glue('plot/03_{prefix}_PCA_TPM.png'), height = 8, width = 12)

###############################################################
############# Plot DGD normalization in PC space ##############
###############################################################

df_pca_dgd %>% 
  ggplot(aes(x = PC1, 
             y = PC2, 
             color = study)) + 
  geom_point() + 
  theme_bw() + 
  labs(title = 'DGD normalization --> PCA')

ggsave(glue('plot/03_{prefix}_PCA_DGD_normalization.png'), height = 8, width = 12)

###############################################################
####################### DGD Latent Space ######################
###############################################################

# See batch effect across dataframe
df_dgd_latent_space %>% 
  ggplot(aes(x = latent_dim_1, 
             y = latent_dim_2, 
             color = study)) + 
  geom_point() + 
  theme_bw() + 
  labs(title = 'DGD Latent variables')

ggsave(glue('plot/03_{prefix}_DGD_latent_variables.png'), height = 8, width = 12)
