setwd('~/Documents/projects/project_DGD/cancer_sample_normalization/')

library(tidyverse)
library(glue)
library(biomaRt)

# Load data
df_clinical <- readRDS('rds/00_metadata.rds')
geneset <- readRDS('rds/geneset.rds')

# Studies
studies <- df_clinical$study %>% unique()

# Initialize empty tibble
df_ES <- tibble(
  study = character(),
  sample = character(),
  TGEP = numeric()
)

for (this_study in studies){
  # this_study <- "Kim_GC"

  # Samples in study
  study_samples <- df_clinical %>% filter(study == this_study) %>% dplyr::select(Sample) %>% unlist()

  for (study_sample in study_samples){

    # Load DEA file from 05_DEA.py output
    dea_file <- suppressMessages(
      read_csv(glue('out/2025_09_09_DEA_{study_sample}.csv'), show_col_types = FALSE) %>% rename('genes' = '...1')
    )

    # Enrichment like in article
    # ES = (enriched cancer marker genes * 16,883 genes) / (significant genes * cancer marker genes)
    # A gene is DE if adjusted p value (q value) is < 0.01 and log2FC > 1 

    # enriched_cancer_marker_genes is the number of genes in the gene set that are DE
    enriched_cancer_marker_genes <- dea_file %>% filter(genes %in% geneset$TGEP & q_value < 0.01 & log2_fold_change > 1) %>% nrow()

    # significant_genes is the number of genes that are DE
    significant_genes <- dea_file %>% filter(q_value < 0.01 & log2_fold_change > 1) %>% nrow()

    # cancer_marker_genes is the number of genes in the gene list
    cancer_marker_genes <- length(geneset$TGEP)

    ES_TGEP <- (enriched_cancer_marker_genes * 16883) / (significant_genes * cancer_marker_genes)

    df_ES <- df_ES %>%
      add_row(study = this_study,
              sample = study_sample,
              TGEP = ES_TGEP)

  }

}

# Add meta data
df_ES_meta <- df_ES %>% left_join(df_clinical, by = c('sample' = 'Sample', 'study'))

# Export 
saveRDS(df_ES_meta, 'rds/04_03_ES_df_plot.rds')

########################## Different p value cutoffs ########################### 

# Initialize empty tibble
df_ES <- tibble(
  study = character(),
  sample = character(),
  TGEP = numeric(),
  cut_off = numeric()
)

for (this_study in studies){
  # this_study <- "Kim_GC"
  
  # Samples in study
  study_samples <- df_clinical %>% filter(study == this_study) %>% dplyr::select(Sample) %>% unlist()
  
  for (study_sample in study_samples){
    
    # Load DEA file from 05_DEA.py output
    dea_file <- suppressMessages(
      read_csv(glue('out/2025_09_09_DEA_{study_sample}.csv'), show_col_types = FALSE) %>% rename('genes' = '...1')
    )
    
    for (cut_off in seq(0.01, 0.1, by = 0.01)){
      
      # Enrichment like in article
      # ES = (enriched cancer marker genes * 16,883 genes) / (significant genes * cancer marker genes)
      # A gene is DE if adjusted p value (q value) is < 0.01 and log2FC > 1 
      
      # enriched_cancer_marker_genes is the number of genes in the gene set that are DE
      enriched_cancer_marker_genes <- dea_file %>% filter(genes %in% geneset$TGEP & q_value < cut_off & log2_fold_change > 1) %>% nrow()
      
      # significant_genes is the number of genes that are DE
      significant_genes <- dea_file %>% filter(q_value < cut_off & log2_fold_change > 1) %>% nrow()
      
      # cancer_marker_genes is the number of genes in the gene list
      cancer_marker_genes <- length(geneset$TGEP)
      
      ES_TGEP <- (enriched_cancer_marker_genes * 16883) / (significant_genes * cancer_marker_genes)
      
      df_ES <- df_ES %>%
        add_row(study = this_study,
                sample = study_sample,
                TGEP = ES_TGEP,
                cut_off = cut_off
                )
      
    }
    
  }
  
}

# Add meta data
df_ES_meta <- df_ES %>% left_join(df_clinical, by = c('sample' = 'Sample', 'study'))

cancer_types <- df_clinical$cancer %>% unique()
# this_cancer <- 'GC'

# Plot
for (this_cancer in cancer_types){
  
  df_ES_meta %>% 
    filter(cancer == this_cancer) %>% 
    # mutate(x_group = interaction(cut_off, Response, sep = "_")) %>% 
    ggplot(aes(x = Response, y = TGEP)) +
    facet_wrap(vars(as.character(cut_off))) + 
    geom_violin() +
    geom_boxplot(width = 0.05, outlier.shape = NA, color = 'grey') + 
    geom_jitter(width = 0.05) +
    labs(
      title = glue('{this_cancer}: Enrichment score the bulkDGD way'), 
      subtitle = 'ES = (enriched cancer marker genes * 16,883 genes) / (significant genes * cancer marker genes)',
      caption = 'ES of 1 means no enrichment',
      y = 'ES of TGEP genes'
    ) + 
    theme(
      plot.subtitle = element_text(size = 10) 
    ) +
    theme_bw() + 
    geom_hline(yintercept = 1, color = "darkred", size = 0.5) 
  
  ggsave(glue('plot/04_03_ES_cutoff_test_TGEP_{this_cancer}.png'), width = 13, height = 9)
  
}

  
  
  
  
  
  
  
  
  
