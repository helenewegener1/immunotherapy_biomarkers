setwd('~/Documents/projects/project_DGD/cancer_sample_normalization/')

library(tidyverse)
library(glue)
library(biomaRt)
library(GSVA)

source('script/functions.R')

############################################################
########################### TGEP ########################### 
############################################################

# Gene list 
TGEP <- c("TIGIT", "CD27", "CD8A", "PDCD1LG2", "LAG3", "CD274",
          "CXCR6", "CMKLR1", "NKG7", "CCL5", "PSMB10", "IDO1", 
          "CXCL9", "HLA-DQA1", "CD276", "STAT1", "HLA-DRB1", "HLA-E")

# Translate genes in dea_file from ensembl to SYMBOLS
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Load one DEA file to get genes 
tmp_dea_file <- read_csv('out/2025_09_09_DEA_Kim_GC_PB-16-045.csv') %>% rename('genes' = '...1')

# Map the ensembl genes in the data to gene symbols
mapping_ensembl_to_symbol <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = tmp_dea_file$genes,
  mart = mart
)

# Find the ensembl mapping of the gene list 
TGEP_ensembl <- mapping_ensembl_to_symbol %>% filter(hgnc_symbol %in% TGEP)

# Use the ensembl IDs in final geneset since those are the IDs of the data
geneset <- list(TGEP = TGEP_ensembl$ensembl_gene_id)

############################################################
############################################################
############################################################

# Load data
df_clinical <- readRDS('rds/00_metadata.rds')

# Studies 
studies <- df_clinical$study %>% unique()

# Initialize empty tibble
# df_DEA_ssGSEA <- tibble(
#   Sample = character(),
#   TGEP = character(),
#   study = numeric()
# )

dea_files <- NULL

# for (this_study in studies){
  # this_study <- "Kim_GC"
  
# All samples
samples <- df_clinical %>% dplyr::select(Sample) %>% unlist()

for (sample in samples){
  
  # sample <- samples[[2]]
  
  # Load DEA file from 05_DEA.py output 
  dea_file <- suppressMessages(
    read_csv(glue('out/2025_09_09_DEA_{sample}.csv'), show_col_types = FALSE) %>% rename('genes' = '...1')
  )
  
  dea_file <- dea_file %>%
    mutate(z_score = sign(log2_fold_change) * qnorm(1 - p_value/2))
  
  colnames(dea_file) <- c('genes', glue("p_value__{sample}"), glue("q_value__{sample}"), 
                          glue("log2_fold_change__{sample}"), glue("z_score__{sample}"))
  
  # dea_file_out <- dea_file %>% dplyr::select(genes, z_score) %>% rename(!!glue("z_score_{sample}") := z_score)
  
  # Merge dea files from all samples
  if (is.null(dea_files)){
    dea_files <- dea_file
  } else {
    dea_files <- dea_files %>% left_join(dea_file, by = 'genes')
  }
  
}

# saveRDS(dea_files, 'rds/07_dea_files.rds')
dea_files <- readRDS('rds/07_dea_files.rds')

# Contains score for all samples 
dea_files_mat <- dea_files %>% dplyr::select(genes, starts_with('z_score')) %>% column_to_rownames('genes') %>% as.matrix()

# Run ssGSEA and save output
ssgsea_on_zscore <- run_ssgsea(dea_files_mat, geneset) 

# Reformat ssGSEA output
ssgsea_on_zscore <- ssgsea_on_zscore %>% t() %>% as.data.frame() %>% rownames_to_column('Sample')
ssgsea_on_zscore$Sample <- str_remove(ssgsea_on_zscore$Sample, 'z_score__')

# # Merge with ES score
# df_ES <- readRDS('rds/06_df_ED.rds')
# ssgsea_on_zscore <- ssgsea_on_zscore %>% left_join(df_ES, by = c('Sample' = 'sample')) %>% dplyr::select(-study)
# ssgsea_on_zscore_long <- ssgsea_on_zscore %>% pivot_longer(cols = -Sample)

# Merge with metadata
df_clinical <- readRDS('rds/00_metadata.rds')
ssgsea_on_zscore_plot <- ssgsea_on_zscore %>% left_join(df_clinical, by = 'Sample')

# Plot 
var <- 'Response'
# this_cancer <- 'GC'
cancer_types <- df_clinical$cancer %>% unique()

for (this_cancer in cancer_types){
  
  var_table <- df_clinical %>%
    filter(cancer == this_cancer) %>%
    # filter(study == this_study) %>%
    group_by(!!sym(var)) %>%
    summarise(n = n())
  
  ssgsea_on_zscore_plot %>% 
    filter(cancer == this_cancer) %>%
    # filter(study == this_study) %>%
    ggplot(aes(y = TGEP, x = !!sym(var))) + 
    geom_violin() + 
    geom_jitter(width = 0.05) + 
    # geom_boxplot() + 
    geom_text(
      data = var_table,
      aes(x = !!sym(var), y = max(ssgsea_on_zscore_plot$TGEP), label = paste0("n = ", n)),
      inherit.aes = FALSE
    ) + 
    # facet_wrap(vars(name)) + 
    theme_bw() + 
    labs(title = glue("{this_cancer}: ssGSEA on z-score"),
         subtitle = 'z-score = sign(log2FC) * qnorm(1 - p_value/2)', 
         y = "ssGSEA score")
  
  ggsave(glue("plot/07_DEA_ssGSEA_{this_cancer}_{var}.png"), width = 10, height = 8)
  
  
}
