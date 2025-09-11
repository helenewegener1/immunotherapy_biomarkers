setwd('~/Documents/projects/project_DGD/cancer_sample_normalization/')

library(tidyverse)
library(glue)
library(biomaRt)

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

############################################################
############################################################
############################################################

# Load data
df_clinical <- readRDS('rds/00_metadata.rds')

# Studies 
studies <- df_clinical$study %>% unique()

# Initialize empty tibble
df_ES <- tibble(
  study = character(),
  sample = character(),
  ES = numeric()
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
    enriched_cancer_marker_genes <- dea_file %>% filter(genes %in% TGEP_ensembl$ensembl_gene_id & q_value < 0.01 & log2_fold_change > 1) %>% nrow()
    
    # significant_genes is the number of genes that are DE
    significant_genes <- dea_file %>% filter(q_value < 0.01 & log2_fold_change > 1) %>% nrow()
    
    # cancer_marker_genes is the number of genes in the gene list 
    cancer_marker_genes <- length(TGEP)
    
    ES <- (enriched_cancer_marker_genes * 16883) / (significant_genes * cancer_marker_genes)
    
    df_ES <- df_ES %>%
      add_row(study = this_study,
              sample = study_sample,
              ES = ES)
    
  }
  
}