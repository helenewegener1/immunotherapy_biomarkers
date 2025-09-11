setwd('~/Documents/projects/project_DGD/cancer_sample_normalization/')

library(tidyverse)
library(glue)
library(biomaRt)

##################################### TGEP ##################################### 
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

################################################################################



# Use the ensembl IDs in final geneset since those are the IDs of the data
geneset <- list(TGEP = TGEP_ensembl$ensembl_gene_id)

saveRDS(geneset, 'rds/geneset.rds')
