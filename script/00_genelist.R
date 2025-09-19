setwd('~/Documents/projects/project_DGD/cancer_sample_normalization/')

library(tidyverse)
library(glue)
library(biomaRt)
library(signifinder)
library(dplyr)

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

######################## Find gene sets on signifinder ######################### 

availableSignatures(description = FALSE) %>% view()
availableSignatures() 

# sysdata.rda downlaoded from signifinder github: https://github.com/CaluraLab/signifinder/blob/main/R/sysdata.rda
load('~/Downloads/sysdata.rda')

# Not sure how to filter these
EMT_Mak %>% head()
EMT_Mak$class %>% table()

EMT_Thompson %>% head()
EMT_Thompson$class %>% table()

PassON_Du %>% head()
PassON_Du$class %>% table()

Hypoxia_Buffa_genes <- Hypoxia_Buffa %>% dplyr::select(SYMBOL) %>% unlist()
ImmuneCyt_Rooney_genes <- ImmuneCyt_Rooney %>% dplyr::select(SYMBOL) %>% unlist()
Tinflam_Ayers_genes <- Tinflam_Ayers %>% filter(class == 'TInflam') %>% dplyr::select(SYMBOL) %>% unlist() 
Tinflam_Thompson_genes <- Tinflam_Thompson %>% dplyr::select(SYMBOL) %>% unlist()
TLS_Cabrita_genes <- TLS_Cabrita %>% dplyr::select(SYMBOL) %>% unlist()
CIN_Carter_25_genes <- CIN_Carter %>% filter(class == 'CIN25') %>% dplyr::select(SYMBOL) %>% unlist()
CIN_Carter_70_genes <- CIN_Carter %>% filter(class == 'CIN70') %>% dplyr::select(SYMBOL) %>% unlist()

ECM_Chakravarthy_up_genes <- ECM_Chakravarthy %>% filter(class == 'ECMup') %>% dplyr::select(SYMBOL) %>% unlist()
ECM_Chakravarthy_down_genes <- ECM_Chakravarthy %>% filter(class == 'ECMdown') %>% dplyr::select(SYMBOL) %>% unlist()
VEGF_Hu_genes <- VEGF_Hu %>% dplyr::select(SYMBOL) %>% unlist()

################################################################################

################# Translate gene lists from SYMBOL to ensembl ################## 

gene_lists_symbol <- list(
  Hypoxia_Buffa_genes, ImmuneCyt_Rooney_genes, Tinflam_Ayers_genes,
  Tinflam_Thompson_genes, TLS_Cabrita_genes, CIN_Carter_25_genes,
  CIN_Carter_70_genes, ECM_Chakravarthy_up_genes, ECM_Chakravarthy_down_genes,
  VEGF_Hu_genes
)

names(gene_lists_symbol) <- c(
  'Hypoxia_Buffa', 'ImmuneCyt_Rooney', 'Tinflam_Ayers',
  'Tinflam_Thompson', 'TLS_Cabrita', 'CIN_Carter_25',
  'CIN_Carter_70', 'ECM_Chakravarthy_up', 'ECM_Chakravarthy_down',
  'VEGF_Hu'
)

# Find the ensembl mapping of the gene list 
gene_lists_symbol_ensembl <- lapply(gene_lists_symbol, function(x) mapping_ensembl_to_symbol %>% filter(hgnc_symbol %in% x))

# Use the ensembl IDs in final geneset since those are the IDs of the data

geneset <- lapply(gene_lists_symbol_ensembl, function(x) x[['ensembl_gene_id']])

################################################################################

# Export genesets
saveRDS(geneset, 'rds/geneset.rds')
