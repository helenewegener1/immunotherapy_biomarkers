setwd('~/Documents/projects/project_DGD/2025_08_28/')

library(tidyverse)
library(glue)
library(biomaRt)
library(GSVA)

source('../2025_08_27/biomarker_analysis/functions.R')
cohorts <- readRDS('data/wetransfer_bulkrnaseq_2025-08-28_1101/cohorts_all_allNorm.rds')

###############################################################
################## 01 Look at bulkDGD output ################## 
###############################################################

date <- '2025_08_28'
cancer <- 'Melanoma'
study <- 'Hugo_MEL'
# study <- 'Riaz_MEL'

# Load results 
preproc <- read_csv(glue('out/{date}_{study}_preprocessed.csv'))
repres <- read_csv(glue('out/{date}_{study}_representations.csv'))
# r_values <- read_csv(glue('out/{date}_{study}_pred_r_values.csv'))
deco_out <- read_csv(glue('out/{date}_{study}_decoder_outputs.csv'))
prob_rep <- read_csv(glue('out/{date}_{study}_prob_rep.csv'))
prob_comp <- read_csv(glue('out/{date}_{study}_prob_comp.csv'))

# How many genes preserved 
genes <- colnames(preproc)[colnames(preproc) %>% startsWith('ENSG')]
length(genes)

# Decoder outputs corresponding to the representations found. 
# The decoder outputs are the rescaled means of the negative binomial distributions 
# used to model the RNA-seq counts for the genes included in the DGD model.
deco_out
# deco_out_clean <- deco_out[, colnames(deco_out) %>% startsWith('ENSG')]

###############################################################
############################ 02 TPM ########################### 
###############################################################

tpm <- cohorts[[cancer]][[study]]$tpm %>% rownames_to_column('genes')

###############################################################
############################ 03 EDA ########################### 
###############################################################

preproc <- data_wrang(preproc)
deco_out <- data_wrang(deco_out)

# write_csv(preproc, glue('out/{date}_{srp_sample}_03_preprocessed.csv'))
# write_csv(deco_out, glue('out/{date}_{srp_sample}_03_decoder_outputs.csv'))

# Quick look 
dim(preproc)
dim(tpm) 
dim(deco_out)

# plot_dist(preproc, value_name = 'Counts')
plot_dist(tpm, value_name = 'TPM')
plot_dist(deco_out, value_name = 'DGD')

###############################################################
########################## 04 ssGSEA ########################## 
###############################################################

# Convert to matrix
mat_tpm <- column_to_rownames(tpm, 'genes') %>% as.matrix()

mat_preproc <- as.matrix(sapply(preproc, as.numeric))
mat_preproc <- mat_preproc[, colnames(mat_preproc) != 'genes']
rownames(mat_preproc) <- preproc$genes

mat_deco_out <- as.matrix(sapply(deco_out, as.numeric))
mat_deco_out <- mat_deco_out[, colnames(mat_deco_out) != 'genes']
rownames(mat_deco_out) <- deco_out$genes

# Translate gene symbols to ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# The 18 genes for T cell–inflamed GEP (gene expression profile)
genes <- c("TIGIT", "CD27", "CD8A", "PDCD1LG2", "LAG3", "CD274",
           "CXCR6", "CMKLR1", "NKG7", "CCL5", "PSMB10", "IDO1", 
           "CXCL9", "HLA-DQA1", "CD276", "STAT1", "HLA-DRB1", "HLA-E")

mapping <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = genes,
  mart = mart
)

# Keep the first Ensembl ID for each gene
mapping_unique <- mapping[!duplicated(mapping$hgnc_symbol), ]

# THE FINAL TGEP GENESET
TGEP_geneset <- list(TGEP = mapping$ensembl_gene_id)
TGEP_geneset_long <- list(TGEP_long = mapping_unique$ensembl_gene_id)

##################################################### 
################# SHORT GENESET LIST ################
##################################################### 
# Run ssGSEA and save output
# ssgsea_preproc <- run_ssgsea(mat_preproc, TGEP_geneset) # Does not make sense to run on raw counts
ssgsea_tpm <- run_ssgsea(mat_tpm, TGEP_geneset)  
ssgsea_deco_out <- run_ssgsea(mat_deco_out, TGEP_geneset)

# Prep for comparing output
tmp_tpm <- ssgsea_tpm %>% t() %>% as.data.frame() %>% 
  rownames_to_column('Sample') %>% mutate(Data = 'TPM', )

tmp_deco <- ssgsea_deco_out %>% t() %>% as.data.frame() %>% 
  rownames_to_column('Sample') %>% mutate(Data = 'DGD')

res <- rbind(tmp_tpm, tmp_deco)

# Plot 
res %>% ggplot(aes(x = Sample, y = TGEP, fill = Data)) + 
  geom_col(position = 'dodge') + 
  scale_fill_manual(values = c('darkgrey', 'black')) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = glue('{study}: 18 gene symbolms --> 18 ensembl genes'))

ggsave(glue('plot/{date}_{study}_TGEP.png'), width = 20, height = 5)

##################################################### 
################# LONG GENESET LIST ################# 
##################################################### 
ssgsea_long_tpm <- run_ssgsea(mat_tpm, TGEP_geneset_long)  
ssgsea_long_deco_out <- run_ssgsea(mat_deco_out, TGEP_geneset_long)

# Prep for comparing output
tmp_tpm_long <- ssgsea_long_tpm %>% t() %>% as.data.frame() %>% 
  rownames_to_column('Sample') %>% mutate(Data = 'TPM', )

tmp_deco_long <- ssgsea_long_deco_out %>% t() %>% as.data.frame() %>% 
  rownames_to_column('Sample') %>% mutate(Data = 'DGD')

res <- rbind(tmp_tpm_long, tmp_deco_long)

# Plot 
res %>% ggplot(aes(x = Sample, y = TGEP_long, fill = Data)) + 
  geom_col(position = 'dodge') + 
  scale_fill_manual(values = c('darkgrey', 'black')) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = glue('{study}: 18 gene symbolms --> 35 ensembl genes'))

ggsave(glue('plot/{date}_{study}_TGEP_long.png'), width = 20, height = 5)

#####################################################  
################ ChatGPT WAY OF GTEP ################ 
#####################################################  

# Subset TPM
gep_mat <- tpm %>%
  filter(genes %in% unlist(TGEP_geneset)) %>%
  column_to_rownames("genes")

# log2(TPM+1)
gep_log <- log2(gep_mat + 1)

# z-score per gene
gep_z <- t(scale(t(gep_log)))

# unweighted GEP score = mean of z-scores per sample
gep_score <- colMeans(gep_z, na.rm = TRUE)

gep_mat <- gep_score %>% as.data.frame() %>% rownames_to_column('Sample')
colnames(gep_mat) <- c('Sample', 'TGEP')
gep_mat$Data <- 'tpm_other_calc'

res <- rbind(tmp_tpm, tmp_deco, gep_mat)

res %>% ggplot(aes(x = Sample, y = TGEP, fill = Data)) + 
  geom_col(position = 'dodge') + 
  scale_fill_manual(values = c('darkgrey', 'black', 'green')) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = '18 gene symbolms --> 18 ensembl genes')


