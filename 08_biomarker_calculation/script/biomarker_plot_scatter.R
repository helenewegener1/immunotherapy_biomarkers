setwd('~/Documents/projects/project_DGD/immunotherapy_biomarkers/')

library(tidyverse)
library(glue)
library(patchwork)

# Load data 
df_clinical <- readRDS('01_metadata/out/00_metadata.rds')
df_plot_TPMssGSEA <- readRDS('08_biomarker_calculation/out/01_TPMssGSEA_df_plot.rds') 
df_plot_DEA_logFC_ssGSEA <- readRDS('08_biomarker_calculation/out/02_DEA_logFC_ssGSEA_df_plot.rds') 
df_plot_DEA_zscore_ssGSEA <- readRDS('08_biomarker_calculation/out/02_DEA_zscore_ssGSEA_df_plot.rds') 
df_plot_ES <- readRDS('08_biomarker_calculation/out/03_ES_df_plot.rds') 
df_plot_signifinder <- readRDS('08_biomarker_calculation/out/04_signifinder_df_plot.rds') 
geneset <- readRDS('05_genesets/out/geneset.rds')

# Get iterables 
cancer_types <- df_clinical$cancer %>% unique()
studies <- df_clinical$study %>% unique()
meta_vars <- c('Response', 'Biopsy_time', 'Treatment', 'Sex', 'RECIST')

for (this_cancer in cancer_types){
  for (this_geneset in names(geneset)){
    
    df_tpm <- df_plot_signifinder %>% 
      filter(cancer == this_cancer) %>% 
      dplyr::select(Sample, !!sym(this_geneset)) %>% 
      rename(!!glue("{this_geneset}_signifinder") := !!sym(this_geneset))
    
    df_plot <- df_plot_DEA_logFC_ssGSEA %>% 
      filter(cancer == this_cancer) %>% 
      dplyr::select(Sample, !!sym(this_geneset), Response, Treatment) %>% 
      rename(!!glue("{this_geneset}_DEA_logFC_ssGSEA") := !!sym(this_geneset)) %>% 
      left_join(df_tpm, by = "Sample")
    
    ggplot(df_plot,
           aes(
             x = !!sym(glue("{this_geneset}_DEA_logFC_ssGSEA")), 
             y = !!sym(glue("{this_geneset}_signifinder")),
             color = Response,
             shape = Treatment
           ) ) + 
      geom_point() +
      theme_bw() + 
      # scale_color_manual(values = c("darkred", "forestgreen"))
      labs(
        title = glue("{this_cancer}: {this_geneset}")
      )
    
    ggsave(glue("08_biomarker_calculation/plot/scatter_{this_cancer}_{this_geneset}.pdf"), width = 8, height = 6)
    
  }
}

