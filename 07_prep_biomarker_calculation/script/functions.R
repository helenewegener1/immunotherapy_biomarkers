

plot_density_tpm <- function(df){
  
  long_format <- df %>%
    pivot_longer(
      cols = -c(genes, study),
      names_to = "Sample",
      values_to = "Normalization"
    )
  
  long_format$Normalization <- as.numeric(long_format$Normalization)
  
  ggplot(long_format, aes(x = log2(Normalization + 1), color = Sample)) +
    geom_density() +
    theme_minimal() +
    labs(x = 'TPM', title = glue("{unique(long_format$study)}: TPM Normalization Density per Sample")) +
    theme_bw() + 
    theme(legend.position = 'none')
  
  ggsave(glue('07_prep_biomarker_calculation/plot/density_{unique(long_format$study)}_TPM.png'), width = 15, height = 8)
  
}

plot_density_dgd <- function(df, study){
  
  df <- df %>% t() %>% as.data.frame()
  colnames(df) <- df['sample', ]
  df <- df[!rownames(df) == 'sample', ]
  df <- df %>% rownames_to_column('genes')
  
  long_format <- df %>%
    pivot_longer(
      cols = -c(genes),
      names_to = "Sample",
      values_to = "Normalization"
    )
  
  long_format$Normalization <- as.numeric(long_format$Normalization)
  
  ggplot(long_format, aes(x = log2(Normalization + 1), color = Sample)) +
    geom_density() +
    theme_minimal() +
    labs(x = 'DGD norm', title = glue("{study}: DGD Normalization Density per Sample")) +
    theme_bw() + 
    theme(legend.position = 'none')
  
  ggsave(glue('07_prep_biomarker_calculation/plot/density_{study}_DGD.png'), width = 15, height = 8)
  
}