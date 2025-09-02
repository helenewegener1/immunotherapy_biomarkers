

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
  
  ggsave(glue('plot/04_01_density_{unique(long_format$study)}_TPM.png'), width = 15, height = 8)
  
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
  
  ggsave(glue('plot/04_01_density_{study}_DGD.png'), width = 15, height = 8)
  
}

##################################################################
#################### Bio marker Calculations #####################
##################################################################

run_ssgsea <- function(mat_expr, geneset){
  
  # Expression matrix: genes x samples: (rownames = Ensembl IDs or gene symbols, colnames = samples)
  gsva_input <- ssgseaParam(exprData = mat_expr,
                            geneSets = geneset)
  
  # Run ssGSEA
  gep_scores <- gsva(gsva_input)
  
  return(gep_scores)
  
}

# Needs predefined:
# response_scores: named list of list of genesets in gene symbols
# run_ssgsea function found in ../2025_08_28/script/functions.R

biomarker_analysis <- function(study_list, score_list) {
  
  study <- 'Hugo_MEL'
  score <- 'TGEP'
  
  for (study in study_list){
    
    # Load data from study
    mat_tpm <- mat_tpm_list[[study]]
    mat_dgd <- mat_dgd_list[[study]]
    
    for (score in score_list){
      
      genes <- response_scores[[score]]
      
      # Map the ensembl genes in the data to gene symbols
      mapping_ensembl_to_symbol <- getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol"),
        filters = "ensembl_gene_id",
        values = rownames(mat_dgd),
        mart = mart
      )
      
      # Find the ensembl mapping of the gene list 
      geneset_mapping <- mapping_ensembl_to_symbol %>% filter(hgnc_symbol %in% genes)
      
      # Use the ensembl IDs in final geneset since those are the IDs of the data
      geneset <- list(score = geneset_mapping$ensembl_gene_id)
      
      # Run ssGSEA and save output
      ssgsea_tpm <- run_ssgsea(mat_tpm, geneset)  
      ssgsea_dgd <- run_ssgsea(mat_dgd, geneset)
      
      # Prep for comparing output
      df_ssgsea_tpm <- ssgsea_tpm %>% t() %>% as.data.frame() %>% 
        rownames_to_column('Sample') %>% mutate(Data = 'TPM', )
      
      df_ssgsea_dgd <- ssgsea_dgd %>% t() %>% as.data.frame() %>% 
        rownames_to_column('Sample') %>% mutate(Data = 'DGD')
      
      df_plot <- rbind(df_ssgsea_tpm, df_ssgsea_dgd)
      
      # Add metadata 
      ## Find the cancer_type that contains your study
      which_caner_type <- names(cohorts)[sapply(cohorts, function(x) study %in% names(x))]
      
      df_clinical <- cohorts[[which_caner_type]][[study]]$Clinical %>% 
        dplyr::select(Sample, Response, Biopsy_time) %>% 
        mutate(Sample = glue('{study}_{Sample}'))
      
      df_plot <- left_join(df_plot, df_clinical, by = 'Sample')
      
      # Plot 
      if (score == 'housekeeping'){
        
        df_plot %>% ggplot(aes(x = Sample, y = score, fill = Data)) +
          geom_col(position = 'dodge') + 
          scale_fill_manual(values = c('grey30', 'darkgrey')) +
          theme_bw() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
          labs(title = glue('{study} - {score} score'),
               subtitle = glue('{score} is based on {length(genes)} genes.'),
               y = glue('{score} score'))
        
        ggsave(glue('plot/04_02_{study}_{score}.png'), width = 20, height = 10)
        
      } else {
        
        for (var_split in c('Response', 'Biopsy_time')){
          
          df_plot %>% ggplot(aes(x = Sample, y = score, fill = Data)) + 
            facet_wrap(vars(!!sym(var_split)), scales = "free_x", ncol = 1, nrow = unique(df_plot[[var_split]]) %>% length()) +
            geom_col(position = 'dodge') + 
            scale_fill_manual(values = c('grey30', 'darkgrey')) +
            theme_bw() + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
            labs(title = glue('{study} - {score} score'),
                 subtitle = glue('{score} is based on {length(genes)} genes.'),
                 y = glue('{score} score'))
          
          ggsave(glue('plot/04_02_{study}_{score}_{var_split}.png'), width = 20, height = 10)
          
        }
        
      }
      
     
      # Table output
      new_name_TPM <- glue("TPM_{score}")
      new_name_DGD <- glue("DGD_{score}")
      
      df_table <- df_plot %>% 
        pivot_wider(names_from = Data,
                    values_from = score) %>% 
        rename(!!new_name_TPM := TPM,
               !!new_name_DGD := DGD)
      
      if(is.null(biomarker_result[[study]])){
        df_table_export <- df_table
      } else {
        df_table %>% dplyr::select(Sample, TPM_TGEP, DGD_TGEP)
        df_table_export <- left_join(biomarker_result[[study]], df_table)
      }
      
      biomarker_result[[study]] <- df_table_export
      
    }
    
  }
  
  return(biomarker_result)
  
}




