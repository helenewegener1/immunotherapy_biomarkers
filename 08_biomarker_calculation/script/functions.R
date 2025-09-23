
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

# biomarker_analysis <- function(study_list, score_list) {
#   
#   # study <- 'Riaz_MEL'
#   # score <- 'TGEP'
#   
#   for (study in study_list){
#     
#     # Load data from study
#     mat_tpm_raw <- mat_tpm_list[[study]]
#     mat_dgd_raw <- mat_dgd_list[[study]]
#     
#     # Subset intersection of genes in mat_tpm and mat_dgd in same study for comparison 
#     overlapping_genes <- intersect(rownames(mat_tpm_raw), rownames(mat_dgd_raw))
#     
#     mat_tpm <- mat_tpm_raw[overlapping_genes, ]
#     mat_dgd <- mat_dgd_raw[overlapping_genes, ]
#     
#     for (score in score_list){
#       
#       genes <- response_scores[[score]]
#       
#       # Map the ensembl genes in the data to gene symbols
#       mapping_ensembl_to_symbol <- getBM(
#         attributes = c("ensembl_gene_id", "hgnc_symbol"),
#         filters = "ensembl_gene_id",
#         values = overlapping_genes,
#         mart = mart
#       )
#       
#       # Find the ensembl mapping of the gene list 
#       geneset_mapping <- mapping_ensembl_to_symbol %>% filter(hgnc_symbol %in% genes)
#       
#       # Use the ensembl IDs in final geneset since those are the IDs of the data
#       geneset <- list(score = geneset_mapping$ensembl_gene_id)
#       
#       # Run ssGSEA and save output
#       ssgsea_tpm <- run_ssgsea(mat_tpm, geneset)  
#       ssgsea_dgd <- run_ssgsea(mat_dgd, geneset)
#       
#       # Prep for comparing output
#       df_ssgsea_tpm <- ssgsea_tpm %>% t() %>% as.data.frame() %>% 
#         rownames_to_column('Sample') %>% mutate(Data = 'TPM', )
#       
#       df_ssgsea_dgd <- ssgsea_dgd %>% t() %>% as.data.frame() %>% 
#         rownames_to_column('Sample') %>% mutate(Data = 'DGD')
#       
#       df_plot <- rbind(df_ssgsea_tpm, df_ssgsea_dgd)
#       
#       # Plot 
#       if (score == 'housekeeping'){
#         
#         df_plot %>% 
#           # ggplot(aes(y = score, fill = Data, x = interaction(!!sym(var_split), Data))) + 
#           ggplot(aes(y = score, fill = Data, x = Data)) + 
#           geom_boxplot() + 
#           # geom_jitter(aes(y = score, x = interaction(!!sym(var_split), Data)), alpha = 0.7, shape = 21) +
#           geom_jitter(aes(y = score, x = Data), alpha = 0.7, shape = 21) +
#           scale_fill_manual(values = c('grey30', 'darkgrey')) +
#           theme_bw() + 
#           theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#           labs(title = glue('{study} - {score} score'),
#                subtitle = glue('{score} is based on {length(genes)} genes.'),
#                y = glue('{score} score'))
#         
#         ggsave(glue('plot/04_02_{study}_{score}.png'), width = 20, height = 10)
#         
#       } else {
#         
#         for (var_split in c('Response', 'Biopsy_time', 'Sex', 'Treatment', 'cancer', 'reference_genome')){
#           
#           # df_plot %>% ggplot(aes(x = Sample, y = score, fill = Data)) + 
#           #   facet_wrap(vars(!!sym(var_split)), scales = "free_x", ncol = 1, nrow = unique(df_plot[[var_split]]) %>% length()) +
#           #   geom_col(position = 'dodge') + 
#           #   scale_fill_manual(values = c('grey30', 'darkgrey')) +
#           #   theme_bw() + 
#           #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#           #   labs(title = glue('{study} - {score} score'),
#           #        subtitle = glue('{score} is based on {length(genes)} genes.'),
#           #        y = glue('{score} score'))
#           # 
#           # ggsave(glue('plot/04_02_{study}_{score}_{var_split}.png'), width = 20, height = 10)
#           
#           df_plot %>% 
#             # filter(Data %in% c('TPM', 'DGD')) %>%
#             ggplot(aes(y = score, fill = Data, x = interaction(!!sym(var_split), Data))) + 
#             geom_boxplot() + 
#             geom_jitter(aes(y = score, x = interaction(!!sym(var_split), Data)), alpha = 0.7, shape = 21) + 
#             # scale_fill_manual(values = c('grey30', 'darkgrey')) +
#             theme_bw() + 
#             theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#             labs(title = glue('{study} - {score} score'),
#                  subtitle = glue('{score} is based on {length(genes)} genes.'),
#                  y = glue('{score} ssGSEA score'), 
#                  x = glue('{var_split} and data type'))
#           
#         ggsave(glue('plot/04_03_{study}_{score}_{var_split}.png'), width = 10, height = 8)
#           
#         }
#         
#       }
#       
#       # # Table output
#       # new_name_TPM <- glue("TPM_{score}")
#       # new_name_DGD <- glue("DGD_{score}")
#       # 
#       # df_table <- df_plot %>% 
#       #   pivot_wider(names_from = Data,
#       #               values_from = score) %>% 
#       #   rename(!!new_name_TPM := TPM,
#       #          !!new_name_DGD := DGD)
#       # 
#       # if(is.null(biomarker_result[[study]])){
#       #   df_table_export <- df_table
#       # } else {
#       #   df_table %>% dplyr::select(Sample, !!sym(new_name_TPM), !!sym(new_name_DGD))
#       #   df_table_export <- left_join(biomarker_result[[study]], df_table)
#       # }
#       # 
#       # biomarker_result[[study]] <- df_table_export
#       
#     }
#     
#   }
#   
#   return(biomarker_result)
#   
# }
# 
# # Cancer type wise biomarker analysis 
# ct_wise_biomarker_analysis <- function(cancer_type, cancer_suffix, score){ 
# 
#   # ct = cancer type
#   ct_studies <- names(mat_dgd_list)[grep(cancer_suffix, names(mat_dgd_list))]
#   
#   # Initialize mergning of expressions  
#   # ct_mat_dgd <- mat_dgd_list[[ct_studies[1]]] %>% as.data.frame() %>% rownames_to_column('genes')
#   ct_mat_tpm <- mat_tpm_list[[ct_studies[1]]] %>% as.data.frame() %>% rownames_to_column('genes')
#   
#   # Add rest of studies 
#   for (ct_study in ct_studies[2:length(ct_studies)]){
#     
#     # tmp_ct_mat_dgd <- mat_dgd_list[[ct_study]] %>% as.data.frame() %>% rownames_to_column('genes')
#     # ct_mat_dgd <- left_join(ct_mat_dgd, tmp_ct_mat_dgd, by = 'genes')
#     
#     tmp_ct_mat_tpm <- mat_tpm_list[[ct_study]] %>% as.data.frame() %>% rownames_to_column('genes')
#     ct_mat_tpm <- left_join(ct_mat_tpm, tmp_ct_mat_tpm, by = 'genes')
#     
#   }
#   
#   # Clean up 
#   rm(tmp_ct_mat_dgd, tmp_ct_mat_tpm)
#   
#   # Clinical data 
#   ct_df_clinical <- df_clinical %>% filter(cancer == cancer_type)
#   rownames(ct_df_clinical) <- ct_df_clinical$Sample
#   
#   # ct_mat_dgd <- ct_mat_dgd %>% column_to_rownames('genes') %>% as.matrix() %>% na.omit()
#   ct_mat_tpm <- ct_mat_tpm %>% column_to_rownames('genes') %>% as.matrix() %>% na.omit()
#   
#   # Subset intersection of genes in mat_tpm and mat_dgd in same study for comparison 
#   overlapping_genes <- intersect(rownames(ct_mat_dgd), rownames(ct_mat_tpm))
#   
#   ct_mat_tpm <- ct_mat_tpm[overlapping_genes, ]
#   # ct_mat_dgd <- ct_mat_dgd[overlapping_genes, ]
#   
#   genes <- response_scores[[score]]
#   
#   # Map the ensembl genes in the data to gene symbols
#   mapping_ensembl_to_symbol <- getBM(
#     attributes = c("ensembl_gene_id", "hgnc_symbol"),
#     filters = "ensembl_gene_id",
#     values = overlapping_genes,
#     mart = mart
#   )
#   
#   # Find the ensembl mapping of the gene list 
#   geneset_mapping <- mapping_ensembl_to_symbol %>% filter(hgnc_symbol %in% genes)
#   
#   # Use the ensembl IDs in final geneset since those are the IDs of the data
#   geneset <- list(score = geneset_mapping$ensembl_gene_id)
#   
#   # Run ssGSEA and save output
#   ssgsea_tpm <- run_ssgsea(ct_mat_tpm, geneset)  
#   ssgsea_dgd <- run_ssgsea(ct_mat_dgd, geneset)
#   
#   # Prep for comparing output
#   df_ssgsea_tpm <- ssgsea_tpm %>% t() %>% as.data.frame() %>% 
#     rownames_to_column('Sample') %>% mutate(Data = 'TPM', )
#   
#   # df_ssgsea_dgd <- ssgsea_dgd %>% t() %>% as.data.frame() %>% 
#   #   rownames_to_column('Sample') %>% mutate(Data = 'DGD')
#   
#   # Merge DGD and TPM run with metadata for plotting 
#   df_plot <- rbind(df_ssgsea_tpm, df_ssgsea_dgd)
#   df_plot_meta <- left_join(df_plot, ct_df_clinical, by = 'Sample')
#   
#   for (var_split in c('Response', 'Biopsy_time', 'Sex', 'Treatment', 'cancer', 'reference_genome')){
#     # var_split <- 'Response'
#     df_plot_meta %>% 
#       # filter(!is.na(!!sym(var_split))) %>%
#       ggplot(aes(y = score, fill = Data, x = interaction(!!sym(var_split), Data))) + 
#       geom_boxplot() + 
#       geom_jitter(aes(y = score, x = interaction(!!sym(var_split), Data)), alpha = 0.7, shape = 21) + 
#       # scale_fill_manual(values = c('grey30', 'darkgrey')) +
#       theme_bw() + 
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#       labs(title = glue('{cancer_type} - {score} score'),
#            subtitle = glue('{score} is based on {length(genes)} genes.'),
#            y = glue('{score} ssGSEA score'), 
#            x = glue('{var_split} and data type'))
#     
#     ggsave(glue('plot/04_04_{cancer_type}_{score}_{var_split}.png'), width = 10, height = 8)
#     
#   }
#   
#   # Using signifinder
#   log_mat_tpm <- log(ct_mat_tpm + 1)
#   log_mat_dgd <- log(ct_mat_dgd + 1)
#   
#   tpm_TinflamSign <- TinflamSign(log_mat_tpm, nametype = 'ENSEMBL')
#   dgd_TinflamSign <- TinflamSign(log_mat_dgd, nametype = 'ENSEMBL')
#   
#   # tpm_TinflamSign <- TinflamSign(mat_tpm, nametype = 'ENSEMBL')
#   # dgd_TinflamSign <- TinflamSign(mat_dgd, nametype = 'ENSEMBL')
#   
#   df_signi <- data.frame(Sample = colnames(ct_mat_tpm),
#                          'TPM_TinflamSign' = tpm_TinflamSign$Tinflam_Ayers,
#                          'DGD_TinflamSign' = dgd_TinflamSign$Tinflam_Ayers
#   )
#   
#   df_signi_long <- df_signi %>% pivot_longer(cols = !Sample,
#                                              names_to = 'Data',
#                                              values_to = 'score')
#   # Merge
#   df_plot_signi <- rbind(df_plot, df_signi_long)
#   
#   df_plot_signi_meta <- left_join(df_plot_signi, df_clinical, by = 'Sample')
#   
#   df_plot_signi_meta %>% 
#     ggplot(aes(y = score, fill = Data, x = Data)) + 
#     geom_boxplot() + 
#     geom_jitter(aes(y = score, x = Data), alpha = 0.7, shape = 21) + 
#     # scale_fill_manual(values = c('grey30', 'darkgrey')) +
#     theme_bw() + 
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(title = glue('{cancer_type} - {score} score'),
#          subtitle = glue('{score} is based on {length(genes)} genes.'),
#          y = glue('{score} ssGSEA score'), 
#          x = glue('{var_split} and data type'))
#   
#   ggsave(glue('plot/04_04_{cancer_type}_signifinder_{score}.png'), width = 10, height = 8)
#   
#   var_split <- 'Response'
#   
#   df_plot_signi_meta %>% 
#     filter(!is.na(!!sym(var_split))) %>%
#     ggplot(aes(y = score, fill = Data, x = interaction(!!sym(var_split), Data))) + 
#     geom_boxplot() + 
#     geom_jitter(aes(y = score, x = interaction(!!sym(var_split), Data)), alpha = 0.7, shape = 21) + 
#     # scale_fill_manual(values = c('grey30', 'darkgrey')) +
#     theme_bw() + 
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(title = glue('{cancer_type} - {score} score'),
#          subtitle = glue('{score} is based on {length(genes)} genes.'),
#          y = glue('{score} ssGSEA score'), 
#          x = glue('{var_split} and data type'))
#   
#   ggsave(glue('plot/04_04_{cancer_type}_signifinder_{score}_{var_split}.png'), width = 10, height = 8)
#   
# }
# 


