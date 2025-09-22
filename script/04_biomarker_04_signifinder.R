setwd('~/Documents/projects/project_DGD/cancer_sample_normalization/')

library(tidyverse)
library(glue)
library(signifinder)

# Load data
df_clinical <- readRDS('rds/00_metadata.rds')
mat_tpm_list <- readRDS('rds/04_mat_tpm_list.rds')
geneset <- readRDS('rds/geneset.rds')
geneset_to_function <- availableSignatures(description = FALSE) 

# Get signifinder function names from geneset name
fun_names <- lapply(names(geneset), function(x) {
  
  # Normalize x if it matches certain patterns
  if (str_detect(x, "CIN_Carter")) {
    x <- "CIN_Carter"
  } else if (str_detect(x, "ECM_Chakravarthy")) {
    x <- "ECM_Chakravarthy"
  }
  
  # Filter geneset_to_function using (possibly updated) x
  geneset_to_function %>%
    filter(signature == x) %>%
    pull(functionName)   # Equivalent to select + unlist
  
})

names(fun_names) <- names(geneset)


# Prep output
df_signifinder <- NULL

# signifinder
for (geneset_name in names(geneset)){
  
  df_one_geneset <- NULL
  
  for (mat_tpm in mat_tpm_list){
    
    fun_name <- fun_names[[geneset_name]]
    
    # Run signifinder function 
    if (geneset_name == 'Tinflam_Ayers'){
      tmp <- do.call(fun_name, list(mat_tpm, nametype = "ENSEMBL", author = "Ayers"))
    } else if (geneset_name == 'Tinflam_Thompson'){
      tmp <- do.call(fun_name, list(mat_tpm, nametype = "ENSEMBL", author = "Thompson"))
    } else {
      tmp <- do.call(fun_name, list(mat_tpm, nametype = "ENSEMBL"))
    }
    
    # Transform to datafrane
    df_tmp <- data.frame(Sample = colnames(mat_tpm),
                         tmp_name = tmp[[1]]
    )
    
    # Rename colnames 
    colnames(df_tmp) <- c('Sample', geneset_name)
    
    # Merge dea files from all datasets 
    if (is.null(df_one_geneset)){
      df_one_geneset <- df_tmp
    } else {
      df_one_geneset <- rbind(df_one_geneset, df_tmp)
    }
    
  }
  
  if (is.null(df_signifinder)){
    df_signifinder <- df_one_geneset
  } else {
    df_signifinder <- left_join(df_signifinder, df_one_geneset, by = 'Sample')
  }
  
  
    
    # # Data is log transformed in function so no need to do it before running the function
    # # TGEP is calculated with household gene normalization and specific weight on the 18 genes - like it was intended to
    # tpm_TinflamSign <- TinflamSign(mat_tpm, nametype = 'ENSEMBL')
    # 
    # df_tpm_TinflamSign_tpm <- data.frame(Sample = colnames(mat_tpm),
    #                                      'TGEP' = tpm_TinflamSign$Tinflam_Ayers
    # )
    # 
    # # Merge dea files from all samples
    # if (is.null(df_tpm_TinflamSign)){
    #   df_tpm_TinflamSign <- df_tpm_TinflamSign_tpm
    # } else {
    #   df_tpm_TinflamSign <- rbind(df_tpm_TinflamSign, df_tpm_TinflamSign_tpm)
    # }
    
  
  
}

# Add metadata 
df_plot <- left_join(df_signifinder, df_clinical, by = 'Sample')

# Export 
saveRDS(df_plot, 'rds/04_04_signifinder_df_plot.rds')


