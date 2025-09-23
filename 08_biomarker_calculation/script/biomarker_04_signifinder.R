setwd('~/Documents/projects/project_DGD/immunotherapy_biomarkers/')

library(tidyverse)
library(glue)
library(signifinder)

# Load data
df_clinical <- readRDS('01_metadata/out/00_metadata.rds')
mat_tpm_list <- readRDS('07_prep_biomarker_calculation/out/mat_tpm_list.rds')
geneset <- readRDS('05_genesets/out/geneset.rds')
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
  
}

# Add metadata 
df_plot <- left_join(df_signifinder, df_clinical, by = 'Sample')

# Export 
saveRDS(df_plot, '08_biomarker_calculation/out/04_signifinder_df_plot.rds')


