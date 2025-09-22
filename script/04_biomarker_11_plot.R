setwd('~/Documents/projects/project_DGD/cancer_sample_normalization/')

library(tidyverse)
library(glue)
library(patchwork)

# Load data 
df_clinical <- readRDS('rds/00_metadata.rds')
df_plot_TPMssGSEA <- readRDS('rds/04_01_TPMssGSEA_df_plot.rds') 
df_plot_DEAssGSEA <- readRDS('rds/04_02_DEAssGSEA_df_plot.rds') 
df_plot_ES <- readRDS('rds/04_03_ES_df_plot.rds') 
df_plot_signifinder <- readRDS('rds/04_04_signifinder_df_plot.rds') 
geneset <- readRDS('rds/geneset.rds')

# Get iterables 
cancer_types <- df_clinical$cancer %>% unique()
studies <- df_clinical$study %>% unique()
meta_vars <- c('Response', 'Biopsy_time', 'Treatment', 'Sex', 'RECIST')

# Plot one var
# this_cancer <- 'GC'
var <- "Response"
# var <- "Treatment"

for (this_cancer in cancer_types){
  
  # Count number in each group of "var"
  var_table <- df_clinical %>%
    filter(cancer == this_cancer) %>% 
    group_by(!!sym(var)) %>%
    summarise(n = n())
  
  for (this_geneset in names(geneset)){
    
    # this_geneset <- names(geneset)[[1]]
    
    # Initialize list with plots
    plot_list <- list(
      TPMssGSEA = df_plot_TPMssGSEA,
      DEAssGSEA = df_plot_DEAssGSEA,
      ES        = df_plot_ES,
      signifinder = df_plot_signifinder
    )
    
    # Plot for each dataset 
    plots <- lapply(names(plot_list), function(name) {
      
      df <- plot_list[[name]]
      
      df %>% 
        filter(cancer == this_cancer) %>% 
        dplyr::select(c(Sample, !!sym(this_geneset), !!sym(var))) %>% 
        ggplot(aes(y = !!sym(this_geneset), x = !!sym(var))) +
          geom_violin() +
          geom_boxplot(width = 0.15, outlier.shape = NA, color = 'grey') + 
          geom_jitter(width = 0.05) +
          geom_text(
            data = var_table,
            aes(x = !!sym(var),
                y = max(df[[this_geneset]][df[['cancer']] == this_cancer]) + 0.05,
                label = paste0("n = ", n)),
            inherit.aes = FALSE
          ) +
          ggtitle(glue('{name}, {this_geneset}, {this_cancer}')) +
          theme_bw()
    })
    
    # Names list of plots 
    names(plots) <- names(plot_list)
    
    # Add specific titles to each plot
    plots$TPMssGSEA <- plots$TPMssGSEA + labs(
      title = glue('{this_cancer}: ssGSEA on TPM'),
      y = glue('ssGSEA on {this_geneset} genes')
    )
    
    plots$DEAssGSEA <- plots$DEAssGSEA + labs(
      title = glue('{this_cancer}: ssGSEA on z-score from bulkDGD DEA'), 
      subtitle = 'z-score = sign(log2FC) * qnorm(1 - p_value/2)',
      y = glue('ssGSEA on {this_geneset} genes')
    )
    
    plots$ES <- plots$ES + labs(
      title = glue('{this_cancer}: Enrichment score the bulkDGD way'), 
      subtitle = 'ES = (enriched cancer marker genes * 16,883 genes) / (significant genes * cancer marker genes)',
      caption = 'ES of 1 means no enrichment',
      y = glue('ES of {this_geneset} genes')
    ) + 
      theme(
        plot.subtitle = element_text(size = 10) 
      ) +
      geom_hline(yintercept = 1, color = "darkred", size = 0.5)
    
    
    plots$signifinder <- plots$signifinder + labs(
      title = glue('{this_cancer}: TinflamSign signifinder = TGEP'), 
      subtitle = 'Normalized by household genes and specific weight on the 18 genes',
      y = glue('{this_geneset} score')
    ) 
    
    # Export plots in one figure
    wrap_plots(plots) + plot_annotation(
      title = glue("{this_cancer} cancer, {this_geneset} genes split by {var}")
      # subtitle = "Optional subtitle here",
      # caption = "Data source: mtcars"
    ) &
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(glue('plot/04_11_TGEP_{this_cancer}_{this_geneset}_{var}.png'), width = 13, height = 9)
    
    
  }
  
}


# Plot two vars
# this_cancer <- 'RCC'
var1 <- "Response"
var2 <- "Treatment"

for (this_cancer in cancer_types){
  
  n_var2 <- df_clinical %>%
    filter(cancer == this_cancer) %>% 
    distinct(!!sym(var2)) %>% 
    nrow()
  
  if (n_var2 > 1){
    
    # Count number in each group of "var"
    var_table <- df_clinical %>%
      filter(cancer == this_cancer) %>% 
      mutate(x_group = interaction(!!sym(var1), !!sym(var2), sep = "_")) %>% 
      group_by(x_group) %>%
      summarise(n = n())
    
    # Initialize list with plots
    plot_list <- list(
      TPMssGSEA = df_plot_TPMssGSEA,
      DEAssGSEA = df_plot_DEAssGSEA,
      ES        = df_plot_ES,
      signifinder = df_plot_signifinder
    )
    
    # Plot for each dataset 
    plots <- lapply(names(plot_list), function(name) {
      df <- plot_list[[name]]
      
      ggplot(
        df %>% filter(cancer == this_cancer) %>% mutate(x_group = interaction(!!sym(var1), !!sym(var2), sep = "_")),
        aes(y = TGEP, x = x_group)
      ) +
        geom_violin() +
        geom_boxplot(width = 0.15, outlier.shape = NA, color = 'grey') + 
        geom_jitter(width = 0.05) +
        geom_text(
          data = var_table,
          aes(x = x_group,
              y = max(df[['TGEP']][df[['cancer']] == this_cancer]) + 0.05,
              label = paste0("n = ", n)),
          inherit.aes = FALSE
        ) +
        ggtitle(name) +
        theme_bw() + 
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1)  # tilt 45 degrees
        )
    })
    
    # Names list of plots 
    names(plots) <- names(plot_list)
    
    # Add specific titles to each plot
    plots$TPMssGSEA <- plots$TPMssGSEA + labs(
      title = glue('{this_cancer}: ssGSEA on TPM'),
      y = 'ssGSEA on TGEP genes'
    )
    
    plots$DEAssGSEA <- plots$DEAssGSEA + labs(
      title = glue('{this_cancer}: ssGSEA on z-score from bulkDGD DEA'), 
      subtitle = 'z-score = sign(log2FC) * qnorm(1 - p_value/2)',
      y = 'ssGSEA on TGEG genes'
    )
    
    plots$ES <- plots$ES + labs(
      title = glue('{this_cancer}: Enrichment score the bulkDGD way'), 
      subtitle = 'ES = (enriched cancer marker genes * 16,883 genes) / (significant genes * cancer marker genes)',
      caption = 'ES of 1 means no enrichment',
      y = 'ES of TGEP genes'
    ) + 
      theme(
        plot.subtitle = element_text(size = 10) 
      ) +
      geom_hline(yintercept = 1, color = "darkred", size = 0.5)
    
    plots$signifinder <- plots$signifinder + labs(
      title = glue('{this_cancer}: TinflamSign signifinder = TGEP'), 
      subtitle = 'Normalized by household genes and specific weight on the 18 genes',
      y = 'TGEP score'
    ) 
    
    # Export plots in one figure
    wrap_plots(plots) 
    ggsave(glue('plot/04_11_TGEP_{this_cancer}_{var2}_{var1}.png'), width = 13, height = 9)
    
    
  }
  
 
}




