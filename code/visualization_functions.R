library(dplyr)

if (!requireNamespace("stringr")) {
  install.packages("stringr") # used in pca function
}


library(ggplot2)
library(ggfortify) # autoplot (PCA)

if (! "ComplexHeatmap" %in% row.names(installed.packages())) {
  devtools::install_github("jokergoo/ComplexHeatmap", force = TRUE)
  # BiocManager::install("ComplexHeatmap")
}
library(ComplexHeatmap)



# ' plotting heatmap of a counts-table for both tissue types with arbitrary datatype 
# '
# ' the counts.df needs to be arranged, so that the gastrocnemius samples comprise
# ' the first 12 columns (6 WT, 6 KO)
plot_combined_CHM <- function(counts.df, datatype="zscore", params) {
  # http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3
  morecols <- colorRampPalette(params$heatmap_pal)
  
  # annotations:
  tissue_vec <- c(rep("gastroc", 12), rep("soleus", 12)) # needs to use the same title as the params$tissue_pal
  condition_vec <- c(rep(c(rep("WT", 6), rep("KO", 6)),2)) # needs to use the same title as the params$condition_pal
  
  top_annot <-
    HeatmapAnnotation(
      tissue = tissue_vec,
      condition = condition_vec,
      col = list(tissue = params$tissue_pal, condition = params$condition_pal),
      gp = gpar(col = "darkgray"),
      show_legend = FALSE
    )
  
  ## sample numbers
  sample_nrs <- colnames(counts.df) %>% 
    stringr::str_extract("(?<=_)\\d{1,2}(?=_)")
  colnames(counts.df) = sample_nrs
  
  chm <- Heatmap(
    as.matrix(counts.df),
    row_title = "most variable genes",
    name = datatype,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_title = "samples",
    col = morecols(50),
    column_title_side = "bottom",
    top_annotation = top_annot,
    column_names_gp = gpar(col = rep(unname(params$tissue_pal), each=12))

  )
  
  # creating custom annotation legend (to obtain the gray border)
  condition_legend = Legend(
    labels = c("WT", "KO"),
    legend_gp = gpar(fill = params$condition_pal),
    border = "darkgray",
    title = "condition"
  )
  tissue_legend = Legend(
    labels = c("gastrocnemius", "soleus"),
    legend_gp = gpar(fill = params$tissue_pal),
    border = "darkgray",
    title = "tissue"
  )
  legend_list <- list(condition_legend, tissue_legend)
  
  draw(chm, annotation_legend_list = legend_list)
}


# ' perform and plot a principle component analysis on a counts table 
plot_combined_pca <- function(counts.df, params){
    
  pca <- prcomp(t(counts.df), scale. = T)
  
  pca.data <- data.frame(Sample = rownames(pca$x),
                         X = pca$x[, 1],
                         Y = pca$x[, 2]) %>%
    mutate(condition = substr(Sample, 1, 2),
           tissue = stringr::str_extract(Sample, "[:alpha:]+$"),
           sample_nr = stringr::str_extract(Sample, "(?<=_)\\d{1,2}(?=_)"))
  
  
  plt <-
    autoplot(
      pca,
      data = pca.data,
      colour = 'tissue',
      shape = 'condition',
      label.label = "sample_nr",
      label = T,
      label.show.legend = FALSE, # (avoiding overlapping labels in legend)
      label.repel = T,
      
    ) +
    scale_color_manual(values = unname(params$tissue_pal)) + 
    theme(legend.text = element_text(size = 13),
          legend.title = element_text(size = 15)) +
    guides(shape = guide_legend(override.aes = list(size = 4)),
           color = guide_legend(override.aes = list(size = 4)))
  
  return(plt)
}