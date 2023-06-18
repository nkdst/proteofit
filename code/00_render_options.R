# libraries ----
library(dplyr)



# parameters ----

## path variabls ----
project.wd <- rprojroot::find_root(".Rhistory")
path.analysis <- file.path(".", "code", "analysis")
path.reports <- file.path(".", "reports")

## document parameters ----
pCutoff <- 0.01
FCutoff <- 1

# colors:
condition_pal <- c("WT" = "turquoise3", "KO" = "indianred1")
tissue_pal <- c("gastroc" = "orange", "soleus" = "purple")
regulated_pal <- list(upregulated = 'red', downregulated = 'royalblue', insignificant = 'gray')
heatmap_pal <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
treemap_pal <- colorRampPalette(RColorBrewer::brewer.pal(11, "Paired"))

# parameter vector containing all general parameters:
all_parameters <-
  list(
    pCutoff = pCutoff,
    FCutoff = FCutoff,
    tissue_pal = tissue_pal,
    condition_pal = condition_pal,
    regulated_pal = regulated_pal,
    heatmap_pal = heatmap_pal,
    treemap_pal = treemap_pal
  )


# render functions ----

# ' this function looks at the yaml of the target file and extracts then only those
# ' given parameters which are needed for this document.
getFittingParametersForDocument <- function(filepath, parameters) {
  document.params <- rmarkdown::yaml_front_matter(filepath)$params
  
  needed_params <- names(parameters) %in% names(document.params) %>%
    parameters[.]
  return(needed_params)
}

# ' rendering the given Rmd file to a `html_document`
customMarkdownRendering <- function(filename,
                                    filename.addition = "",
                                    parameters = all_parameters,
                                    path.input = path.analysis,
                                    path.output = path.reports,
                                    path.root = project.wd) {
  filepath <-
    file.path(path.input, filename) # no file-extension yet
  filepath.rmd <- paste0(filepath, ".Rmd")
  rmd.params <-
    getFittingParametersForDocument(filepath = filepath.rmd, parameters)
  
  
  rmarkdown::render(
    input = paste0(filepath, ".Rmd"),
    output_file = paste0(filepath, filename.addition, ".html"),
    output_dir = path.output,
    output_format = "html_document",
    params = rmd.params,
    knit_root_dir = path.root
  )
}