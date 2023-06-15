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

# parameter vector containing all general parameters:
all_parameters <-
  list(
    pCutoff = pCutoff,
    FCutoff = FCutoff,
    tissue_pal = tissue_pal,
    condition_pal = condition_pal,
    regulated_pal = regulated_pal,
    heatmap_pal = heatmap_pal
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


# Documents --------------------------------------------------------------------
# TODO: make workflow diagram!


## 01 First Analysis and  pre-processing ----
customMarkdownRendering(filename = "01_FirstAnalysis")

## 03 DESeq ----
customMarkdownRendering(filename = "03_DESeq")


## 04 fgsea ----
customMarkdownRendering(
  filename = "04_fgsea",
  filename.addition = "_maxSize200_",
  parameters = c(all_parameters, reevaluate = F, fgsea.maxSize = 200)
)


## 05 fgseaRes pathways ----
customMarkdownRendering(filename = "05_fgseaRes_pathways",
                        parameters = c(all_parameters, reevaluate = FALSE))


## 06 Proteomics ----
customMarkdownRendering(filename = "06_Proteomics")


## 07 currentPlots ----
customMarkdownRendering(filename = "07_currentPlots")


## 08 GO (instead of fgsea) ----
customMarkdownRendering(filename = "08_GO")




## testing ----
# rmarkdown::render(
#   file.path("./code/testing.Rmd"),
#   output_file = "./code/testing.html",
#   output_dir = "./code",
#   output_format = "html_document",
#   # output_format = rmarkdown::html_notebook(toc = TRUE, toc_float = TRUE),
#   knit_root_dir = project.wd
# )