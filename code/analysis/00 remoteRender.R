source(file.path("code", "analysis", "99_render_options.R"))

# Documents --------------------------------------------------------------------

## 01 First Analysis and  pre-processing ----
customMarkdownRendering(filename = "01_FirstAnalysis")

## 02 DEA ----
customMarkdownRendering(filename = "02_DEA")


## 03 pathway annotations ----
customMarkdownRendering(filename = "03_annotation")

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