project.wd <- rprojroot::find_root(".Rhistory")

# parameters
pCutoff <- 0.01
FCutoff <- 1


# path.analysis <- file.path("./code/analysis")
path.analysis <- file.path(".", "code", "analysis")
path.reports <- file.path(".", "reports")


# 03

rmd.params <- list(pCutoff = pCutoff, FCutoff = FCutoff)
filename.addition <- ""
rmarkdown::render(
  file.path(path.analysis, "03_DESeq.Rmd"),
  output_file = file.path(path.analysis, paste0("03_DESeq", filename.addition, ".html")),
  output_dir = path.reports,
  output_format = "html_document",
  params = rmd.params,
  knit_root_dir = project.wd
)

# 04

rmd.params <- list(pCutoff = pCutoff, reevaluate = F, fgsea.maxSize = 200)
filename.addition <- "_maxSize400_"
rmarkdown::render(
  file.path(path.analysis, "04_fgsea.Rmd"),
  output_file = file.path(path.analysis, paste0("04_fgsea", filename.addition, ".html")),
  output_dir = path.reports,
  output_format = "html_document",
  params = rmd.params,
  knit_root_dir = project.wd
)

# 05

rmd.params <- list(reevaluate = F)
filename.addition <- ""
rmarkdown::render(
  file.path(path.analysis, "05_fgseaRes_pathways.Rmd"),
  output_file = file.path(path.analysis, paste0("05_fgseaRes_pathways", filename.addition, ".html")),
  output_dir = path.reports,
  output_format = "html_document",
  params = rmd.params,
  knit_root_dir = project.wd
)



# 07

rmd.params <- list(pCutoff = pCutoff, FCutoff = FCutoff)
filename.addition <- ""
rmarkdown::render(
  file.path(path.analysis, "07_currentPlots.Rmd"),
  output_file = file.path(path.analysis, paste0("07_currentPlots", filename.addition, ".html")),
  output_dir = path.reports,
  output_format = "html_document",
  params = rmd.params,
  knit_root_dir = project.wd
)



# testing

rmarkdown::render(
  file.path("./code/testing.Rmd"),
  output_file = "./code/testing.html",
  output_dir = "./code",
  output_format = "html_document",
  # output_format = rmarkdown::html_notebook(toc = TRUE, toc_float = TRUE),
  knit_root_dir = project.wd
)
