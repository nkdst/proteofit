source(file.path("code", "00_render_options.R"))
# use this ^ file to adjust parameters like cutoffs, color palettes…

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

## 05 GO (instead of fgsea) ----
customMarkdownRendering(filename = "05_GO")
