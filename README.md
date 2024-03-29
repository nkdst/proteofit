- [About this project](#about-this-project)
	- [performed analysis steps](#performed-analysis-steps)
- [Using this repository](#using-this-repository)
- [Structure](#structure)
	- [Directory structure](#directory-structure)
	- [Code structure](#code-structure)


# About this project

This repository is part of my bachelor thesis, which was done in collaboration
with Helmholtz Zentrum München (Menden Lab) and the Bartels Lab at the Institute
for Cardiovascular Prevention (IPEK) part of LMU. The thesis is connected to the
Proteofit project, hence the preliminary name of this repository (it does not represent it!).

## performed analysis steps

- differential expression analysis (**DEA**) with `DESeq2`
- gene set enrichment analysis (**GSEA**) with `fgsea`
- **GO enrichment analysis** with `clusterProfiler`





# Using this repository

The RNA-Seq data (count tables) is not available in this repository and needs to be manually added before running the code.
The count tables are expected to be in the `/data/counts` directory and named `readcount_genename_gastroc.xls` and `readcount_genename_soleus.xls` respectively.



# Structure

## Directory structure


- `/code/analysis` : .Rmd files for the analysis steps
- `/data` : counts data, metadata, RData objects, etc.
- `/plots`: plots generated by the code.
- `/reports`: knitted reports of the .Rmd files


## Code structure

The code is split into several .Rmd files, which are numbered in the order they should be run. They can be dynamically knitted into html notebooks, which are saved in the same directory. All `.Rmd` files can similar also be knitted into html reports.

Best way to render/knit the files is to use the `00_remoteRender.R` script, which is also located in the `/code` directory. It provides the possibility to render all `.Rmd` files in the `/code/analysis` directory and save them as html documents in the `/reports` directory. The main reason for it to be used is that it allows for better control over the rendering process: I noticed that rendering the files directly from within RStudio uses graphic settings like resolution depending on which monitor is open or which device is used. 
This lead to inconsistent plots. To remedy this the `rmarkdown::render()` function is used to produce more consistent results.

Generally the pictures saved in the `/plots` directory are the ones used in the thesis, as they are optimized to have the right format and resolution. Even with rendering the plots with the dedicated script, some figures are still displayed suboptimal in the knitted files.

Additionally `00_remoteRender.R` contains parameters to for example decide the color of the condition, or adjust the p-value cutoff for all documents in one place.
