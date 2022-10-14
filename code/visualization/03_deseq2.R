## --------------------------- ##
##
## Purpose of script: DESeq2 Analysis
##
## Author: Nick Diercksen
##
## Date Created: 2022-09-23
##
## --------------------------- ##

# libraries ----
# BiocManager::install("genefilter")
# BiocManager::install("DESeq2")
library(DESeq2)
library(dplyr)
library(tibble)
library(readr)

# including other file(s) ----
# source("code/load_data.R") # is already loaded by 01_overview.R


# Gastroc tissue,for paired analysis use patient ----

## loading data ----
current_folder <- "data/x204sc20081494-z01-f004_mus_musculus_result/3.quantification/count"
count <- file.path(".", current_folder,"readcount_genename_gastroc.xls") %>%
  readr::read_tsv() %>%
  tibble::column_to_rownames(var = "gene_id")

metadata <- file.path(".", current_folder,"sample_table_gastroc.csv") %>%
  readr::read_csv2() %>%
  tibble::column_to_rownames("run")

# select sample columns
reorder_index <- match(rownames(metadata), colnames(count))
count <- count[,reorder_index]


# Check metadata consistency
all(rownames(metadata) %in% colnames(count)) %>%
  assertthat::assert_that(., msg = "metadata and count table do not match")


## DESeq2 object for  tissue, for paired analysis use patient
dds <- DESeqDataSetFromMatrix(countData =  count,
                              colData =  metadata,
                              design = ~ genotype)


### filter low expressed genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# main DESeq
ddsDE <- DESeq(dds)

# normalized counts
normCounts = counts(ddsDE, normalized = T)


# DESeq results
res = results(ddsDE, alpha = 0.05) # FDR cutoff

# summary(res)
resOrdered <- res[order(res$padj),]
plotMA(ddsDE)


# shrink (improves discoverability of genes with lower count rate)
res_shrink <- lfcShrink(ddsDE, coef = 2)
plotMA(res_shrink)





dds_sf <- estimateSizeFactors(dds)


# default analysis
dds_da <- DESeq(dds)
disp.plt <- plotDispEsts(dds_da, ylim = c(1e-6, 1e2) )
res <- results(dds_da)










# shrink
res_shrink <- lfcShrink(dds_da, coef = 2)



## preparing data ----
# making count table matrices
gastroc.cm <- counts.gastroc[, 1:13] %>%
  tibble::column_to_rownames(var = 'gene_id')
#   janitor::row_to_names(row_number = 1)
# class(gastroc.cm) <- "numeric"

gastroc.coldata <- counts.gastroc[, -c(2:13)] %>%
  tibble::column_to_rownames(var = 'gene_id')
