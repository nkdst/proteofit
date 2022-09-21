## --------------------------- ##
##
## Purpose of script: creating a pca plot(s)
##
## Author: Nick Diercksen
##
## Date Created: 2022-09-19
##
## --------------------------- ##

## libraries
library(dplyr)
library(janitor)
library(ggplot2)
library(ggfortify)
library(ggrepel)


### including other file(s)
# source("code/load_data.R") # is already loaded by 01_overview.R

# making count table matrices
counts.gastroc.m <- counts.gastroc[, 1:13] %>%
  t() %>%
  janitor::row_to_names(row_number = 1)
class(counts.gastroc.m) <- "numeric"

# normalized
counts.gastroc.norm <- sweep(counts.gastroc.m,1,total_readcounts.gastroc$reads,FUN="/")

pca <- prcomp(counts.gastroc.norm)

pca.data <- data.frame(Sample = rownames(pca$x),
                       X = pca$x[, 1],
                       Y = pca$x[, 2]) %>%
  mutate(type = substr(Sample, 1,2))

# TODO: how to get the label 'a' removed from the legend?
plt <- autoplot(pca, data = pca.data, colour = 'type', label.repel=TRUE) +
  ggtitle("PCA gastroc normalized")
ggsave(filename = "./plots/02_pca_gastroc_normalized.png", plt, dpi=300)


# loading scores (ranked genes)
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores)
gene_scores_ranked <- sort(gene_scores, decreasing = TRUE)
top_10_genes <- names(gene_scores_ranked[1:10])

pca$rotation[top_10_genes, 1] # with +/- sign