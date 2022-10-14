## --------------------------- ##
##
## Purpose of script: creating a pca plot(s)
## 
## Note: c
##
## Author: Nick Diercksen
##
## Date Created: 2022-09-19
##
## --------------------------- ##

# libraries ----
library(dplyr)
library(janitor)
library(ggplot2)
library(ggfortify)
library(ggrepel)



# loading data ----
count <- readr::read_tsv("./data/readcount_genename_gastroc.xls") %>%
  tibble::column_to_rownames(var = "gene_id")

# making count table matrices
counts.gastroc.m <- counts.gastroc[, 1:13] %>%
  t() %>%
  janitor::row_to_names(row_number = 1)
class(counts.gastroc.m) <- "numeric"



rpm.gastroc <- colSums(counts.gastroc[, 2:13]) %>%
  as.data.frame() %>%
  rename(reads = ".") %>%
  tibble::rownames_to_column() %>%
  mutate(reads = reads / 10^6) # RPM (Reads Per)


# Normalization ----
# normalized per total readcount (sample wide)
counts.gastroc.norm <- sweep(counts.gastroc.m,1,rpm.gastroc$reads,FUN="/")
# TODO: do the counts need to be normalized by the read depth (gene length)?
# needs to be done before read count normalization to get TPM metric


# PCA ----
pca <- prcomp(counts.gastroc.norm)

pca.data <- data.frame(Sample = rownames(pca$x),
                       X = pca$x[, 1],
                       Y = pca$x[, 2]) %>%
  mutate(type = substr(Sample, 1,2))

# TODO: how to get the label 'a' removed from the legend?
plt <- autoplot(pca, data = pca.data, colour = 'type', label.repel=TRUE) +
  ggtitle("PCA gastroc normalized")
ggsave(filename = "./plots/02_pca_gastroc_normalized.png", plt, dpi=300)


## loading scores (ranked genes) ----
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores)
gene_scores_ranked <- sort(gene_scores, decreasing = TRUE)
top_10_genes <- names(gene_scores_ranked[1:10])


pca$rotation[top_10_genes, 1] # with +/- sign


# MDS ----
counts.gastroc.norm

## first, calculate the distance matrix using the Euclidian distance.
## NOTE: We are transposing, scaling and centering the data just like PCA.
distance.matrix <- dist(scale(counts.gastroc.norm, center=TRUE, scale=TRUE),
                        method="euclidean")

## do the MDS math (this is basically eigen value decomposition)
mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE) # classical multidimensional scaling

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
mds.var.per

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values <- mds.stuff$points
mds.data <- data.frame(Sample=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2])
mds.data


ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  theme_bw() +
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  ggtitle("MDS plot using Euclidean distance")

## TODO: outsource for cluster computing ----
# bigger the RAM and computing power are needed
# https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached
# usethis::edit_r_environ()

# average of the absolute value of the log fold change
# plotMDS()
# log2.data.matrix <- log2(counts.gastroc.norm)
# log2.distance.matrix <- matrix(
#   0,
#   nrow = ncol(log2.data.matrix),
#   ncol = ncol(log2.data.matrix)
# )
# 
# 
# 
# for(i in 1:ncol(log2.distance.matrix)) {
#   for(j in 1:i) {
#     log2.distance.matrix[i, j] <-
#       mean(abs(log2.data.matrix[,i] - log2.data.matrix[,j]))
#   }
# }
# 
# 
# ## do the MDS math (this is basically eigen value decomposition)
# ## cmdscale() is the function for "Classical Multi-Dimensional Scaling"
# mds.stuff <- cmdscale(as.dist(log2.distance.matrix),
#                       eig = TRUE,
#                       x.ret = TRUE)
# 
# ## calculate the percentage of variation that each MDS axis accounts for...
# mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
# mds.var.per