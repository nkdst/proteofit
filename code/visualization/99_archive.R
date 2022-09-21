## --------------------------- ##
##
## Purpose of script:
##
## Author: Nick Diercksen
##
## Date Created: 2022-09-21
##
## --------------------------- ##


#### PCA ####
# normalized
counts.gastroc.norm <- sweep(counts.gastroc.m,1,total_readcounts.gastroc$reads,FUN="/")

pca <- prcomp(counts.gastroc.norm)

# plot(pca$x[,1], pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

pca.data <- data.frame(Sample = rownames(pca$x),
                       X = pca$x[, 1],
                       Y = pca$x[, 2]) %>%
  mutate(type = substr(Sample, 1,2))

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) + 
geom_text() + 
  geom_point() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_bw() + 
  ggtitle("PCA gastroc normalized")