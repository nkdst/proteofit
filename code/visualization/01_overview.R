## --------------------------- ##
##
## Purpose of script: creates overview visualizations for the data
##    currently includes
##      - readcounts
##
## Author: Nick Diercksen
##
## Date Created: 2022-09-15
##
## --------------------------- ##


## libraries
library(ggplot2)
library(treemap)
library(dplyr)

## including other file(s)
source("code/load_data.R")




gene_biotypes <- table(counts.gastroc$gene_biotype) %>%
  as.data.frame() %>%
  arrange(desc(Freq))

gene_biotypes[,gene_biotypes$Freq > 100]


## treemap
treemap::treemap(gene_biotypes,
                 index = "Var1",
                 vSize = "Freq",
                 title = "Gene biotypes")

# new label
# iris%>%
#   group_by(Species)%>%
#   summarise(Sum.Sepal.Length=sum(Sepal.Length))%>%
#   mutate(Species.Index=paste(Species, Sum.Sepal.Length, sep ="\n"))%>%
#   treemap(index="Species.Index", vSize="Sum.Sepal.Length")

## bar plots
# gene types
ggplot(gene_biotypes, aes(x = reorder(Var1, -Freq), y = Freq)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x=element_text(color = "black", angle=90, vjust=.8, hjust=0.8)) + 
  labs(x = "gene biotype", y = "frequency")


# count table matrices
counts.gastroc.m <- counts.gastroc[, 1:13] %>%
  t() %>%
  janitor::row_to_names(row_number = 1)

total_readcounts.gastroc <- colSums(counts.gastroc[, 2:13]) %>%
  as.data.frame() %>%
  rename(reads = ".") %>%
  tibble::rownames_to_column()

# total read counts
ggplot(total_readcounts.gastroc, aes(x = rowname, y = reads)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x=element_text(color = "black", angle=45, vjust=.8, hjust=0.8)) + 
  labs(x = "sample", y = "reads", title = "total read counts gastroc")


