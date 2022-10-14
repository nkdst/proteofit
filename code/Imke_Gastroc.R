# analyse dataset seperately. 
# combine 2 of them 
install.packages("readxl")
install.packages("DESeq2")
install.packages("RColorBrewer")
install.packages("pheatmap")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("tidyr")

# load library ####
library(readxl)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(clusterProfiler)
library(org.Mm.eg.db)



# load data ####
count <- read_tsv("readcount_genename.xls")%>%
  column_to_rownames(var = "gene_id")

# head
head(count, 6)

# Create data frame
metadata <- read_csv2("sample_table.csv")%>%
  column_to_rownames("run")

reorder_index <- match(rownames(metadata), colnames(count))

count <- count[,reorder_index]

# Check that all of the samples are in the same order in the metadata and count data
all(rownames(metadata) %in% colnames(count))

# Create a DESeq2 object for  tissue, for paired analysis use patient
dds <- DESeqDataSetFromMatrix(countData =  count,
                              colData =  metadata,
                              design = ~ genotype)

# Determine the size factors to use for normalization
dds <- estimateSizeFactors(dds)

# Extract the normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
# OR this will give you unnormalized values
unNormalized_counts <- counts(dds, normalized = FALSE)

# unsupervised clustering analyses with heatmap and PCA of normalized counts
# before visulization you have to log transform it. 

# Transform the normalized counts (vst is blind to sample information)
vsd <- vst(dds, blind = T)

# Create heatmap of sample correlation values
vsd %>% 
  assay() %>% # Extract the matrix of transformed counts
  cor() %>% # Compute the correlation values between samples
  pheatmap(annotation = dplyr::select(metadata, c("genotype")))

# Plot the PCA of PC1 and PC2 and color by siRNA 
library(RColorBrewer)
n_cols <- 11
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(n_cols)

pcaData <- plotPCA(vsd, intgroup=  c("genotype"),
                   returnData = T)
# normal PCA
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  scale_color_manual(values =  mycolors)
# grouped PCA
library(ggfortify)
pcaData_2 <- prcomp(pcaData[,1:2])

autoplot(pcaData_2, data =pcaData , colour = "genotype",frame = TRUE)+
  geom_text(aes(label = genotype), stat = "unique", size = 3, angle =0)+
  scale_fill_manual(values =  mycolors)+
  scale_color_manual(values =  mycolors)+
  theme_classic()+
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.position = "none"
  )


# If outlier in PCA and heatmap, remove it and do the DESeq2 object again

# Run the DESeq2 analysis
dds <- DESeq(dds)

# Plot dispersions
plotDispEsts(dds)
resultsNames(dds)

# result Pmepa1_Scramble ####
result <- results(dds,
                           contrast = c("genotype", "KO", "WT")) # The results function of the DESeq2 package performs independent filtering by default using the mean of normalized counts as a filter statistic.


# MA plot
plotMA(result, ylim=c(-8,8))

# MA plot after schrinkage
# Shrink the log2 fold change estimates to be more accurate
result <- lfcShrink(dds, 
                             coef = "genotype_WT_vs_KO")


plotMA(result, ylim=c(-5,5))

mcols(result)

head(result, n=10)

# summary of up and down genes
summary(result)


#
#library(annotable) 
grcm38 <- readRDS("grcm38.rds")

# Save results as a data frame
result_df <- data.frame(result) %>% rownames_to_column(var = "ensgene") %>%
  left_join(grcm38[, c("ensgene", "symbol", "description")], by = "ensgene")


# Subset the results to only return the significant genes with p-adjusted values less than 0.05

result_sig <- subset(result, padj < 0.1) %>%
  data.frame() %>%
  arrange(padj)%>%
  rownames_to_column()%>%
  left_join(grcm38 [, c("ensgene", "symbol", "description")],
            by = c("rowname" = "ensgene") ) %>% 
  arrange(log2FoldChange)

writexl::write_xlsx(result_sig, "result.xlsx")


# visulization of the result
# Subset normalized counts to E3 ligase genes
sig_normalized_counts <- normalized_counts[result_sig$rowname, ] %>% 
  data.frame() %>% 
  dplyr::select(starts_with(c("WT","KO")))

# Choose a color palette from RColorBrewer
library(RColorBrewer)
heat_colors <- brewer.pal(6, "YlOrRd")

display.brewer.all()

# Run pheatmap
pheatmap(sig_normalized_counts, 
         color = heat_colors,
         cluster_rows = T,
         show_rownames = F,
         annotation = select(metadata, genotype), 
         scale = "row")



# Obtain logical vector regarding whether padj values are less than 0.05
result <- result %>%
  data.frame()%>%
  rownames_to_column(var = "ensmbl") %>%
  mutate(threshold = padj < 0.1)%>%
  na.omit()


# Volcano plot

ggplot(result) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold))+
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25)))

# plot top 20 down
myColor <- RColorBrewer::brewer.pal(8, "Dark2")

result_up <- subset(result, padj < 0.1) %>%
  data.frame() %>%
  arrange(desc(log2FoldChange))%>%
  filter(baseMean > 100)%>%
  top_n( 20, log2FoldChange)

up_normalized_counts <- normalized_counts[result_up$ensmbl, ]


result_up <- data.frame(up_normalized_counts)%>% rownames_to_column(var = "ensembl")%>%
  dplyr::select(ensembl,starts_with(c("WT", "KO"))) %>% 
  gather(key = "samplename",
         value = "normalized_counts", 2:12)%>%
  inner_join(rownames_to_column(metadata, var = "samplename"),
             by = "samplename")%>%
  left_join(grcm38[, c("ensgene", "symbol", "description")],
            by = c("ensembl" = "ensgene")) %>% 
  mutate(genotype = factor(genotype)) %>% 
  mutate(genotype = relevel(genotype, "WT")) %>% 
  mutate(symbol = fct_reorder(symbol, desc(normalized_counts))) %>% 
  na.omit()


ggplot(result_up, aes(x = symbol, y = normalized_counts, fill = genotype)) +
  geom_dotplot(binaxis='y', stackdir='center',
               dotsize = 0.5, position=position_dodge(0.8))+
  geom_boxplot()+
  scale_y_log10() +
  xlab("Genes") +
  scale_fill_manual(values = myColor)+
  ylab("Normalized Counts") +
  ggtitle("Top 20 Down-regulated Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

# plot up 20 genes
result_down <- subset(result, padj < 0.1) %>%
  data.frame() %>%
  arrange(desc(log2FoldChange))%>%
  filter(baseMean > 100)%>%
  top_n(- 20, log2FoldChange)

down_normalized_counts <- normalized_counts[result_down$ensmbl, ]


result_down <- data.frame(down_normalized_counts)%>% rownames_to_column(var = "ensembl")%>%
  dplyr::select(ensembl,starts_with(c("WT", "KO"))) %>% 
  gather(key = "samplename",
         value = "normalized_counts", 2:12)%>%
  inner_join(rownames_to_column(metadata, var = "samplename"),
             by = "samplename")%>%
  left_join(grcm38[, c("ensgene", "symbol", "description")],
            by = c("ensembl" = "ensgene")) %>% 
  mutate(genotype = factor(genotype)) %>% 
  mutate(genotype = relevel(genotype, "WT")) %>% 
  mutate(symbol = fct_reorder(symbol, desc(normalized_counts))) %>% 
  na.omit()



ggplot(result_down, aes(x = symbol, y = normalized_counts, fill = genotype)) +
  geom_dotplot(binaxis='y', stackdir='center',
               dotsize = 0.5, position=position_dodge(0.8))+
  geom_boxplot()+
  scale_y_log10() +
  xlab("Genes") +
  ylab("Normalized Counts") +
  scale_fill_manual(values = myColor)+
  ggtitle("Top 20 Up-regulated Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))
  

# GO analysis with top 20 UP and Down genes upon Pmepa1
result_GO <- result %>%
  filter(padj < 0.01 &
           baseMean > 100) %>% 
  left_join(grcm38[, c("ensgene", "symbol", "description", "entrez")], by = c("ensmbl" ="ensgene"))

result_GO<-enrichGO(result_GO$entrez, "org.Mm.eg.db", ont= "BP", pAdjustMethod = "bonferroni", readable = T)


GO_simplified  <- simplify(
  result_GO,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL)

dotplot(result_GO,showCategory=15, font.size = 12,title= "GO_Gastroc")

result_GO_df <- data.frame(result_GO)


# Proteosome pathway genes upon Nfe2l1 KD
myColor <- RColorBrewer::brewer.pal(8, "Dark2")[3]

UPS_result <-result_GO_df %>% 
  dplyr::filter(Description == "mitochondrion organization") %>% 
  dplyr::select(geneID) %>% 
  separate(geneID,into = paste("gene", 1:220, sep = "/")) %>% 
  gather(key = key, value =  symbol, 1:220)%>%
  right_join(grcm38[,c("symbol", "ensgene")], by = "symbol" ) %>% 
  na.omit() %>% 
  distinct()

normalized_counts_df <- data.frame(normalized_counts)
# plot UPS
UPS_normalized_counts <- normalized_counts_df[UPS_result$ensgene, ] %>% 
  data.frame() %>%
  rownames_to_column(var = "ensembl")%>%
  gather(key = "samplename",
         value = "normalized_counts", 2:12)%>%
  inner_join(rownames_to_column(metadata, var = "samplename"),
             by = "samplename")%>%
  left_join(grcm38[, c("ensgene", "symbol", "description")],
            by = c("ensembl" = "ensgene")) %>% 
  mutate(genotype = factor(genotype)) %>% 
  mutate(genotype = relevel(genotype, "WT")) %>% 
  mutate(symbol = fct_reorder(symbol, desc(normalized_counts))) %>% 
  group_by_at(vars(symbol)) %>% 
  left_join(
    x = filter(., genotype != ""),
    y = filter(., genotype == "WT")%>% 
      group_by(symbol) %>% 
      summarise(normalized_counts = mean(normalized_counts)),
    by = "symbol") %>% 
  mutate(fold_change = normalized_counts.x / normalized_counts.y) %>% 
  dplyr::select(- ends_with(".y")) %>% 
  filter (fold_change > 2.5 | fold_change < 0.35)    %>% 
  filter(genotype!= "WT") %>% 
  ungroup() %>% 
  mutate(symbol = fct_reorder(symbol, desc(fold_change))) %>% 
  ggplot(aes(x = symbol, y = fold_change, fill = genotype)) +
  geom_dotplot(binaxis='y', stackdir='center',
               dotsize = 0.5, position=position_dodge(0.8))+
  geom_boxplot(alpha = 0.4)+
  scale_y_log10() +
  scale_fill_manual(values = myColor)+
  xlab("Genes") +
  ylab("Fold Change") +
  ggtitle("mitochondrion organization") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  geom_text(aes(21, 1.2), label = "WT")

UPS_normalized_counts



# Protein catabolic process ####

myColor <- RColorBrewer::brewer.pal(8, "Dark2")[3]

UPS_result <-result_GO_df %>% 
  dplyr::filter(Description == "proteasomal protein catabolic process") %>% 
  dplyr::select(geneID) %>% 
  separate(geneID,into = paste("gene", 1:220, sep = "/")) %>% 
  gather(key = key, value =  symbol, 1:220)%>%
  right_join(grcm38[,c("symbol", "ensgene")], by = "symbol" ) %>% 
  na.omit() %>% 
  distinct()
# plot UPS
UPS_normalized_counts <- normalized_counts[UPS_result$ensgene, ] %>% 
  data.frame() %>%
  rownames_to_column(var = "ensembl")%>%
  gather(key = "samplename",
         value = "normalized_counts", 2:13)%>%
  inner_join(rownames_to_column(metadata, var = "samplename"),
             by = "samplename")%>%
  left_join(grcm38[, c("ensgene", "symbol", "description")],
            by = c("ensembl" = "ensgene")) %>% 
  mutate(genotype = factor(genotype)) %>% 
  mutate(genotype = relevel(genotype, "WT")) %>% 
  mutate(symbol = fct_reorder(symbol, desc(normalized_counts))) %>% 
  group_by_at(vars(symbol)) %>% 
  left_join(
    x = filter(., genotype != ""),
    y = filter(., genotype == "WT")%>% 
      group_by(symbol) %>% 
      summarise(normalized_counts = mean(normalized_counts)),
    by = "symbol") %>% 
  mutate(fold_change = normalized_counts.x / normalized_counts.y) %>% 
  dplyr::select(- ends_with(".y")) %>% 
  filter (fold_change > 2.5 | fold_change < 0.325)    %>% 
  filter(genotype!= "WT") %>% 
  ungroup() %>% 
  mutate(symbol = fct_reorder(symbol, desc(fold_change))) %>% 
  ggplot(aes(x = symbol, y = fold_change, fill = genotype)) +
  geom_dotplot(binaxis='y', stackdir='center',
               dotsize = 0.5, position=position_dodge(0.8))+
  geom_boxplot(alpha = 0.4)+
  scale_y_log10() +
  scale_fill_manual(values = myColor)+
  xlab("Genes") +
  ylab("Fold Change") +
  ggtitle("proteasomal protein catabolic process") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  geom_text(aes(30, 1.2), label = "WT")

UPS_normalized_counts


# Proteasome ####


myColor <- RColorBrewer::brewer.pal(8, "Dark2")[3]

UPS_result <-result_GO_df %>% 
  dplyr::filter(Description == "proteasome-mediated ubiquitin-dependent protein catabolic process") %>% 
  dplyr::select(geneID) %>% 
  separate(geneID,into = paste("gene", 1:220, sep = "/")) %>% 
  gather(key = key, value =  symbol, 1:220)%>%
  right_join(grcm38[,c("symbol", "ensgene")], by = "symbol" ) %>% 
  na.omit() %>% 
  distinct()
# plot UPS
UPS_normalized_counts <- normalized_counts[UPS_result$ensgene, ] %>% 
  data.frame() %>%
  rownames_to_column(var = "ensembl")%>%
  gather(key = "samplename",
         value = "normalized_counts", 2:13)%>%
  inner_join(rownames_to_column(metadata, var = "samplename"),
             by = "samplename")%>%
  left_join(grcm38[, c("ensgene", "symbol", "description")],
            by = c("ensembl" = "ensgene")) %>% 
  mutate(genotype = factor(genotype)) %>% 
  mutate(genotype = relevel(genotype, "WT")) %>% 
  mutate(symbol = fct_reorder(symbol, desc(normalized_counts))) %>% 
  group_by_at(vars(symbol)) %>% 
  left_join(
    x = filter(., genotype != ""),
    y = filter(., genotype == "WT")%>% 
      group_by(symbol) %>% 
      summarise(normalized_counts = mean(normalized_counts)),
    by = "symbol") %>% 
  mutate(fold_change = normalized_counts.x / normalized_counts.y) %>% 
  dplyr::select(- ends_with(".y")) %>% 
  filter (fold_change > 2.5 | fold_change < 0.325)    %>% 
  filter(genotype!= "WT") %>% 
  ungroup() %>% 
  mutate(symbol = fct_reorder(symbol, desc(fold_change))) %>% 
  ggplot(aes(x = symbol, y = fold_change, fill = genotype)) +
  geom_dotplot(binaxis='y', stackdir='center',
               dotsize = 0.5, position=position_dodge(0.8))+
  geom_boxplot(alpha = 0.4)+
  scale_y_log10() +
  scale_fill_manual(values = myColor)+
  xlab("Genes") +
  ylab("Fold Change") +
  ggtitle("proteasome-mediated ubiquitin-dependent protein catabolic process") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  geom_text(aes(27, 1.2), label = "WT")

UPS_normalized_counts

# energy metabolite ####


myColor <- RColorBrewer::brewer.pal(8, "Dark2")[3]

UPS_result <-result_GO_df %>% 
  dplyr::filter(Description == "generation of precursor metabolites and energy") %>% 
  dplyr::select(geneID) %>% 
 # dplyr::filter(geneID != "Aldoa")%>%
  separate(geneID,into = paste("gene", 1:130, sep = "/")) %>% 
  gather(key = key, value =  symbol, 1:130)%>%
  right_join(grcm38[,c("symbol", "ensgene")], by = "symbol" ) %>% 
  na.omit() %>% 
  distinct()

# plot UPS
UPS_normalized_counts <- normalized_counts[UPS_result$ensgene, ] %>% 
  data.frame() %>%
  rownames_to_column(var = "ensembl")%>%
  gather(key = "samplename",
         value = "normalized_counts", 2:13)%>%
  inner_join(rownames_to_column(metadata, var = "samplename"),
             by = "samplename")%>%
  left_join(grcm38[, c("ensgene", "symbol", "description")],
            by = c("ensembl" = "ensgene")) %>% 
  mutate(genotype = factor(genotype)) %>% 
  mutate(genotype = relevel(genotype, "WT")) %>% 
  mutate(symbol = fct_reorder(symbol, desc(normalized_counts))) %>% 
  group_by_at(vars(symbol)) %>% 
  left_join(
    x = filter(., genotype != ""),
    y = filter(., genotype == "WT")%>% 
      group_by(symbol) %>% 
      summarise(normalized_counts = mean(normalized_counts)),
    by = "symbol") %>% 
  mutate(fold_change = normalized_counts.x / normalized_counts.y) %>% 
  dplyr::select(- ends_with(".y")) %>% 
  filter (fold_change > 1.5 | fold_change < 0.375)    %>% 
  filter(genotype!= "WT") %>% 
  ungroup() %>% 
  mutate(symbol = fct_reorder(symbol, desc(fold_change))) %>% 
  ggplot(aes(x = symbol, y = fold_change, fill = genotype)) +
  geom_dotplot(binaxis='y', stackdir='center',
               dotsize = 0.5, position=position_dodge(0.8))+
  geom_boxplot(alpha = 0.4)+
  scale_y_log10() +
  scale_fill_manual(values = myColor)+
  xlab("Genes") +
  ylab("Fold Change") +
  ggtitle("generation of precursor metabolites and energy") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  geom_text(aes(24, 1.5), label = "WT")

UPS_normalized_counts

# add gene symbol
Pmepa1_geneID <- data.frame(Pmepa1_Scramble) %>% rownames_to_column(var = "ensgene") %>%
  left_join(grcm38[, c("ensgene", "symbol", "description", "entrez")], by = "ensgene")%>%
  na.omit(log2FoldChange)%>%
  filter(padj < 0.1)%>%
  filter(baseMean >= 6000 )%>%
  filter(abs(log2FoldChange)>= 0.2)

# KEGG
Pmepa1_KEGG<-enrichKEGG(Pmepa1_geneID$entrez, organism = "mmu", pAdjustMethod = "bonferroni",pvalueCutoff = 0.05)
dotplot(Pmepa1_KEGG, title= "KEGG analysis: Pmepa1_Scramble")


