# RNA-seq-analysis
RNA-seq analysis for differential gene expression (for R)
#prepare packages
library(DESeq2)
library(tidyverse)
install.packages("readxl")
library(readxl)

#load data

Gene_counts <- read.csv("Gene_counts", header= T, row.names = 1)

#set gene_name as row
Gene_counts[ duplicated(Gene_counts[,1]), 1]
row.names(Gene_counts) <- Gene_counts$gene_name
Gene_counts <- Gene_counts[!duplicated(Gene_counts$gene_name), ]

#import coldata

library (readxl)
sample_info <- read_excel("Coldata.xlsx")

sample_info$`Condition`<-factor(sample_info$`Condition`)


#Convert data frame into matrix count
true_matrix <- as.matrix(Gene_counts[, -1])  # omit the first column which contains Group names
rownames(true_matrix) <- Gene_counts$gene_name
true_matrix

#Round matrix count
true_matrix <- round(true_matrix)

#set Condition as design formula
design_formula <- ~ Condition

#Set factor levels
design_formula$Condition <- factor(design_formula$Condition, levels = c("Control", "Case"))

dds <- DESeqDataSetFromMatrix(countData = true_matrix, colData = sample_info, design = design_formula)

keep <- rowSums(counts(dds)) >=10

ddsDE <- DESeq(dds)

#Export Normalised counts

normCounts <- counts(ddsDE, normalized = T)

write.csv(normCounts,"first.normalised.case.v.control.csv")

#Solve formatting issue of normcount data

#DESeq results
res <- results (ddsDE, alpha = 0.05)

resOrdered <-res[order(res$padj),]
write.csv(resOrdered, "case.v.control.csv")

#summary res


#make plot in console "plotMA(ddsDE,ylim= c(-5,5))"


#PART TWO PLOTTING

library(ggplot2)
library(pheatmap)
install.packages("pheatmap")

#Load DeSeq data (deSeqRes and normCounts) - normCounts already readily available in R

deSeqRes <- read.csv("case.v.control.csv",row.names = 1)

View(normCounts)

deSeqRes$sig <- ifelse(deSeqRes$padj <=0.05, "yes", "no")


#make plot with ggplot
ggplot(deSeqRes, aes(x= log10(baseMean), y= log2FoldChange, color= sig)) +
  geom_point()


#remove the NA from ggplot
deSeqRes <- na.omit(deSeqRes)
ggplot(deSeqRes, aes(x= log10(baseMean), y= log2FoldChange, color= sig)) +
  geom_point()

#pheatmap

signi <- subset(deSeqRes, padj <=0.05)

allSig <-merge(normCounts, signi,by=0)

sigCounts <- allSig[,2:12]
row.names(sigCounts) <-allSig$Row.names

#pheatmap - write in console pheatmap(log2(sigCounts +1), scale= 'row', show_rownames = F)
