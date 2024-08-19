library(Biobase)
library(GenomicRanges)
library(SummarizedExperiment)
library(DESeq2)
library(htmltools)
library(ggplot2)
library(EnhancedVolcano)
library(magrittr)
library(BiocParallel)
register(SnowParam(4)) #MulticoreParam() not supported on Windows, use SnowParam()
library(IHW)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)
library(readr)


countData <- read.csv('rawcounts.csv', header = TRUE, sep = ",")
metaData <- read.csv('metadata.csv', header = TRUE, sep = ",")

head(countData)
metaData


dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~genotype, tidy = TRUE)
dds

## Pre-filtering ##
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$condition <- factor(dds$genotype, levels = c("Wild_type","sEH_KO"))

## Differential expression analysis ##
dds <- DESeq(dds)
res <- results(dds)
res
resultsNames(dds)
res <- results(dds, name="genotype_Wild_type_vs_sEH_KO")
res <- results(dds, contrast=c("genotype","Wild_type","sEH_KO"))

## Log fold change shrinkage for visualization and ranking ##
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="genotype_Wild_type_vs_sEH_KO", type="apeglm")
resLFC

## p-values and adjusted p-values ##
resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

## Independent hypothesis weighting ##
# (unevaluated code chunk)
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult
## Exploring and exporting results ##
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")

par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

## Plot counts ##
plotCounts(dds, gene=which.min(res$padj), intgroup="genotype")
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="genotype", 
                returnData=TRUE)

ggplot(d, aes(x=genotype, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

mcols(res)$description

## Exporting results to CSV files ##
write.csv(as.data.frame(resOrdered), 
          file="Wild_type_vs_sEH_KO_results.csv")
resSig <- subset(resOrdered, padj < 0.1)
resSig
## Data transformations and visualization ##
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

## Effects of transformations on the variance ##
# this gives log2(n + 1)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

## Data quality assessment by sample clustering and visualization ##

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)[("genotype")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

## Heatmap of the sample-to-sample distances ##
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

## Principal component plot of the samples ##
plotPCA(vsd, intgroup=("genotype"))



################################################################################
################################################################################
######################### Volcano Plots ########################################
################################################################################
################################################################################

ens <- rownames(dds)
## Enhanced Volcano ##
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'WT vs sEH-KO',
                pCutoff = 10e-2,
                FCcutoff = 0.49,
                pointSize = 3.0,
                labSize = 6.0)
rownames(res)
EnhancedVolcano(res,
                lab = rownames(res),
                selectLab = c("Nlrp3", "Ephx2", "Cxcl3", "Cxcl1", 
                              "Crip2", "Inka2", "Cx3cr1", "Cxcl1", 
                              "Zmat3", "Gpr84", "Ptgfrn", "Sorl1",
                              "Il1b", "Slc19a2"),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 10e-1,
                FCcutoff = 0.49,
                pointSize = 1.5,
                labSize = 5,
                title = 'DESeq2 results',
                subtitle = 'Differential expression',
                caption = 'FC cutoff, 0.49; p-value cutoff, 10e-2',
                legendPosition = "none",
                legendLabSize = 14,
                col = c('grey30', 'grey30', 'grey30', 'red2'),
                colAlpha = 0.9,
                drawConnectors = TRUE,
                widthConnectors = 0.3,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                hline = NULL,
                cutoffLineType = "blank")
pacman::p_load(here,  
               tidyverse, 
               janitor, # Cleaning column names  
               scales, # Transform axis scales   
               ggrepel)

vol_plot <- as.data.frame(res) %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(pvalue))) + 
  geom_point() 
vol_plot
# Create new categorical column ------------------------------------------------ 
data_all<- as.data.frame(res) %>%
  mutate(gene_type = case_when(log2FoldChange >= 0.49 & pvalue <= 0.1 ~ "Wild-type",
                               log2FoldChange <= -0.49 & pvalue <= 0.1 ~ "sEH-KO",
                               TRUE ~ "ns"))
data_all %>%
  count(gene_type)
wild_type <- data_all %>%
  filter(gene_type %in% c("Wild-type"))
sEH_KO <- data_all %>%
  filter(gene_type %in% c("sEH-KO"))

# Check gene_type categories ---------------------------------------------------
data_all %>%
  distinct(gene_type) %>%
  pull()  
# Add colour, size and alpha (transparency) to volcano plot --------------------
cols <- c("sEH-KO" = "darkblue", "Wild-type" = "darkgreen", "ns" = "grey") 
sizes <- c("sEH_KO" = 1, "Wild-type" = 1, "ns" = 0.5) 
alphas <- c("sEH_KO" = 1, "Wild-type" = 1, "ns" = 1)

# Layer more subplots ----------------------------------------------------------

spec_wt_genes <- data_all %>%
  filter(rownames(data_all) %in% c("Ephx2", "Crip2", "Ptgfrn", "Zmagt3", "Inka2", "Slc18a2"))

spec_sehko_genes <- data_all %>%
  filter(rownames(data_all) %in% c("Sorl1", "Gpr84", "Nlrp3", "Il1b", "Cxcl3", "Cx3cr1", "Cxcl1"))



final_plot <- ggplot(data = data_all,
                     aes(x = log2FoldChange ,
                         y = -log10(pvalue))) + 
  geom_point(aes(colour = gene_type), 
             alpha = 0.6, 
             shape = 16,
             size = 1) +
  geom_point(data = spec_wt_genes,
             shape = 21,
             size = 2,
             alpha = 1, 
             fill = "darkgreen", 
             colour = "black") +
  geom_text_repel(data = spec_wt_genes,   
                  aes(label = rownames(spec_wt_genes)),
                  force = 2, nudge_y = 0.6, nudge_x = 0.5) +
  geom_point(data = spec_sehko_genes,
             shape = 21,
             size = 2,
             alpha = 1, 
             fill = "darkblue", 
             colour = "black") +
  geom_text_repel(data = spec_sehko_genes,   
                  aes(label = rownames(spec_sehko_genes)),
                  force = 2, nudge_x = -0.5) +
  scale_colour_manual(values = cols) + 
  scale_x_continuous(breaks = c(seq(-10, 10, 2)),     
                     limits = c(-6, 6)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text = element_text(colour = "black"))

final_plot


## Exporting results to CSV files ##
write.csv(as.data.frame(res), 
          file="WT_vs_sEH2.csv")

## HeatMap ##

vsd <- vst(dds, blind=FALSE, fitType = "local", nsub = 100)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("genotype","id")])
pheatmap(mat, 
         annotation_col = anno, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))
         (50),
         fontsize = 10)
