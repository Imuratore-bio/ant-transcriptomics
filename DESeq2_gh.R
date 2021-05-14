#partly adapted from http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

library(tximport)
library(readr)

cts <- as.matrix(read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/kallisto/no_low_counts_ensembl_cDNA.csv",row.names=1, check.names=FALSE))
coldata <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/kallisto/no_low_info_size_colonies.csv", row.names=1)
#examine the count matrix and column data to see if they are consistent in terms of sample order

library(DESeq2)

cts2 <- round(cts)

dds <- DESeqDataSetFromMatrix(countData = cts2,
                              colData = coldata,
                              design = ~ type + condition)
#first term in formula is controlled for, second is tested

dds$condition <- factor(dds$condition, levels = c("minum","media","major"))
dds$type <- factor(dds$type, levels = c("Ac22","M1","M2"))

dds <- DESeq(dds)

normcounts <- counts(dds, normalized=T)

write.csv(normcounts,'C:/Users/imura/Documents/grad_4/seq_analysis/deseq2_normcounts.csv')

res <- results(dds)

#change name(s) based on desired comparison, order matters
res <- results(dds, name="condition_major_vs_minum")
res <- results(dds, contrast=c("condition","major","minum"))
#res
#Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="condition_major_vs_minum", type="apeglm")
#resLFC

#Independent hypothesis weighting
library(IHW)
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)

#MA plot
#plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-3,3))

#Alternative shrinkage estimators
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")

#par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

#Plot counts
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
#plotCounts(dds, gene="XM_012208190.1", intgroup="condition")
#XM_012198788.1 - inscutable
#plotCounts(dds, gene="XM_012207394.1", intgroup="condition")
#growth factor activity

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

#Data transformations and visualization
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

ntd <- normTransform(dds)

library(vsn)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

library(pheatmap)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#for a subset of genes, in this case sensory-related only
'dds_order <- dds[, c("ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S1_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S2_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S6_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S10_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S13_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S16_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S20_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S21_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S25_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S29_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S3_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S4_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S7_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S8_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S11_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S18_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S19_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S26_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S28_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S5_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S9_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S12_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S22_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S23_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S24_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S27_L001_R1_23", 
                     +                                    "ensembl_cDNA_output_IDs_merged_merged_merged_LIB042574_TRA00159635_S30_L001_R1_23")]
'

#most significant genes
mat <- assay(vsd)[ head(order(res$padj),30), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[,c("type","condition")])
pheatmap(mat, fontsize=7, annotation_col=df)

#Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
library(RColorBrewer)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#Principal component plot of the samples
plotPCA(vsd, intgroup=c("condition", "type"))

pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_polygon(aes(fill = condition), alpha = 0.3) +
  coord_fixed()
#+   geom_text(label = rownames(pcaData), size = 2) 
#to see sample IDs

#JF style, with ellipses
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3, aes(shape=type)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  stat_ellipse(geom="polygon",level=0.65, alpha = 1/5, aes(fill = condition),show.legend = NA) +
  coord_fixed()

#write.csv(resLFC,'C:/Users/imura/Documents/grad_4/seq_analysis/no_low_padj_maj_min.csv'