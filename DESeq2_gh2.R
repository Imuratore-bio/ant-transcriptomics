library(tximport)
library(readr)

cts <- as.matrix(read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/kallisto/no_low_counts_ensembl_cDNA.csv",row.names=1, check.names=FALSE))
coldata <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/kallisto/no_low_info_size_colonies.csv", row.names=1)
#examine the count matrix and column data to see if they are consistent in terms of sample order
library(DESeq2)
#if it is complaining about integers:
cts2 <- round(cts)

dds <- DESeqDataSetFromMatrix(countData = cts2,
                              colData = coldata,
                              design = ~ type + condition)

dds$condition <- factor(dds$condition, levels = c("minum","media","major"))
dds$type <- factor(dds$type, levels = c("Ac22","M1","M2"))

dds <- DESeq(dds)

normcounts <- counts(dds, normalized=T)
#write.csv(normcounts,'C:/Users/imura/Documents/grad_4/seq_analysis/deseq2_normcounts.csv')
res <- results(dds)
#change name(s) based on desired comparison, order matters
res <- results(dds, name="condition_major_vs_minum")
res <- results(dds, contrast=c("condition","major","minum"))
#Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="condition_major_vs_minum", type="apeglm")

#Independent hypothesis weighting
library(IHW)
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)

#MA plot
plotMA(resLFC, ylim=c(-3,3))

#Alternative shrinkage estimators
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")

xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

#Plot counts
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

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

#most sig genes
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

#optional, to use all genes use these lines instead
#pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE, ntop=19000)
#percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_polygon(aes(fill = condition), alpha = 0.3) +
  coord_fixed()
#+   geom_text(label = rownames(pcaData), size = 2) 
#to see sample IDs

#different style, with ellipses
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3, aes(shape=type)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  stat_ellipse(geom="polygon",level=0.65, alpha = 1/5, aes(fill = condition),show.legend = NA) +
  coord_fixed() +
  theme_bw() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20)) 


library(ggfortify)

#loadings plot adapted from:
#https://tavareshugo.github.io/data-carpentry-rnaseq/03_rnaseq_pca.html
#https://tavareshugo.github.io/data-carpentry-rnaseq/00_exercises.html#33_Visualise_variable_loadings
library(tidyr)
library(dplyr)

rv <- rowVars(assay(vsd))
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(assay(vsd)[select,]))

pc_loadings <- pca$rotation

pc_loadings <- pc_loadings %>% 
  as_tibble(rownames = "gene")

top_genes <- pc_loadings %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_loadings <- pc_loadings %>% 
  filter(gene %in% top_genes)

#to replace with actual names from biomart
#top_loadings[top_loadings == "XM_012198828.1"] <- "12 kDa FK506-binding protein, transcript variant X2"
#top_loadings[top_loadings == "XM_012199022.1"] <- "PREDICTED: uncharacterized protein LOC105617468"
#top_loadings[top_loadings == "XM_012199072.1"] <- "cytochrome P450 4g15-like"
#top_loadings[top_loadings == "XM_012200229.1"] <- "probable sodium/potassium/calcium exchanger CG1090"
#top_loadings[top_loadings == "XM_012200257.1"] <- "15-hydroxyprostaglandin dehydrogenase [NAD(+)]-like"
#top_loadings[top_loadings == "XM_012200567.1"] <- "GRAM domain-containing protein 2-like"
#top_loadings[top_loadings == "XM_012200775.1"] <- "PREDICTED: uncharacterized protein LOC105619251"
#top_loadings[top_loadings == "XM_012204254.1"] <- "calsyntenin-1"
#top_loadings[top_loadings == "XM_012206379.1"] <- "probable tRNA(His) guanylyltransferase, transcript variant X2"
#top_loadings[top_loadings == "XM_012206719.1"] <- "pheromone-binding protein Gp-9-like"
#top_loadings[top_loadings == "XM_012207376.1"] <- "transcription elongation factor SPT5-like"
#top_loadings[top_loadings == "XM_012207686.1"] <- "fatty acid synthase"
#top_loadings[top_loadings == "XM_012207883.1"] <- "fatty acid synthase-like"
#top_loadings[top_loadings == "XM_012208118.1"] <- "probable cytochrome P450 305a1"
#top_loadings[top_loadings == "XM_012208245.1"] <- "leucine-rich repeats and immunoglobulin-like domains protein 2"
#top_loadings[top_loadings == "XM_012208925.1"] <- "puromycin-sensitive aminopeptidase-like"
#top_loadings[top_loadings == "XM_012209016.1"] <- "neuroglian, transcript variant X1"
#top_loadings[top_loadings == "XM_012209044.1"] <- "glycoprotein 3-alpha-L-fucosyltransferase A"
#top_loadings[top_loadings == "XM_012209210.1"] <- "prohormone-3"
#top_loadings[top_loadings == "XR_001046523.1"] <- "nuclear pore complex protein DDB_G0274915 homolog, transcript variant X2"

#abbreviated name versions
top_loadings[top_loadings == "XM_012198828.1"] <- "12 kDa FK506"
top_loadings[top_loadings == "XM_012199022.1"] <- "uncharacterized (LOC105617468)"
top_loadings[top_loadings == "XM_012199072.1"] <- "cytochrome"
top_loadings[top_loadings == "XM_012200229.1"] <- "CG1090"
top_loadings[top_loadings == "XM_012200257.1"] <- "15-hydroxyprostaglandin dehydrogenase"
top_loadings[top_loadings == "XM_012200567.1"] <- "GRAM domain"
top_loadings[top_loadings == "XM_012200775.1"] <- "uncharacterized (LOC105619251)"
top_loadings[top_loadings == "XM_012204254.1"] <- "calsyntenin-1"
top_loadings[top_loadings == "XM_012206379.1"] <- "guanylyltransferase"
top_loadings[top_loadings == "XM_012206719.1"] <- "Gp-9"
top_loadings[top_loadings == "XM_012207376.1"] <- "SPT5"
top_loadings[top_loadings == "XM_012207686.1"] <- "fatty acid synthase"
top_loadings[top_loadings == "XM_012207883.1"] <- "fatty acid synthase-like"
top_loadings[top_loadings == "XM_012208118.1"] <- "probable cytochrome"
top_loadings[top_loadings == "XM_012208245.1"] <- "leucine repeats/immunoglobulin"
top_loadings[top_loadings == "XM_012208925.1"] <- "aminopeptidase"
top_loadings[top_loadings == "XM_012209016.1"] <- "neuroglian"
top_loadings[top_loadings == "XM_012209044.1"] <- "glycoprotein"
top_loadings[top_loadings == "XM_012209210.1"] <- "prohormone-3"
top_loadings[top_loadings == "XR_001046523.1"] <- "DDB_G0274915"

library(ggrepel)
ggplot(data = top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               color = "gray") +
  geom_text_repel(aes(x = PC1, y = PC2, label = gene), force_pull = 0.1, size = 5, max.overlaps = Inf) +
  theme(text = element_text(size = 40)) +
  scale_x_continuous(expand = c(0.02, 0.02))  +
  xlab(paste0("PC1")) +
  ylab(paste0("PC2")) +
  theme_bw() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20))
