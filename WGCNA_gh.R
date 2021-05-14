#WGCNA
#adapted partly from http://pklab.med.harvard.edu/scw2014/WGCNA.html
#install.packages(c("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel") ) 
source("http://bioconductor.org/biocLite.R") 
#biocLite("impute")
#install.packages("WGCNA")

library(WGCNA)

library(tximport)
library(readr)
library(DESeq2)


counts = read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/kallisto/no_low_counts_ensembl_cDNA.csv",row.names=1, check.names=FALSE)

cts <- as.matrix(read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/kallisto/no_low_counts_ensembl_cDNA.csv",row.names=1, check.names=FALSE))
cts2 <- round(cts)

coldata <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/kallisto/no_low_info_size_colonies.csv", row.names=1)
#examine the count matrix and column data to see if they are consistent in terms of sample order


dds <- DESeqDataSetFromMatrix(countData = cts2,
                              colData = coldata,
                              design = ~ type + condition)

dds$condition <- factor(dds$condition, levels = c("minim","media","major"))
dds$type <- factor(dds$type, levels = c("Ac22","M1","M2"))

dds <- DESeq(dds)
ddv = vst(dds)
write.csv(assay(ddv), "C:/Users/imura/Documents/grad_4/seq_analysis/no_low_merged_transformed_counts.csv")

#need to re-run above lines if switching counts data
dat = read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/no_low_merged_transformed_counts.csv")
rownames(dat) = dat$X
dat$X = NULL
tdat = as.data.frame(t(dat))

gsg = goodSamplesGenes(tdat)
gsg$allOK

#if FALSE
if (!gsg$allOK)
{if (sum(!gsg$goodGenes)>0)
  printFlush(paste("Removing genes:", paste(names(tdat)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(tdat)[!gsg$goodSamples], collapse=", ")))
  tdat= tdat[gsg$goodSamples, gsg$goodGenes]
}

#check
gsg=goodSamplesGenes(tdat, verbose = 3)
gsg$allOK 

traitData <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/kallisto/no_low_WGCNA_info_size_colonies.csv", row.names=1)
knitr::kable(head(traitData), caption = "Head of experimental Traits") 

rownames(traitData) = traitData$sample
traitData$sample = NULL

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

gene.names=rownames(traitData)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft=pickSoftThreshold(tdat,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")

cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence")) 
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")

abline(h=0.80,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#based on above graphs
softPower = 4;

#calclute the adjacency matrix
adj= adjacency(tdat,type = "unsigned", power = softPower);

cor <- WGCNA::cor
bwnet = blockwiseModules(adj, maxBlockSize = 12000,
                         power = 4, TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "Atta_blockwise",
                         verbose = 3)

minModuleSize = 30;

dynamicMods = cutreeDynamic(dendro = bwnet$dendrograms[[1]],  method="tree", minClusterSize = minModuleSize);

dynamicColors = labels2colors(dynamicMods)

plotDendroAndColors(bwnet$dendrograms[[1]], dynamicColors, main = "Gene dendrogram and module colors in block 1", dendroLabels = FALSE)

MEList = moduleEigengenes(tdat, colors = dynamicColors)
MEs = MEList$eigengenes

#Define number of genes and samples
nGenes = ncol(tdat)
MEDiss = 1-cor(MEs);
nSamples = nrow(tdat)
#Recalculate MEs with color labels
MEs0 = moduleEigengenes(tdat, dynamicColors)$eigengenes
MEs = orderMEs(MEs0)

#cmodule merging
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(tdat, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors

plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))

#as.data.frame(MEs0) to see contents ME for samples, moduleTraitCor for groups
write.csv(as.data.frame(MEs0),'C:/Users/imura/Documents/grad_4/seq_analysis/no_low_25merged_sample_module_correlation.csv')

moduleTraitCor = cor(MEs, traitData, use= "p")
write.csv(moduleTraitCor,'C:/Users/imura/Documents/grad_4/seq_analysis/no_low_25merged_module_trait_correlation.csv')

#module trait significance
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
write.csv(moduleTraitPvalue,'C:/Users/imura/Documents/grad_4/seq_analysis/no_low_25merged_module_trait_pvals.csv')


#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(2, 8, 2, 2))

#display the corelation values with a heatmap plot
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(traitData),
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.5,
               zlim= c(-.1,.1),
               main= paste("Module-trait relationships"))

#gene correlation to modules
datKME=signedKME(tdat, MEs, outputColumnName="ME.")
write.csv(datKME,'C:/Users/imura/Documents/grad_4/seq_analysis/no_low_25merged_gene_module_correlation.csv')

#print out genes in each module, do this for each color module of interest
color_names <- row.names(datKME[which(abs(datKME$ME.greenyellow) > 0.8),])
write.csv(color_names,'C:/Users/imura/Documents/grad_4/seq_analysis/no_low_25merged_lightgreen.csv')

#module scatter plot
library(plyr)
module = "greenyellow"
modNames = substring(names(MEs), 3)
column = match(module, modNames);
moduleColors = labels2colors(bwnet$colors)
moduleGenes = moduleColors==module;
geneModuleMembership = as.data.frame(cor(tdat, MEs, use = "p"));

weight = as.data.frame(revalue(coldata$condition, c("minim"=1, "media"=2, "major"=3)));
geneTraitSignificance = as.data.frame(cor(tdat, weight, use = "p"));
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for subcaste identity",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")

#gene trait significance
write.csv(geneTraitSignificance,'C:/Users/imura/Documents/grad_4/seq_analysis/no_low_25merged_geneTrait_pvals.csv')