library(ggplot2)
library(dplyr)

#from DESeq2
#all genes
nrmcts <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/for_scatter_no_low_deseq2_normcounts.csv", row.names=1)

ggplot(nrmcts,aes(condition, AVERAGE)) +
  geom_smooth(method='lm') +
  theme_bw() +   theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Gene expression")

#restricted to significant (non-colony) genes only
sig <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/no_low_COLONYREMOVED_SIGONLY_padj_all_nonredundant_mixed.csv", row.names=1)

sig_nrmcts <- nrmcts[which(row.names(t(nrmcts)) %in% row.names(sig) | row.names(t(nrmcts)) %in% c("sample", "condition", "type") )]
sig_nrmcts$AVERAGE <- c(rowMeans(sig_nrmcts[ , -which(names(sig_nrmcts) %in% c("condition","type"))]))

ggplot(sig_nrmcts,aes(condition, AVERAGE)) +
  geom_smooth(method='lm') +
  theme_bw() +   theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Significant expression")

#restricted to positively size-correlated genes only

majmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/padj_maj_min_no_NA.csv", stringsAsFactors = FALSE)
majmed <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/padj_maj_med_no_NA.csv", stringsAsFactors = FALSE)
medmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/padj_med_min_no_NA.csv", stringsAsFactors = FALSE)

growth <- c()

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] > 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] > 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] > 0)) {
    
    growth = c(growth, majmin$X[row])
    
  }
}

growth_nrmcts <- nrmcts[which(row.names(t(nrmcts)) %in% growth | row.names(t(nrmcts)) %in% c("sample", "condition", "type") )]

growth_nrmcts$AVERAGE <- c(rowMeans(growth_nrmcts[ , -which(names(growth_nrmcts) %in% c("condition","type"))]))

ggplot(growth_nrmcts,aes(condition, AVERAGE)) +
  geom_smooth(method='lm') +
  theme_bw() +   theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Growth-linked expression")


#restricted to negatively size-correlated genes only
majmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/padj_maj_min_no_NA.csv", stringsAsFactors = FALSE)
majmed <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/padj_maj_med_no_NA.csv", stringsAsFactors = FALSE)
medmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/padj_med_min_no_NA.csv", stringsAsFactors = FALSE)

down <- c()

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] < 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] < 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] < 0)) {
    
    down = c(down, majmin$X[row])
    
  }
}

down_nrmcts <- nrmcts[which(row.names(t(nrmcts)) %in% down | row.names(t(nrmcts)) %in% c("sample", "condition", "type") )]

down_nrmcts$AVERAGE <- c(rowMeans(down_nrmcts[ , -which(names(down_nrmcts) %in% c("condition","type"))]))

ggplot(down_nrmcts,aes(condition, AVERAGE)) +
  geom_smooth(method='lm') +
  theme_bw() +   theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Inverse growth expression")

#restricted to media-linked genes only
majmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/padj_maj_min_no_NA.csv", stringsAsFactors = FALSE)
majmed <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/padj_maj_med_no_NA.csv", stringsAsFactors = FALSE)
medmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/padj_med_min_no_NA.csv", stringsAsFactors = FALSE)

ushape = c()

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] < 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] > 0)) {
    
    ushape = c(ushape, majmin$X[row])
    
  }
}

u_nrmcts <- nrmcts[which(row.names(t(nrmcts)) %in% ushape | row.names(t(nrmcts)) %in% c("sample", "condition", "type") )]

u_nrmcts$AVERAGE <- c(rowMeans(u_nrmcts[ , -which(names(u_nrmcts) %in% c("condition","type"))]))

ggplot(u_nrmcts,aes(condition, AVERAGE)) +
  #geom_smooth(method='lm') +
  geom_smooth(method='loess') +
  theme_bw() +   theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Media-biased expression")

#restricted to size-linked genes only
majmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/no_low_padj_maj_min_no_NA.csv")
majmed <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/no_low_padj_maj_med_no_NA.csv")
medmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/no_low_padj_med_min_no_NA.csv")

majmin_nogrowth = majmin
majmed_nogrowth = majmed
medmin_nogrowth = medmin

size <- c()

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] > 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] > 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] > 0)) {
    
    size <- c(size, majmin$X[row])
    
  }
}

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] < 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] < 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] < 0)) {
    
    size <- c(size, majmin$X[row])
    
  }
}

size_nrmcts <- nrmcts[which(row.names(t(nrmcts)) %in% size | row.names(t(nrmcts)) %in% c("sample", "condition", "type") )]

size_nrmcts$AVERAGE <- c(rowMeans(size_nrmcts[ , -which(names(size_nrmcts) %in% c("condition","type"))]))

ggplot(size_nrmcts,aes(condition, AVERAGE)) +
  geom_smooth(method='lm') +
  theme_bw() +   theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Size-linked expression")

#restricted to non-size-linked genes only
majmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/no_low_padj_maj_min_no_NA.csv")
majmed <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/no_low_padj_maj_med_no_NA.csv")
medmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/no_low_padj_med_min_no_NA.csv")

majmin_nogrowth = majmin
majmed_nogrowth = majmed
medmin_nogrowth = medmin

nonsize <- c()

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] > 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] > 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] > 0)) {
    
    nonsize <- c(nonsize, majmin$X[row])
    
  }
}

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] < 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] < 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] < 0)) {
    
    nonsize <- c(nonsize, majmin$X[row])
    
  }
}

nonsize_nrmcts <- nrmcts[which(!(row.names(t(nrmcts)) %in% nonsize) & (!row.names(t(nrmcts)) %in% c("AVERAGE")) )]

nonsize_nrmcts$AVERAGE <- c(rowMeans(nonsize_nrmcts[ , -which(names(nonsize_nrmcts) %in% c("condition","type"))]))

ggplot(nonsize_nrmcts,aes(condition, AVERAGE)) +
  geom_smooth(method='lm') +
  theme_bw() +   theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Non-size linked expression")

#SIGNIFICANT sensory genes only
sig_sensory_nrmcts <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/SENSORY_ONLY_for_scatter_no_low_deseq2_normcounts.csv", row.names=1)

ggplot(sig_sensory_nrmcts,aes(condition, AVERAGE)) +
  geom_smooth(method='lm') +
  theme_bw() +   theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Significant sensory expression")

#sensory genes only, regardless of significance
sensory <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/all_sensory_genes.csv", row.names=1)

sensory_nrmcts <- nrmcts[which(row.names(t(nrmcts)) %in% row.names(sensory) | row.names(t(nrmcts)) %in% c("sample", "condition", "type") )]
sensory_nrmcts$AVERAGE <- c(rowMeans(sensory_nrmcts[ , -which(names(sensory_nrmcts) %in% c("condition","type"))]))

ggplot(sensory_nrmcts,aes(condition, AVERAGE)) +
  geom_smooth(method='lm') +
  theme_bw() +   theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Sensory-related expression")

#sensory genes only, regardless of significance
sensory <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/all_sensory_genes.csv", row.names=1)

sensory_nrmcts <- nrmcts[which(row.names(t(nrmcts)) %in% row.names(sensory) | row.names(t(nrmcts)) %in% c("sample", "condition", "type") )]
sensory_nrmcts$AVERAGE <- c(rowMeans(sensory_nrmcts[ , -which(names(sensory_nrmcts) %in% c("condition","type"))]))

ggplot(sensory_nrmcts,aes(condition, AVERAGE)) +
  geom_smooth(method='lm') +
  #geom_smooth(method='loess') +
  theme_bw() +   theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Sensory-related expression")

#metabolism genes only, regardless of significance
meta <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/all_metabolism_etaboli_genes.csv", row.names=1)

meta_nrmcts <- nrmcts[which(row.names(t(nrmcts)) %in% row.names(meta) | row.names(t(nrmcts)) %in% c("sample", "condition", "type") )]
meta_nrmcts$AVERAGE <- c(rowMeans(meta_nrmcts[ , -which(row.names(meta_nrmcts) %in% c("condition","type"))]))

ggplot(meta_nrmcts,aes(condition, AVERAGE)) +
  geom_smooth(method='lm') +
  theme_bw() +   theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Metabolism-related expression")

#metabolism and sensory together
ggplot() +
  geom_smooth(data=sensory_nrmcts, aes(x=condition, y=AVERAGE), color= "forestgreen", fill = "green", method='lm') +
  geom_smooth(data=meta_nrmcts, aes(x=condition, y=AVERAGE), color= "blue", fill = "dodgerblue", method='lm') +
  theme_bw() +
  theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Metabolism vs. sensory-related expression")

#size and non-size together
ggplot() +
  geom_smooth(data=size_nrmcts, aes(x=condition, y=AVERAGE), color= "blue", fill = "dodgerblue", method='lm') +
  geom_smooth(data=nonsize_nrmcts, aes(x=condition, y=AVERAGE), color= "forestgreen", fill = "green", method='lm') +
  theme_bw() +
  theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Size-linked vs. non-size-linked expression")

#growth and down together
ggplot() +
  geom_smooth(data=growth_nrmcts, aes(x=condition, y=AVERAGE), color= "blue", fill = "dodgerblue", method='lm') +
  geom_smooth(data=down_nrmcts, aes(x=condition, y=AVERAGE), color= "darkred", fill = "coral", method='lm') +
  theme_bw() +
  theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Growth-correlated vs. growth negatively-correlated expression")

#green yellow module
gy <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/no_low_25merged_greenyellow.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
gy_nrmcts <- nrmcts[which(row.names(t(nrmcts)) %in% gy$x | row.names(t(nrmcts)) %in% c("sample", "condition", "type") )]
gy_nrmcts$AVERAGE <- c(rowMeans(gy_nrmcts[ , -which(names(gy_nrmcts) %in% c("condition","type"))]))

ggplot(gy_nrmcts,aes(condition, AVERAGE)) +
  geom_smooth(method='lm', color= "green", fill = "greenyellow") +
  theme_bw() +   theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Green yellow module expression")

#PC1 PC2 top loading genes
PC1 <- read.csv('C:/Users/imura/Documents/grad_5/RNA_paper/top_500_PC1_loadings.csv', sep = ",", header = TRUE, stringsAsFactors = FALSE)
PC2 <- read.csv('C:/Users/imura/Documents/grad_5/RNA_paper/top_500_PC2_loadings.csv', sep = ",", header = TRUE, stringsAsFactors = FALSE)

PC1_nrmcts <- nrmcts[which(row.names(t(nrmcts)) %in% PC1$gene | row.names(t(nrmcts)) %in% c("sample", "condition", "type") )]
PC1_nrmcts$AVERAGE <- c(rowMeans(PC1_nrmcts[ , -which(row.names(PC1_nrmcts) %in% c("condition","type"))]))

PC2_nrmcts <- nrmcts[which(row.names(t(nrmcts)) %in% PC2$gene | row.names(t(nrmcts)) %in% c("sample", "condition", "type") )]
PC2_nrmcts$AVERAGE <- c(rowMeans(PC2_nrmcts[ , -which(names(PC2_nrmcts) %in% c("condition","type"))]))

ggplot() +
  geom_smooth(data=PC1_nrmcts, aes(x=condition, y=AVERAGE), color= "blue", fill = "dodgerblue", method='lm') +
  geom_smooth(data=PC2_nrmcts, aes(x=condition, y=AVERAGE), color= "purple", fill = "violet", method='lm') +
  theme_bw() +
  theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("PC1-correlated vs. PC2-correlated expression")

#Media-linked expression and mushroom body proportional volume together
df <- read.csv("C:/Users/imura/Documents/old_csv_downloads/Eva_Atta.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

ggplot() +
  geom_smooth(data=df, aes(x=Bin, y=MBS), color= "blue", fill = "dodgerblue", method='loess') +
  geom_smooth(data=u_nrmcts, aes(x=condition, y=AVERAGE), color= "forestgreen", fill = "green", method='loess') +
  theme_bw() +
  theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Media-linked expression vs. mushroom body proportional volume")

#to instead combine in post
parOrig <- par()
par(bg=NA)
ggplot() +
  geom_smooth(data=df, aes(x=Bin, y=MBS), color= "green", fill = "greenyellow", method='loess') +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Media-linked expression vs. MB proportional volume")
par(parOrig)

ggplot() +
  geom_smooth(data=u_nrmcts, aes(x=condition, y=AVERAGE), color= "forestgreen", fill = "green", method='loess') +
  theme_bw() +
  theme(text = element_text(size=20)) +
  xlab("Worker size") + ylab("Media-linked expression vs. MB proportional volume")




