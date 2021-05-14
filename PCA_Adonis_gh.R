normcounts <- read.csv(file="C:/Users/imura/Documents/grad_4/seq_analysis/adonis4.csv", row.names=1)
head(normcounts)

treat=c("pointsix", "pointsix", "onepointeight", "onepointeight", "three", "pointsix", "onepointeight", "onepointeight", "three", "pointsix", "onepointeight", "three", "pointsix", "pointsix", "onepointeight", "onepointeight", "pointsix", "pointsix", "three", "three", "three", "pointsix", "onepointeight", "three", "onepointeight", "pointsix", "three")
colony=as.factor(c("Ac22", "Ac22", "Ac22", "Ac22", "Ac22", "M1", "M1", "M1", "M1", "M2", "M2", "M2", "Ac22", "M1", "M1", "M1", "M2", "M2", "M2", "M2", "M2", "Ac22", "Ac22", "Ac22", "M1", "M2", "M2"))
g=data.frame(treat, colony)
head(g)
colData<- g
normcounts_t=t(normcounts)
which(apply(normcounts_t, 2, var)==0) #identify which columns have zero variance, PCA will not work on these
normcounts_t = normcounts_t[ , apply(normcounts_t, 2, var) != 0] #removes columns with zero variance

pca <- prcomp(normcounts_t,center = TRUE, scale. = TRUE)
#head(pca)

#Introduce principle components
li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
#grabs first two/top principal components
pca_s <- pca_s[,c(1,2)]
pca_s$Samples = row.names(pca_s)
pca_s$treat=colData$treat
head(pca_s)

#biplot to show eigenvalues
#https://www.datacamp.com/community/tutorials/pca-analysis-r
ggbiplot(prcomp(normcounts_t[,1:20]),ellipse=TRUE, labels=rownames(normcounts_t[,1:20]), groups=treat)

#Run adonis for significance
library(vegan)
adonis(pca$x ~ treat, data = pca_s, method='man', na.rm = TRUE)