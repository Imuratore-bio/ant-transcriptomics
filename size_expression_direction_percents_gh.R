data <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/GSEA/amil_defog_iso2go4.csv")

i <- sapply(data, is.factor)
data[i] <- lapply(data[i], as.character)

majmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/no_low_padj_maj_min_no_NA.csv")
majmed <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/no_low_padj_maj_med_no_NA.csv")
medmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/no_low_padj_med_min_no_NA.csv")

maj_up = 0

med_up = 0

min_up = 0


#counts as upregulated if up relative to only one other group
for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] > 0 & majmin$padj[row] < 0.05) | (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] > 0 & majmed$padj[row] < 0.05)) {
    
    maj_up = maj_up + 1
    
  }
}

maj_up

for (row in 1:nrow(majmed)) {
  
  if  ((length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] < 0 & majmed$padj[row] < 0.05) | (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] > 0 & medmin$padj[row] < 0.05)) {
    
    med_up = med_up + 1
    
  }
}

med_up

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] < 0 & majmin$padj[row] < 0.05) | (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] < 0 & medmin$padj[row] < 0.05)) {
    
    min_up = min_up + 1
    
  }
}

min_up

####
olfactory_genes = 0
for (row in 1:nrow(data)) {
  
  if  (grepl("GO:0050911|GO:0042048|GO:0007608|GO:0021988|GO:0021989|GO:0021772|GO:0004984|GO:1990834|GO:0038022|GO:0004930", data$GO.term.accession[row])){
    
    olfactory_genes = olfactory_genes + 1
    
  }
}

olfactory_genes

olfactory_genes/nrow(data)

####combining both parts

olfact_majup = 0

olfact_medup = 0

olfact_minup = 0

for (row in 1:nrow(majmin)) {
  
  if  (((!is.na(majmin$log2FoldChange[row]) & majmin$log2FoldChange[row] > 2.5 & majmin$padj[row] < 0.05) | (!is.na(majmed$log2FoldChange[row]) & majmed$log2FoldChange[row] > 2.5 & majmed$padj[row] < 0.05)) & ( majmin$X[row] %in% data$ï..Transcript.stable.ID & !is.na(data$GO.term.accession[which(data$ï..Transcript.stable.ID == majmin$X[row])]) & grepl("GO:0050911|GO:0042048|GO:0007608|GO:0021988|GO:0021989|GO:0021772|GO:0004984|GO:1990834|GO:0038022|GO:0004930", data$GO.term.accession[grep(majmin$X[row], data$ï..Transcript.stable.ID)]))) {
    
    olfact_majup = olfact_majup + 1
    
  }
}

olfact_majup/olfactory_genes

for (row in 1:nrow(majmed)) {
  
  if  (((!is.na(majmed$log2FoldChange[row]) & majmed$log2FoldChange[row] < -2.5 & majmed$padj[row] < 0.05) | (!is.na(medmin$log2FoldChange[row])) & medmin$log2FoldChange[row] > 2.5 & medmin$padj[row] < 0.05) & ( majmed$X[row] %in% data$ï..Transcript.stable.ID & !is.na(data$GO.term.accession[which(data$ï..Transcript.stable.ID == majmed$X[row])]) & grepl("GO:0050911|GO:0042048|GO:0007608|GO:0021988|GO:0021989|GO:0021772|GO:0004984|GO:1990834|GO:0038022|GO:0004930", data$GO.term.accession[grep(majmed$X[row], data$ï..Transcript.stable.ID)]))) {
    
    olfact_medup = olfact_medup + 1
    
  }
}

olfact_medup/olfactory_genes

for (row in 1:nrow(majmin)) {
  
  if  (((!is.na(majmin$log2FoldChange[row]) & majmin$log2FoldChange[row] < -2.5 & majmin$padj[row] < 0.05) | (!is.na(medmin$log2FoldChange[row]) & medmin$log2FoldChange[row] < -2.5 & medmin$padj[row] < 0.05)) & ( majmin$X[row] %in% data$ï..Transcript.stable.ID & !is.na(data$GO.term.accession[which(data$ï..Transcript.stable.ID == majmin$X[row])]) & grepl("GO:0050911|GO:0042048|GO:0007608|GO:0021988|GO:0021989|GO:0021772|GO:0004984|GO:1990834|GO:0038022|GO:0004930", data$GO.term.accession[grep(majmin$X[row], data$ï..Transcript.stable.ID)]))) {
    
    olfact_minup = olfact_minup + 1
    
  }
}

olfact_minup/olfactory_genes

#isolate DEGs not linked to growth 
data <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/GSEA/amil_defog_iso2go4.csv", stringsAsFactors = FALSE)

i <- sapply(data, is.factor)
data[i] <- lapply(data[i], as.character)

majmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/padj_maj_min_no_NA.csv", stringsAsFactors = FALSE)
majmed <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/padj_maj_med_no_NA.csv", stringsAsFactors = FALSE)
medmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/padj_med_min_no_NA.csv", stringsAsFactors = FALSE)

growth = 0

#using significance
for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] > 0 & majmin$padj[row] < 0.05) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] > 0 & majmed$padj[row] < 0.05) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] > 0 & medmin$padj[row] < 0.05)) {
    
    growth = growth + 1
    
  }
}

growth

#not using significance
growth = 0

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] > 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] > 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] > 0)) {
    
    growth = growth + 1
    
  }
}

growth

#to get u-shaped pattern genes
ushape = 0

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] < 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] > 0)) {
    
    ushape = ushape + 1
    
  }
}

ushape

#downward trend genes
down = 0

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] < 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] < 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] < 0)) {
    
    down = down + 1
    
  }
}

down

#exclude growth genes
majmin_nogrowth = majmin
majmed_nogrowth = majmed
medmin_nogrowth = medmin

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] > 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] > 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] > 0)) {
    
    majmin_nogrowth <- majmin_nogrowth[which(majmin_nogrowth$X != majmin$X[row]),] 

  }
}

write.csv(majmin_nogrowth,'C:/Users/imura/Documents/grad_4/seq_analysis/no_low_maj_min_padj_growth_excluded.csv')

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] > 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] > 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] > 0)) {
    
    majmed_nogrowth <- majmed_nogrowth[which(majmed_nogrowth$X != majmed$X[row]),] 
    
  }
}

write.csv(majmed_nogrowth,'C:/Users/imura/Documents/grad_4/seq_analysis/no_low_maj_med_padj_growth_excluded.csv')

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] > 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] > 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] > 0)) {
    
    medmin_nogrowth <- medmin_nogrowth[which(medmin_nogrowth$X != medmin$X[row]),] 
    
  }
}

write.csv(medmin_nogrowth,'C:/Users/imura/Documents/grad_4/seq_analysis/no_low_med_min_padj_growth_excluded.csv')

#exclude genes positively OR negatively correlated to worker size
majmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/no_low_padj_maj_min_no_NA.csv")
majmed <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/no_low_padj_maj_med_no_NA.csv")
medmin <- read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/no_low_padj_med_min_no_NA.csv")

majmin_nogrowth = majmin
majmed_nogrowth = majmed
medmin_nogrowth = medmin

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] > 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] > 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] > 0)) {
    
    majmin_nogrowth <- majmin_nogrowth[which(majmin_nogrowth$X != majmin$X[row]),] 
    
  }
}

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] < 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] < 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] < 0)) {
    
    majmin_nogrowth <- majmin_nogrowth[which(majmin_nogrowth$X != majmin$X[row]),]
    
  }
}

write.csv(majmin_nogrowth,'C:/Users/imura/Documents/grad_4/seq_analysis/no_low_maj_min_padj_size_linked_excluded.csv')

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] > 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] > 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] > 0)) {
    
    majmed_nogrowth <- majmed_nogrowth[which(majmed_nogrowth$X != majmed$X[row]),] 
    
  }
}

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] < 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] < 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] < 0)) {
    
    majmed_nogrowth <- majmed_nogrowth[which(majmed_nogrowth$X != majmed$X[row]),] 
    
  }
}

write.csv(majmed_nogrowth,'C:/Users/imura/Documents/grad_4/seq_analysis/no_low_maj_med_padj_size_linked_excluded.csv')

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] > 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] > 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] > 0)) {
    
    medmin_nogrowth <- medmin_nogrowth[which(medmin_nogrowth$X != medmin$X[row]),] 
    
  }
}

for (row in 1:nrow(majmin)) {
  
  if  ((length(na.omit(majmin$log2FoldChange[row])) > 0 & majmin$log2FoldChange[row] < 0) & (length(na.omit(majmed$log2FoldChange[row])) > 0 & majmed$log2FoldChange[row] < 0) & (length(na.omit(medmin$log2FoldChange[row])) > 0 & medmin$log2FoldChange[row] < 0)) {
    
    medmin_nogrowth <- medmin_nogrowth[which(medmin_nogrowth$X != medmin$X[row]),]
    
  }
}

write.csv(medmin_nogrowth,'C:/Users/imura/Documents/grad_4/seq_analysis/no_low_med_min_padj_size_linked_excluded.csv')