library(qvalue)

#p values from deseq2
p_vals <- read.csv('C:/Users/imura/Documents/grad_4/seq_analysis/padj_all_redundant.csv')

#use unadjusted p vals, select column
pvalues <- p_vals$pvalue

qobj <- qvalue(p = pvalues)

qobj$pi0