##adapted from: https://github.com/z0on/GO_MWU
modulemembership = read.csv("C:/Users/imura/Documents/grad_4/seq_analysis/GSEA/no_low_25merged_gene_module_correlation.csv", check.names=FALSE)
modulemembership[which(modulemembership$ME.greenyellow >= 0.8), "ME.greenyellow"] <- 0
modulemembership[which(abs(modulemembership$ME.greenyellow) >= 0.8), "ME.greenyellow"] <- 0

library(plyr)
library(dplyr)

setwd("/Users/imura/Documents/grad_4/seq_analysis/GSEA")

# Edit these to match your data file names: 
input="no_low_signlogp_med_min.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="amil_defog_iso2go4.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go-basic.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")


# Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Strawberry/perl/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.5,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=30,   # a GO category should contain at least this many genes to be considered
           #higher gives more overlap
           clusterCutHeight=0.75, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=0.001,
                  level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.01, # FDR cutoff to print in regular (not italic) font.
                  level3=0.001, # FDR cutoff to print in large bold font.
                  txtsize=2.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.6, # height of the hierarchical clustering tree
                  colors=c("blue4","red4", "dodgerblue2","firebrick1") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results

