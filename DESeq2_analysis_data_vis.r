####  DESeq2 Script for Differential Gene Expression Analysis in 
      # Functional Genomics BIOL: 6850
### Resources and Citations:
# Love et al 2016 DESeq2 GenomeBiology
# https://bioconductor.riken.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

### set working directory
setwd("~/Desktop/FunGen/Final project")

#### Install the DESeq2 package if you have not already
library(DESeq2)

## Input data ##

### Input the count data, the gene(/transcript) count matrix and labels
  ### How you inport this will depend on what your final output was from the mapper/counter that you used.
  ## this works with output from PrepDE.py from Ballgown folder.
countdata <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
dim(countdata)
head(countdata)

#remove gene- and  everything after | symbol to clean up
rownames(countdata) <- gsub("^gene-", "", rownames(countdata))  # Remove "gene-" prefix
rownames(countdata) <- gsub("\\|.*", "", rownames(countdata))  # Remove "|" and everything after

#make sure it worked
head(rownames(countdata))


### Input the phenotype data
# Note: The 4_inputDESeq2_Final file contains information on each sample, e.g., size and SRA# The exact way to import this depends on the format of the file.
##  Make sure the individual names match between the count data and the metadata in the correct order
coldata <-(read.table("4_inputDEseq2_Final.txt", header=TRUE, row.names=1))
dim(coldata)
head(coldata)

#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(coldata) %in% colnames(countdata))
countdata <- countdata[, rownames(coldata)]
all(rownames(coldata) == colnames(countdata))
#should output TRUE for both if matched

## Create the DESEQ dataset and define the statistical model (page 6 of the manual)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData=coldata,  design = ~treatment)
dds
#40478 total genes

#####   Prefiltering    Manual - starting at  1.3.6 
# Here we perform a minimal pre-filtering to remove rows that have less than 20 reads mapped.
dds <- dds[ rowSums(counts(dds)) > 20, ]
dds
#24408 total genes so lost ~15000 genes after filtering

## set factors for statistical analyses
###### Note you need to change condition to treatment (to match our design above)
#  and levels to our treatment names in the 4_inputDESeq2_Final: Big and small size
# example:
dds$condition <- factor(dds$treatment, levels=c("Big","Small"))

######   Differential expression analysis - running DESeq2
dds <- DESeq(dds)
res <- results(dds)
res

# order our results table by the smallest adjusted p value:
  resOrdered <- res[order(res$padj),]
  resOrdered
# summary of basic information, default is p<0.1
  summary(res)
#How many adjusted p-values were less than 0.1?
  sum(res$padj < 0.1, na.rm=TRUE)
#If the adjusted p value will be a value other than 0.1, alpha should be set to that value:
  res05 <- results(dds, alpha=0.05)
  summary(res05)
  sum(res05$padj < 0.05, na.rm=TRUE)
  
###  MA-plot
  ##plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts. 
  ## Points which fall out of the window are plotted as open triangles pointing either up or down
  plotMA(res, main="DESeq2 Results", ylim=c(-8,8), xlab="Mean of Normalized Counts", ylab="Log2 Fold Change")
  
  # Add highlighting for even more significant genes (if highlighted in red padj < 0.05, if in blue originally p-val < 0.1)
  with(res, points(res$baseMean[res$padj < 0.05], res$log2FoldChange[res$padj < 0.05], 
                   col="red", pch=16, cex=0.5))
  #dev.off()
  
  ##  Write your results to a file 
  write.csv(as.data.frame(resOrdered), file="DGESeq_results.csv")  
  
  ## Extracting transformed values
  rld <- rlog(dds)
  vsd <- varianceStabilizingTransformation(dds)
  head(assay(rld), 3)
  
### Heatmap of the count matrix
  #library("genefilter")
  topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
  
  library("pheatmap")
  mat  <- assay(vsd)[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
  anno <- as.data.frame(colData(vsd)[, c("treatment", "type")])
  df <- as.data.frame(colData(dds)[,c("treatment","type")])
    pheatmap(mat, annotation_col = anno)
  
  # Heatmap of the sample-to-sample distances
  sampleDists <- dist(t(assay(rld)))
  library("RColorBrewer")
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld$treatment)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  
 # 2.2.3 Principal component plot of the samples
  plotPCA(rld, intgroup=c("treatment"))
  
  
## Find the sample names of the outliers for analysis ##
  # Extract and transpose expression matrix since PCA plotting before is using rld
  expr_mat <- t(assay(rld))  # transpose so samples are rows
  
  # Run PCA on the above expression matrix
  pca <- prcomp(expr_mat, scale. = TRUE)
  
  # Turn into data frame
  pca_df <- as.data.frame(pca$x)
  pca_df$sample <- rownames(pca_df)
  
  # Add group info (e.g., Big vs Small)
  pca_df$group <- colData(rld)$treatment
  
  library(ggplot2)
  library(ggrepel)
  
  # Define what you want outliers to be (here, anything > 1.8 SD from the mean is indicated as an outlier)
  outliers <- subset(pca_df, abs(PC1) > 1.8 * sd(PC1) | abs(PC2) > 1.8 * sd(PC2))
  
  #plot PCA and label outliers with sample IDs
  ggplot(pca_df, aes(x = PC1, y = PC2)) + #make a plot
    geom_point(aes(color = group)) +  #add the points to the plot, different colors for Big/small
    geom_text_repel(data = outliers, aes(label = sample), color = "black") + #add labels to outlier points by pulling sample ID from sample column
    theme_minimal()
  
  #print outlier IDs
  print(outliers$sample)
  
## Make ranked list for GSEA ###
  ### Merge 'gene names' with DGE results by Gene Model
  
  ## Import the DGE results file make sure the gene model name is 'gene_id' to match annotation file
  DGEresults <- read.csv("DGESeq_results.csv", stringsAsFactors = FALSE)
  summary(DGEresults)
  dim(DGEresults)
  
  ## Rename first column so it matches "gene_id" in annotation file
  names(DGEresults)[1]<- "gene_id" 
  
## example aa <- within(resOrdered, z <- x + y - 2)
  #creating rank column within DGEresults
DGE_Rank <-  within(DGEresults, rank <- sign(log2FoldChange) * -log10(pvalue))
DGE_Rank 

#subset the results so only Gene Name and rank
DGErank <- DGE_Rank[, c("gene_id","rank")]
DGErank

#subset the results so only Gene Name and rank
DGErank_withName <- na.omit(DGErank)
DGErank_withName
dim(DGErank_withName)

#export as .rnk file for GSEA analysis
write.table(as.data.frame(DGErank_withName), file="DGErankName_final.rnk", quote=FALSE, row.names=FALSE, sep = "\t")  


####  We also need the normalized expression DATA
nt <- normTransform(dds) # defaults to log2(x+1)
head(assay(nt))
# compare to original count data
head(countdata)

# make it a new dataframe
NormTransExp<-assay(nt)
summary(NormTransExp)
head(NormTransExp)

# Convert row names into a column
NormTransExp <- data.frame(gene_id = rownames(NormTransExp), NormTransExp, row.names = NULL)
# make sure it worked
head(NormTransExp)

#subset the results so only Genes with Name
NormTransExp_withName <- na.omit(NormTransExp)
head(NormTransExp_withName)
dim(NormTransExp_withName)

#export as txt file
write.table(as.data.frame(NormTransExp_withName), file="NormTransExp_withName.txt", quote=FALSE, row.names=FALSE, sep = "\t")  


