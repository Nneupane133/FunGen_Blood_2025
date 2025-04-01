####  DESeq2 Script for Differential Gene Expression Analysis in 
      # Functional Genomics BIOL: 6850
### Resources and Citations:
# Love et al 2016 DESeq2 GenomeBiology
# https://bioconductor.riken.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html


### You will need to set your working directory to the location you have your data.
# You can do this by using  the Session menu to set working directory To Source File Directory
setwd("~/Desktop/FunGen/Final project")

#### Install the DESeq2 package if you have not already
## try http:// if https:// URLs are not supported
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("DESeq2")

## Load the DESeq2 library 
library(DESeq2)

## Use the Session menu to set working directory To Source File Directory

##########   1.3 Input data   ##############

### Input the count data, the gene(/transcript) count matrix and labels
  ### How you inport this will depend on what your final output was from the mapper/counter that you used.
  ## this works with output from PrepDE.py from Ballgown folder.
countdata <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
dim(countdata)
head(countdata)

### Input the meta data or phenotype data
# Note: The PHENO_DATA file contains information on each sample, e.g., sex or population. The exact way to import this depends on the format of the file.
##  Make sure the individual names match between the count data and the metadata
coldata <-(read.table("4_inputDEseq2_Final.txt", header=TRUE, row.names=1))
dim(coldata)
head(coldata)


#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(coldata) %in% colnames(countdata))
countdata <- countdata[, rownames(coldata)]
all(rownames(coldata) == colnames(countdata))


## Create the DESEQ dataset and define the statistical model (page 6 of the manual)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData=coldata,  design = ~treatment)
#look at it
dds

#####   Prefiltering    Manual - starting at  1.3.6 
# Here we perform a minimal pre-filtering to remove rows that have less than 20 reads mapped.
## You can play around with this number to see how it affects your results!
dds <- dds[ rowSums(counts(dds)) > 20, ]
# look.  How many genes were filtered out?
dds

## set factors for statistical analyses
###### Note you need to change condition to treatment (to match our design above)
#  and levels to our treatment names in the PHENO_DATA: Ad_lib is the control, Caloric_Restriction is the treatment group
# example:
#dds$condition <- factor(dds$condition, levels=c("untreated","treated"))
#dds$condition <- factor(dds$treatment, levels=c("Ad_lib","Caloric_restriction"))

######     1.4 Differential expression analysis
### Question 2. Look at the manual - what is happening at this point?
dds <- DESeq(dds)
res <- results(dds)
res


###  Question 3. What does each column mean?
# We can order our results table by the smallest adjusted p value:
  resOrdered <- res[order(res$padj),]
  resOrdered
# We can summarize some basic tallies using the summary function the default is p<0.1.
  summary(res)
#How many adjusted p-values were less than 0.1?
  sum(res$padj < 0.1, na.rm=TRUE)
#If the adjusted p value will be a value other than 0.1, alpha should be set to that value:
  res05 <- results(dds, alpha=0.05)
  summary(res05)
  sum(res05$padj < 0.05, na.rm=TRUE)

  
###    1.5.1 MA-plot
  ##plotMA shows the log2 fold changes attributable to a given variable over the meanof normalized counts. 
  ## Points will be colored red if the adjusted p value is less than 0.1. 
  ## Points which fall out of the window are plotted as open triangles pointing either up or down
  plotMA(res, main="DESeq2", ylim=c(-8,8))
  
  #After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. 
  # One can then recover the gene identiers by saving the resulting indices:
  #idx <- identify(res$baseMean, res$log2FoldChange)
    # after selecting a gene. You need to press escape to move on
 # rownames(res)[idx]

    
##  1.5.2 Plot counts - sanity check!
  
  ##  Write your results to a file 
  write.csv(as.data.frame(resOrdered), file="DGESeq_results.csv")  
  
  ## 2.1.2 Extracting transformed values
  rld <- rlog(dds)
  vsd <- varianceStabilizingTransformation(dds)
  head(assay(rld), 3)
  
  
### Heatmap of the count matrix
  #library("genefilter")
  topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
  
  library("pheatmap")
  mat  <- assay(vsd)[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
  anno <- as.data.frame(colData(vsd)[, c("treatment", "type")])
  df <- as.data.frame(colData(dds)[,c("treatment","type")])
    pheatmap(mat, annotation_col = anno)
  
  #2.2.2 Heatmap of the sample-to-sample distances
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
  

############################# Make ranked list for GSEA ####################
  
  ### Merge 'gene names' with DGE results by Gene Model
  
  ## Import the DGE results file make sure the gene model name is 'gene_id' to match annotation file
  DGEresults <- read.csv("DGESeq_results.csv", stringsAsFactors = FALSE)
  summary(DGEresults)
  dim(DGEresults)
  
  ## Rename first column so it matches "gene_id" in annotation file
  names(DGEresults)[1]<- "gene_id" 
  
## example aa <- within(resOrdered, z <- x + y - 2)
DGE_Rank <-  within(DGEresults, rank <- sign(log2FoldChange) * -log10(pvalue))
DGE_Rank 

#had extra column of numbering genes so remove row names and set to default
#rownames(DGE_Rank) <- NULL
# make sure it worked
#head(DGE_Rank)

# Remove "gene-" from the gene_id column
DGE_Rank$gene_id <- gsub("^gene-", "", DGE_Rank$gene_id)
# Make sure it worked
head(DGE_Rank$gene_id)

# Remove everything after the pipe (including the pipe)
DGE_Rank$cleaned_gene_id <- gsub("\\|.*", "", DGE_Rank$gene_id)
# make sure it worked
head(DGE_Rank$cleaned_gene_id)

#combine gene and ranked list (DGE_Rank)
DGErankIDs  <-cbind(DGE_Rank$cleaned_gene_id,DGE_Rank)
head(DGErankIDs)
summary(DGErankIDs)  

#write.csv(as.data.frame(DGErankIDs), file="DGErankIDs.csv", row.names=FALSE)  
write.table(as.data.frame(DGErankIDs), file="DGErankName.rnk", quote=FALSE, row.names=FALSE, sep = "\t")  


####  We also need the normalized expression DATA
nt <- normTransform(dds) # defaults to log2(x+1)
head(assay(nt))
# compare to original count data
head(countdata)

# make it a new dataframe
NormTransExp<-assay(nt)
summary(NormTransExp)
head(NormTransExp)

# Remove "gene-" from the row names
rownames(NormTransExp) <- gsub("^gene-", "", rownames(NormTransExp))
NormTransExpIDs  <-cbind(cleaned_gene_id,NormTransExp)
head(NormTransExpIDs)

# Remove everything after the pipe symbol (|) in the gene ID
NormTransExp_cleaned <- gsub("\\|.*", "", rownames(NormTransExp))
# make sure it worked
head(NormTransExp_cleaned)

## Rename first column so it matches "gene_id" in annotation file
names(NormTransExp_cleaned)[1]<- "gene_id"
head(NormTransExp_cleaned)

#subset the results so only Genes with Name
NormTransExp_cleaned_withName <- na.omit(NormTransExp_cleaned)
head(NormTransExp_cleaned_withName)
dim(NormTransExp_cleaned_withName)
#output NULL because character vector

# Convert the character vector into a data frame
NormTransExp_cleaned_withName <- data.frame(gene_id = NormTransExp_cleaned_withName)

# Check the dimensions
dim(NormTransExp_cleaned_withName)
#should output 2 numbers NOT NULL

#write.csv(as.data.frame(NormTransExpIDs), file="NormTransExpressionData.csv", row.names=FALSE)  
write.table(as.data.frame(NormTransExp_cleaned_withName), file="NormTransExp_cleaned_Names.txt", quote=FALSE, row.names=FALSE, sep = "\t")  


