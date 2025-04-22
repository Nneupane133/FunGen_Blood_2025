# FunGen_Blood_2025

## Introduction
This is the team project analyzing the RNAseq dataset for blood tissue from dogs. In this project, we used the RNAseq data to compare the body size difference between the big and small dog. 

## Hypothesis and Research question
We hypothesized that blood-based gene expression profiles would reveal key regulatory pathways associated with canine body size. We predicted that if transcriptomic variation in blood is linked to physiological traits such as growth, energy metabolism, and cellular maintenance, then we would expect to observe differential expression of genes involved in pathways like insulin signaling, oxidative metabolism, ECM remodeling, and stress response between large and small dog breeds.

## Goal
To investigate differential gene expression associated with canine body size by analyzing RNA-seq data from blood samples collected from 25 healthy individuals representing 14 distinct dog breeds (12 large, 13 small). we aim to:

- Identify differentially expressed genes (DEGs) between big and small dogs.
- Determine whether metabolic, growth, and ECM-related pathways are enriched among these DEGs.
- Explore the function of the insulin signaling pathway and related metabolic networks in relation to bodu size difference

## Methodology
1. **Sample Selection**
- Paired-end RNA-seq data from 25 healthy dogs across 14 breeds.
- Balanced distribution: 12 large dogs and 13 small dogs.
- Public datasets retrieved from NCBI SRA under BioProject IDs:
    - PRJNA629466, PRJNA803741, PRJNA823683, PRJNA1011250, PRJNA1028815.
- Inclusion criteria: known age, sex, disease-free, and sequencing method consistency.

2. **Data Download & Preprocessing**
- Used SRA Toolkit (fastq-dump) to download FASTQ files.
- Quality check: FastQC v0.10.1.
- Trimming: Trimmomatic v0.39 to remove adapters and low-quality bases.

3. **Read Mapping**
- Aligned reads to Canis lupus familiaris reference genome (Dog_Tasha_GCF_000002285.5) using HISAT2 v2.2.0.

4. **Quantification**
- Used SAMtools and StringTie v2.2.1 for gene-level quantification and transcript assembly.
- Gene count matrix generated with prepDE.py3.

5. **Statistical Modeling (DESeq2)**: (for Differential Gene Expression Analysis (DGE))
- Software: R 4.4.2 and DESeq2 v1.46.0.
- Filtered out genes with < 20 total reads.
- Significance threshold: adjusted p-value < 0.05.
- Identified 495 significantly differentially expressed genes (DEGs).

6. **Data Visualization**
- MA plot: visualized gene fold changes vs mean expression.
Heatmap: top 20 variable genes (via pheatmap).
PCA plot: used to assess clustering and detect outliers (ggplot2, ggrepel).
Outliers: SRR18645116, SRR18645117, SRR18645125 (due to library prep method).

7. **GSEA Preparation**
- DEGs ranked using: rank = sign(log2FC) * -log10(p-value).
- Ranked list saved as .rnk file.
- Used full gene set to detect enriched pathways.
- Calculated Enrichment Score (ES) and Normalized Enrichment Score (NES).

8. **Protein Variant Analysis**
- Targeted Gene Mapping:
    - Focused on insulin/IGF pathway genes: INS, IGF1, IGF2, INSR, IGF1R, IGFBPs, IRS1–4.
    - Reads aligned with HISAT2, counted with StringTie, converted using SAMtools.

9. **Variant Calling**
- Performed using bcftools.
- Protein variant analysis conducted with Geneious Prime.

## Key findings
- 858 DEGs were identified at p < 0.1, with 495 genes significant at p < 0.05.
- Pathways upregulated in big dogs included:
    - ECM-receptor interaction
    - Focal adhesion 
- Pathways downregulated in big dogs included:
    - Oxidative phosphorylation
    - TCA cycle
    - Ribosome
    - Glycolysis/gluconeogenesis
- Insulin signaling pathway was not significantly enriched, but several related genes were explored for expression and variant analysis.

-  Variant analysis of IIS genes (e.g., INS, IGF1, IGF2R, INSR) identified SNP-level changes potentially related to breed-specific physiological differences.
