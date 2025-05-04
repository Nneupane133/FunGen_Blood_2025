# Load necessary packages
if (!require("VennDiagram")) install.packages("VennDiagram", dependencies=TRUE)
if (!require("ggVennDiagram")) install.packages("ggVennDiagram", dependencies=TRUE)

library(VennDiagram)
library(ggVennDiagram)
library(grid)

# Read DEG CSVs (change 'gene' to your actual column name if different)
tr_df <- read.csv("TrimmomaticDEGsortedfile.csv")
bb_df <- read.csv("BBdukDEGsortedfile.csv")

# Extract gene names (assuming column is called 'gene')
# If your column name is different, adjust here
tr_genes <- unique(as.character(tr_df$gene))
bb_genes <- unique(as.character(bb_df$gene))

# Base Venn Diagram
venn.plot <- venn.diagram(
  x = list(Trimmomatic = tr_genes, BBduk = bb_genes),
  filename = NULL,
  fill = c("cornflowerblue", "orangered"),
  alpha = 0.5,
  cex = 2.5,                  # Number size
  fontface = "bold",
  cat.cex = 2,                # Label size
  cat.fontface = "bold",
  cat.col = c("cornflowerblue", "orangered"),
  cat.pos = c(-20, 20),       # Positioning text labels
  cat.dist = 0.05,
  margin = 0.1,
  main = "Overlap of DEGs from Two Pipelines",
  main.cex = 2.2,
  main.fontface = "bold"
)
grid.newpage()
grid.draw(venn.plot)


# ggVennDiagram (optional prettier version)
gene_lists <- list(Trimmomatic = tr_genes, BBduk = bb_genes)
ggVennDiagram(gene_lists, label_alpha = 0.1) +
  ggplot2::labs(title = "DEG Overlap") +
  ggplot2::theme_minimal()

# Extract overlap and unique sets
common_genes <- intersect(tr_genes, bb_genes)
only_trimmomatic <- setdiff(tr_genes, bb_genes)
only_bbduk <- setdiff(bb_genes, tr_genes)

# Save results
writeLines(common_genes, "Common_Genes.txt")
writeLines(only_trimmomatic, "Only_Trimmomatic.txt")
writeLines(only_bbduk, "Only_BBduk.txt")

# Summary output
cat("✅ Venn Diagram & Lists Generated!\n")
cat("Only Trimmomatic:", length(only_trimmomatic), "\n")
cat("Only BBduk:", length(only_bbduk), "\n")
cat("Common:", length(common_genes), "\n")

