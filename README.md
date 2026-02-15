# RNAseq Liver DEG Analysis

This repository contains an R script (`analysis.R`) for performing Differential Expression Analysis (DEA) on RNA-seq data. It automates the process of filtering differentially expressed genes (DEGs), generating statistical plots, and analyzing overlaps between comparison groups.

## Overview

The analysis pipeline performs the following steps:

1. **Data Loading**: Reads DESeq2 output files (CSV format).
2. **Filtering**: Identifies significant DEGs using standard thresholds:
   - Adjusted p-value < 0.05
   - |log2 Fold Change| > 1
3. **File Generation**: Exports lists of upregulated and downregulated genes to CSV.
4. **Visualization**:
   - **Volcano Plots**: Visualizes significance (-log10 padj) vs. fold change (log2FC).
   - **Heatmaps**: Displays the expression profiles of the top 30 significant genes.
   - **Venn Diagram**: Shows the overlap of DEGs across three specific comparison groups.

## Dependencies

The script requires the following R packages:

- `pheatmap` (for heatmap generation)
- `VennDiagram` (for overlap analysis)

The script checks for `pheatmap` and attempts to install it if missing. You may need to install `VennDiagram` manually if it is not present.

```r
install.packages("pheatmap")
install.packages("VennDiagram")
```

## Input Data

The script expects CSV files containing differential expression results. The required columns are:

- `padj`: Adjusted p-value
- `log2FoldChange`: Log2 Fold Change
- *(Row Names)*: Gene identifiers

The default configuration processes the following input files located in the working directory:

- `deseq2_kallisto_liver.group1_vs_2_minus_outliers.csv`
- `deseq2_kallisto_liver.group2_vs_3_minus_outliers.csv`
- `deseq2_kallisto_liver.group1_vs_3_minus_outliers.csv`

## Usage

To run the analysis, execute the script in your R environment:

```bash
Rscript analysis.R
```

Ensure your input CSV files are in the same directory as the script.

## Output Files

For each analyzed comparison (e.g., "Group1 vs Group2"), the script generates:

### Data Tables

- `upregulated_[Label].csv`: Genes with padj < 0.05 and log2FC > 1.
- `downregulated_[Label].csv`: Genes with padj < 0.05 and log2FC < -1.

### Plots

- **Volcano Plot** (`volcano_[Label].png`):
  - Red: Upregulated
  - Blue: Downregulated
  - Grey: Not significant
- **Heatmap** (`heatmap_[Label].png`): Clustered heatmap of the top 30 significant DEGs.

### Summary

- `venn_DEG_overlap.png`: A Venn diagram illustrating the intersection of significant genes between the three comparison groups.
