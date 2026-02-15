# ============================================================
# 1. Load Libraries
# ============================================================

# Install only if not already installed
if (!require(pheatmap)) install.packages("pheatmap", dependencies = TRUE)

library(pheatmap)

# ============================================================
# 2. Define DEG Analysis Function
# ============================================================

analyze_file <- function(filepath, label) {
  
  cat("\n====================================\n")
  cat("Analyzing:", label, "\n")
  cat("====================================\n")
  
  # Load data
  data <- read.csv(filepath)
  
  # Remove rows where adjusted p-value is NA
  data <- data[!is.na(data$padj), ]
  
  # Define DEGs
  upregulated <- subset(data, padj < 0.05 & log2FoldChange > 1)
  downregulated <- subset(data, padj < 0.05 & log2FoldChange < -1)
  
  # Print summary
  cat("Total genes analyzed:", nrow(data), "\n")
  cat("Upregulated genes:", nrow(upregulated), "\n")
  cat("Downregulated genes:", nrow(downregulated), "\n")
  
  # Save DEG tables
  write.csv(upregulated,
            paste0("upregulated_", gsub(" ", "_", label), ".csv"),
            row.names = FALSE)
  
  write.csv(downregulated,
            paste0("downregulated_", gsub(" ", "_", label), ".csv"),
            row.names = FALSE)
  
  # Create plots
  create_volcano(data, label)
  create_heatmap(data, label)
}

# ============================================================
# 3. Volcano Plot Function
# ============================================================

create_volcano <- function(data, label) {
  
  data$Significance <- "Not Significant"
  data$Significance[data$padj < 0.05 & data$log2FoldChange > 1] <- "Upregulated"
  data$Significance[data$padj < 0.05 & data$log2FoldChange < -1] <- "Downregulated"
  
  png(paste0("volcano_", gsub(" ", "_", label), ".png"),
      width = 1000, height = 800)
  
  plot(data$log2FoldChange,
       -log10(data$padj),
       pch = 20,
       col = ifelse(data$Significance == "Upregulated", "#D62728",
                    ifelse(data$Significance == "Downregulated", "#1F77B4", "grey80")),
       xlim = c(-2, 2),
       ylim = c(0, 10),
       main = paste("Volcano Plot:", label),
       xlab = "log2 Fold Change",
       ylab = "-log10 adjusted p-value",
       cex.main = 1.4,
       cex.lab = 1.2)
  
  abline(v = c(-1, 1), lty = 2)
  abline(h = -log10(0.05), lty = 2)
  
  legend("topright",
         legend = c("Upregulated", "Downregulated"),
         col = c("#D62728", "#1F77B4"),
         pch = 20,
         bty = "n")
  
  dev.off()
}

# ============================================================
# 4. Heatmap Function
# ============================================================

create_heatmap <- function(data, label) {
  
  sig <- subset(data, padj < 0.05 & abs(log2FoldChange) > 1)
  
  if (nrow(sig) < 5) {
    cat("Not enough DEGs for heatmap:", label, "\n")
    return()
  }
  
  sig <- sig[order(sig$padj), ]
  top_genes <- head(sig, 30)
  
  matrix_data <- as.matrix(top_genes$log2FoldChange)
  rownames(matrix_data) <- top_genes[,1]
  
  # Scale for visualization
  matrix_data <- scale(matrix_data)
  
  png(paste0("heatmap_", gsub(" ", "_", label), ".png"),
      width = 800, height = 1000)
  
  pheatmap(matrix_data,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           border_color = NA,
           fontsize_row = 8,
           main = paste("Top 30 DEGs:", label))
  
  dev.off()
}

install.packages("VennDiagram")
library(VennDiagram)

create_venn <- function(file1, file2, file3) {
  
  d1 <- read.csv(file1)
  d2 <- read.csv(file2)
  d3 <- read.csv(file3)
  
  d1 <- subset(d1, padj < 0.05 & abs(log2FoldChange) > 1)
  d2 <- subset(d2, padj < 0.05 & abs(log2FoldChange) > 1)
  d3 <- subset(d3, padj < 0.05 & abs(log2FoldChange) > 1)
  
  genes1 <- d1[,1]
  genes2 <- d2[,1]
  genes3 <- d3[,1]
  
  venn.plot <- venn.diagram(
    x = list(
      "Group1 vs Group2" = genes1,
      "Group2 vs Group3" = genes2,
      "Group1 vs Group3" = genes3
    ),
    filename = "venn_DEG_overlap.png",
    fill = c("#FF9999", "#99CCFF", "#99FF99"),
    alpha = 0.6,
    cex = 1.5,
    cat.cex = 1.2,
    cat.pos = 0
  )
}
create_venn(
  "deseq2_kallisto_liver.group1_vs_2_minus_outliers.csv",
  "deseq2_kallisto_liver.group2_vs_3_minus_outliers.csv",
  "deseq2_kallisto_liver.group1_vs_3_minus_outliers.csv"
)


# ============================================================
# 5. Run Analysis
# ============================================================

analyze_file("deseq2_kallisto_liver.group1_vs_2_minus_outliers.csv",
             "Group1 vs Group2")

analyze_file("deseq2_kallisto_liver.group2_vs_3_minus_outliers.csv",
             "Group2 vs Group3")

analyze_file("deseq2_kallisto_liver.group1_vs_3_minus_outliers.csv",
             "Group1 vs Group3")

