# scripts/pca_analysis.R

library(Seurat)
library(ggplot2)

# Load processed Seurat object
hb_combined <- readRDS("results/hb_combined_seurat.rds")

# Run PCA on the highly variable genes
hb_combined <- RunPCA(hb_combined, features = VariableFeatures(object = hb_combined))

# Save PCA results to RDS
saveRDS(hb_combined, file = "results/hb_combined_seurat.rds")

# Plot PCA with groupings by sample
pca_plot <- DimPlot(hb_combined, reduction = "pca", group.by = "sample") +
  ggtitle("PCA Plot Grouped by Sample")
ggsave("results/pca_plot_by_sample.png", pca_plot, width = 8, height = 6)

# Elbow plot to determine how many PCs to keep
elbow_plot <- ElbowPlot(hb_combined, ndims = 50) +
  ggtitle("Elbow Plot of PCs")
ggsave("results/elbow_plot.png", elbow_plot, width = 8, height = 6)

# Summary output
cat("\n✅ PCA complete. PCA and Elbow plots saved in results/.\n")
