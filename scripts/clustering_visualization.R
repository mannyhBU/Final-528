# clustering_visualization.R

library(Seurat)
library(ggplot2)

# Load combined Seurat object
hb_combined <- readRDS("results/hb_combined_seurat.rds")

# Run clustering
hb_combined <- FindNeighbors(hb_combined, dims = 1:30, reduction = "harmony")
hb_combined <- FindClusters(hb_combined, resolution = 0.5)

# Run UMAP (if not already done)
hb_combined <- RunUMAP(hb_combined, dims = 1:30, reduction = "harmony")

# Plot UMAP colored by clusters
umap_clusters <- DimPlot(hb_combined, reduction = "umap", group.by = "seurat_clusters", 
                         label = TRUE, pt.size = 0.4) + 
  ggtitle("UMAP: Cell Clusters")
ggsave("results/umap_clusters.png", umap_clusters, width = 8, height = 6)

# Plot UMAP colored by sample
umap_samples <- DimPlot(hb_combined, reduction = "umap", group.by = "sample", 
                        label = TRUE, pt.size = 0.4) + 
  ggtitle("UMAP: Sample Origin")
ggsave("results/umap_samples.png", umap_samples, width = 8, height = 6)

# Save updated Seurat object
saveRDS(hb_combined, file = "results/hb_combined_seurat_clustered.rds")

cat("\u2705 Clustering and visualization complete. UMAPs and updated Seurat object saved.\n")
