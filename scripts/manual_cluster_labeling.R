# Manual Cluster Labeling Script (Updated with Paper-Informed Labels)

library(Seurat)
library(dplyr)
library(ggplot2)

# Load your clustered object
hb <- readRDS("results/hb_combined_seurat.rds")

# STEP 1: Find all marker genes per cluster
markers <- FindAllMarkers(hb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(markers, file = "results/all_cluster_markers.rds")
write.csv(markers, file = "results/all_cluster_markers.csv", row.names = FALSE)

# STEP 2: Top 5 markers per cluster
top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# STEP 3: Heatmap of top marker genes
top_genes <- unique(top5$gene)
heatmap <- DoHeatmap(hb, features = top_genes) + ggtitle("Top 5 Marker Genes per Cluster")
ggsave("results/manual_top5_heatmap.png", heatmap, width = 10, height = 8)

# STEP 4: Identify top 3 largest clusters
cluster_sizes <- table(Idents(hb))
top_clusters <- names(sort(cluster_sizes, decreasing = TRUE)[1:3])

# Plot violin and feature plots for top 3 clusters
for (cluster_id in top_clusters) {
  cluster_genes <- top5 %>% filter(cluster == cluster_id) %>% pull(gene) %>% head(5)

  vln <- VlnPlot(hb, features = cluster_genes, pt.size = 0.1) + ggtitle(paste("Violin Plot - Cluster", cluster_id))
  ggsave(paste0("results/vln_cluster_", cluster_id, ".png"), vln, width = 10, height = 6)

  feat <- FeaturePlot(hb, features = cluster_genes, cols = c("lightgrey", "blue")) + ggtitle(paste("Feature Plot - Cluster", cluster_id))
  ggsave(paste0("results/feature_cluster_", cluster_id, ".png"), feat, width = 10, height = 6)
}

# STEP 5: Assign paper-informed manual labels to clusters
manual_labels <- c(
  "0" = "Tumor",
  "1" = "Hepatocyte",
  "2" = "Endothelial",
  "3" = "Immune",
  "4" = "Kupffer",
  "5" = "Stellate",
  "6" = "Cholangiocyte",
  "7" = "Tumor",
  "8" = "Hepatocyte",
  "9" = "Immune",
  "10" = "Endothelial",
  "11" = "Tumor",
  "12" = "Immune",
  "13" = "Tumor",
  "14" = "Stellate",
  "15" = "NK/T Cell",
  "16" = "Tumor",
  "17" = "B Cell",
  "18" = "Unknown"
)

hb$cell_class_concise <- manual_labels[as.character(hb$seurat_clusters)]
Idents(hb) <- "cell_class_concise"

# STEP 6: Plot with manual labels
umap_labeled <- DimPlot(hb, group.by = "cell_class_concise", label = TRUE, repel = TRUE) + 
  ggtitle("UMAP with Manual Labels")
ggsave("results/umap_manual_labels.png", umap_labeled, width = 10, height = 8)

# STEP 7: Save final labeled object
saveRDS(hb, file = "results/hb_combined_seurat_labeled.rds")
cat("\u2705 Manual labeling complete with literature-based annotations. All outputs saved in results/.\n")