# scripts/marker_gene_analysis.R

library(Seurat)
library(dplyr)
library(ggplot2)

# Load object
hb_combined <- readRDS("results/hb_combined_seurat.rds")
DefaultAssay(hb_combined) <- "RNA"

# 🔧 Join layers (fixes the error!)
hb_combined <- JoinLayers(hb_combined)

# Run DE analysis
markers <- FindAllMarkers(hb_combined,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)

# Make sure markers were found
if (nrow(markers) == 0) {
  stop("❌ No differentially expressed genes were identified.")
}

# Save full marker table
write.csv(markers, file = "results/all_markers_table.csv", row.names = FALSE)

# Get top 5 markers per cluster
top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

write.csv(top5, file = "results/top5_markers_per_cluster.csv", row.names = FALSE)

# Heatmap of top 5 markers per cluster
top_genes <- unique(top5$gene)

pdf("results/heatmap_top5_genes.pdf", width = 10, height = 8)
DoHeatmap(hb_combined, features = top_genes) +
  ggtitle("Top 5 Marker Genes per Cluster")
dev.off()

cat("✅ Marker gene analysis complete.\n")