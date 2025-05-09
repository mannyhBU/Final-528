library(Seurat)
library(ggplot2)

# Load Seurat object
hb <- readRDS("results/hb_combined_seurat.rds")


# Load built-in cell cycle gene lists
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Score cells based on expression of S and G2/M markers
hb <- CellCycleScoring(hb, s.features = s.genes, g2m.features = g2m.genes)

# Save UMAP colored by phase
umap_phase <- DimPlot(hb, group.by = "Phase", pt.size = 0.5) +
  ggtitle("Cell Cycle Phase per Cell")
ggsave("results/umap_cell_cycle_phase.png", umap_phase, width = 8, height = 6)

# Save violin plot of scores
vln <- VlnPlot(hb, features = c("S.Score", "G2M.Score"), group.by = "seurat_clusters", pt.size = 0.1) +
  ggtitle("Cell Cycle Scores by Cluster")
ggsave("results/vln_cell_cycle_scores.png", vln, width = 8, height = 6)

# Save updated object
saveRDS(hb, file = "results/hb_combined_seurat_with_cellcycle.rds")
