# Pseudotime Analysis using Monocle 3

# Load libraries
library(Seurat)
library(monocle3)
library(ggplot2)

# Load labeled Seurat object
hb <- readRDS("results/hb_combined_seurat_labeled.rds")

# Convert to Monocle CellDataSet
cds <- as.cell_data_set(hb)

# Carry over cluster identities
cds@clusters$UMAP$clusters <- hb$seurat_clusters
colData(cds)$cell_type <- hb$cell_class_concise

# Use UMAP as the dimensional reduction
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)

# Select a root population (e.g., Hepatoblasts as in the paper)
# Automatically pick root cells from the cluster labeled "Hepatoblasts"
root_cells <- colnames(hb)[which(hb$cell_class_concise == "Hepatoblasts")]

# Order cells in pseudotime
cds <- order_cells(cds, root_cells = root_cells)

# Plot pseudotime trajectory with cell labels
p1 <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
) + ggtitle("Pseudotime Trajectory")
ggsave("results/monocle_pseudotime_umap.png", p1, width = 9, height = 7)

# Plot by manual cell types
p2 <- plot_cells(
  cds,
  color_cells_by = "cell_type",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
) + ggtitle("Trajectory Colored by Cell Type")
ggsave("results/monocle_celltype_umap.png", p2, width = 9, height = 7)

# Plot pseudotime along trajectory per cell type group
p3 <- plot_cells(
  cds,
  color_cells_by = "cell_type",
  show_trajectory_graph = FALSE
) + ggtitle("UMAP by Cell Type (No Trajectory Overlay)")
ggsave("results/monocle_celltype_static.png", p3, width = 9, height = 7)