library(Seurat)
library(ggplot2)
library(plyr)
library(harmony)

# Optional: remove or adjust if you're running from the correct working directory
# setwd("/projectnb/bf528/students/mannyh/final-project")

# Get the list of sample directories
sample_subdirs <- list.dirs("data/processed", recursive = FALSE)

count <- 1

for (folder in sample_subdirs) {
  # Look for the subfolder that contains the actual matrix files
  matrix_folder <- list.files(folder, pattern = "filtered_feature_bc_matrix$", full.names = TRUE)

  if (length(matrix_folder) == 0 || !dir.exists(matrix_folder)) {
    cat("⚠️ Skipping:", folder, "- no matrix folder found\n")
    next
  }

  sample_name <- basename(folder)
  cat("🔄 Processing:", sample_name, "\n")

  # Load data
  hb_data <- Read10X(data.dir = matrix_folder)
  hb <- CreateSeuratObject(counts = hb_data, min.cells = 0, min.features = 0)

  # Calculate percent mitochondrial
  hb[["percent.mt"]] <- PercentageFeatureSet(hb, pattern = "^MT-")

  # Filtering
  hb <- subset(hb, subset = nFeature_RNA > 500 & nCount_RNA > 800 & percent.mt < 10)

  # Rename cells and add metadata
  hb <- RenameCells(hb, add.cell.id = sample_name)
  hb$sample <- sample_name

  # Merge into combined object
  if (count == 1) {
    hb_combined <- hb
  } else {
    hb_combined <- merge(hb_combined, y = hb)
  }

  count <- count + 1
}

# Normalize, find variable features, and scale
hb_combined <- NormalizeData(hb_combined)
hb_combined <- FindVariableFeatures(hb_combined)
hb_combined <- ScaleData(hb_combined)

# PCA and Harmony integration
hb_combined <- RunPCA(hb_combined, verbose = FALSE)
hb_combined <- RunHarmony(hb_combined, group.by.vars = "sample")

# UMAP and clustering
hb_combined <- RunUMAP(hb_combined, reduction = "harmony", dims = 1:30)
hb_combined <- FindNeighbors(hb_combined, reduction = "harmony", dims = 1:30)
hb_combined <- FindClusters(hb_combined, resolution = 0.5)

# Save the final Seurat object
saveRDS(hb_combined, file = "results/hb_combined_seurat.rds")

# Plot UMAP
umap_plot <- DimPlot(hb_combined, group.by = "sample", pt.size = 0.5, label = TRUE) +
  ggtitle("UMAP Colored by Sample")
ggsave("results/umap_by_sample.png", umap_plot, width = 8, height = 6)

cat("✅ Analysis complete. Seurat object saved.\n")
