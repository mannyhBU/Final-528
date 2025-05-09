# Rebuild Seurat object in v4 from raw or lightly processed object

library(Seurat)

# Load original object (update filename if needed)
hb <- readRDS("results/hb_combined_seurat.rds")  # or use hb_combined_raw.rds

# Set assay (if needed)
DefaultAssay(hb) <- "RNA"

# Run preprocessing steps
hb <- NormalizeData(hb, verbose = FALSE)
hb <- FindVariableFeatures(hb, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
hb <- ScaleData(hb, verbose = FALSE)
hb <- RunPCA(hb, npcs = 30, verbose = FALSE)
hb <- FindNeighbors(hb, dims = 1:20, verbose = FALSE)
hb <- FindClusters(hb, resolution = 0.3, verbose = FALSE)
hb <- RunUMAP(hb, dims = 1:20, verbose = FALSE)

# Save output
saveRDS(hb, file = "results/hb_combined_seurat2.rds")
cat("✅ Rebuild complete. Saved as hb_combined_seurat2.rds\n")
