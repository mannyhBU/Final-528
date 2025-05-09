library(Seurat)

hb <- readRDS("results/hb_combined_seurat_clustered.rds")

# Try rebuilding the RNA assay
DefaultAssay(hb) <- "RNA"
hb <- NormalizeData(hb)
hb <- FindVariableFeatures(hb)
hb <- ScaleData(hb)

# Save the fixed object
saveRDS(hb, file = "results/hb_combined_seurat2.rds")
cat("✅ RNA assay repaired and saved as hb_combined_seurat2.rds\n")
