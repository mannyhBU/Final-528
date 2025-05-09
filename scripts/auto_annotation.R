#!/usr/bin/env Rscript

# Load libraries
suppressMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(SingleCellExperiment)
})

# Load Seurat object
seurat_file <- "results/hb_combined_seurat_clustered.rds"
hb_combined <- readRDS(seurat_file)

# Check available layers
assay <- hb_combined[["RNA"]]
layer_names <- Layers(assay)
cat("Available layers in RNA assay:\n")
print(layer_names)

# Pick a reasonable layer (prefer 'counts', else fallback to first available)
chosen_layer <- if ("counts" %in% layer_names) "counts" else layer_names[1]
cat(paste("Using layer:", chosen_layer, "\n"))

# Extract counts matrix
counts_mat <- GetAssayData(assay, layer = chosen_layer)

# Check matrix shape
cat("✅ Layer chosen:", chosen_layer, "\n")
cat("📐 counts_mat dims:", dim(counts_mat), "\n")
cat("📐 full meta.data dims:", dim(hb_combined@meta.data), "\n")

# Subset metadata to match counts matrix cells
# Subset metadata to match cells
meta <- hb_combined@meta.data[colnames(counts_mat), ]

# Build SCE
sce <- SingleCellExperiment(
  assays = list(counts = counts_mat),
  colData = meta
)

# Create logcounts for SingleR
logcounts(sce) <- log1p(counts(sce))

# Load local reference
ref <- readRDS("data/HumanPrimaryCellAtlasData.rds")

# Run SingleR
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)



# Add labels to Seurat object
hb_combined$SingleR_label <- pred$labels

# Save annotated object
saveRDS(hb_combined, file = "results/hb_combined_seurat_annotated.rds")
