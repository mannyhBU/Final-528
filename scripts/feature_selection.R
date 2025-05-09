# feature_selection.R

library(Seurat)
library(ggplot2)

# Load the combined Seurat object
hb_combined <- readRDS("results/hb_combined_seurat.rds")

# Identify highly variable features using the vst method
hb_combined <- FindVariableFeatures(hb_combined, selection.method = "vst", nfeatures = 2000)

# Save the top variable features plot
variable_feature_plot <- VariableFeaturePlot(hb_combined)
variable_feature_plot_labeled <- LabelPoints(plot = variable_feature_plot, points = head(VariableFeatures(hb_combined), 10), repel = TRUE)

# Save plot to file
ggsave("results/highly_variable_features.png", variable_feature_plot_labeled, width = 8, height = 6)

# Print basic summary to terminal
total_genes <- nrow(hb_combined)
num_variable_features <- length(VariableFeatures(hb_combined))
num_non_variable <- total_genes - num_variable_features

cat("\u2705 Feature selection complete:\n")
cat("  Total genes: ", total_genes, "\n")
cat("  Highly variable genes: ", num_variable_features, "\n")
cat("  Non-variable genes: ", num_non_variable, "\n")

# Save variable gene info as a CSV table for reference
variable_genes <- VariableFeatures(hb_combined)
write.csv(data.frame(gene = variable_genes), file = "results/variable_genes.csv", row.names = FALSE)