library(Seurat)
library(dplyr)
library(ggplot2)

# Load your labeled Seurat object
hb_combined <- readRDS("results/hb_combined_seurat_labeled.rds")

# Make sure manual labels were applied correctly
print(table(Idents(hb_combined)))  # This should show multiple cell types

# Create data frame with sample and cell type info
meta_df <- hb_combined@meta.data
meta_df$cell_type <- Idents(hb_combined)

# Group by sample and cell type
prop_df <- meta_df %>%
  group_by(sample = orig.ident, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count))

# Plot cell type proportions by sample
p <- ggplot(prop_df, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Cell Type Proportions by Sample", x = "Sample", y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("results/cell_type_proportions_by_sample.png", p, width = 10, height = 6)
