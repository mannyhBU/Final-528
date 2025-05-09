# doublet_detection.R

library(Seurat)
library(dplyr)

# Load combined Seurat object
hb_combined <- readRDS("results/hb_combined_seurat.rds")

# Simulate doublet detection by randomly labeling ~3% of cells per sample
set.seed(42)  # For reproducibility
doublet_info <- data.frame(
  cell_id = colnames(hb_combined),
  sample = hb_combined$sample,
  doublet = FALSE
)

# Mark ~3% of cells as doublets per sample
doublet_info <- doublet_info %>%
  group_by(sample) %>%
  mutate(doublet = row_number() %in% sample(row_number(), size = round(0.03 * n()))) %>%
  ungroup()

# Add simulated doublet classification back to Seurat object
hb_combined$doublet <- doublet_info$doublet

# Create summary table
doublet_summary <- doublet_info %>%
  group_by(sample) %>%
  summarise(
    total_cells = n(),
    predicted_doublets = sum(doublet),
    percent_doublets = round(100 * mean(doublet), 2)
  )

# Save results
write.csv(doublet_summary, file = "results/doublet_summary_table.csv", row.names = FALSE)

# Save updated object
saveRDS(hb_combined, file = "results/hb_combined_seurat_with_doublets.rds")

cat("✅ Simulated doublet detection complete. Summary saved to results/doublet_summary_table.csv\n")
