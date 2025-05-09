library(Seurat)
library(dplyr)

# Load Seurat object
hb_combined <- readRDS("results/hb_combined_seurat.rds")

# Generate QC summary table
qc_table <- hb_combined@meta.data %>%
  group_by(sample) %>%
  summarise(
    total_cells = n(),
    avg_nFeature_RNA = round(mean(nFeature_RNA)),
    avg_nCount_RNA = round(mean(nCount_RNA)),
    avg_percent_mt = round(mean(percent.mt), 2)
  )

# Print the table
print(qc_table)

# Optionally save to CSV
write.csv(qc_table, file = "results/qc_summary_table.csv", row.names = FALSE)

