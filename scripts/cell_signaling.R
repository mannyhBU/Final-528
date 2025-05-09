#!/usr/bin/env Rscript

# Load libraries
library(igraph)
library(ggraph)
library(ggplot2)

# Define cell types (nodes)
cell_types <- c("Tumor Cells", "Cholangiocytes", "Hepatoblasts", 
                "Endothelial Cells", "Stromal Cells", "Immune Cells")

# Define directed edges (signaling relationships)
edges <- data.frame(
  from = c("Tumor Cells", "Tumor Cells", "Endothelial Cells", "Stromal Cells", 
           "Hepatoblasts", "Immune Cells", "Tumor Cells"),
  to   = c("Immune Cells", "Endothelial Cells", "Stromal Cells", "Tumor Cells", 
           "Cholangiocytes", "Tumor Cells", "Cholangiocytes")
)

# Create graph object
g <- graph_from_data_frame(edges, vertices = cell_types, directed = TRUE)

# Plot using ggraph
p <- ggraph(g, layout = "circle") +
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')), 
                 end_cap = circle(6, 'mm'), edge_colour = "gray60") +
  geom_node_circle(aes(r = 0.1), fill = "lightblue", color = "black") +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  theme_void() +
  ggtitle("Simplified Cell-Cell Signaling Network (HB Tumor Context)")

# Save plot
ggsave("results/simple_cell_signaling_network.png", p, width = 8, height = 6)

# Show message
cat("✅ Signaling network plot saved to results/simple_cell_signaling_network.png\n")

