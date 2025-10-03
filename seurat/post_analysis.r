library(Seurat)
library(ggplot2)
library(patchwork)

# Reload all objects
seurat_list <- lapply(list.files(path = "filtered_datasets_20250521_0146", # Adjust folder
                                pattern = "*.rds", 
                                full.names = TRUE), 
                     readRDS)

# Get filenames (without path and .rds extension)
file_names <- c("cerebral_cortex",
  "hippocampus", "midbrain")
# Assign names to the list
names(seurat_list) <- file_names

seurat_list$cerebral_cortex

# Assign project.name based on filenames
for (i in seq_along(seurat_list)) {
  seurat_list[[i]]$project.name <- names(seurat_list)[i]
}
seurat_list

# View known markers
known_markers <- c(
'Gja1', 'Mfge8', # Astrocytes
'Laptm5', 'Abca9', # Microglia
'Cldn5', # Endothelial
'Tbr1', 'Gda', # Glutamatergic neurons
'Iptr2', # Oligodendrocytes
'Cd79a', 'H2-Eb1', 'H2-Ab1', 'Mef2c', 'Cr2', 'Fcer2a', 'Jchain', # B-cells
'Ly6d', 'Mzb1', # Memory-B-cells
'Klf9', # Macrophages
'Ctsd', 'Fcer1a', # Mast-cells
'Fos', 'Cd44', 'Fcer1g', # Monocytes
'S100a8', 'S100a9', 'S100a6', 'Retn', # Neutrophils
'Slc40a1', 'Vcam1', # NK-cell & Red-pulp-macrophage
'Ccl5', 'Cd27', 'Cd8b1', 'Cd3g', 'Lat', # T-cells
'Slc46a3' # Timeless-high-dendritic-cell
) 
source('ratopin_functions.R')
output_dir <- paste0("plots_", format(Sys.time(), "%Y%m%d_%H%M"))
dir.create(output_dir)
view_marker_expression_fixed_list(seurat_list, known_markers, output_dir = output_dir, file_name = 'KnownMarkers', omit_violin = TRUE )

## HVG expression

# Find HVGs
top20_list = list()
for (i in seq_along(seurat_list)){
  obj <- seurat_list[[i]]
  HVGs <-VariableFeatures(obj)
  write.csv(HVGs, paste0(output_dir, '/HVG_',obj@project.name,'.csv'))
  top20 <- head(HVGs,20)
  top20_list[[obj@project.name]] <- top20
}
top20_list
view_marker_expression(seurat_list, top20_list, output_dir = output_dir, file_name = 'HVG_Expression', omit_violin = TRUE )




#Junk code

library(tidyverse)
library(celldex)
library(SingleR)
ref <- fetchReference("hpca", "2024-02-26")



# Plot enrichment scores
enrichment_scores <-ggplot(top_pathways, 
  aes(x = reorder(pathway, NES), y = NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() +
  labs(x = "Pathway", y = "NES", title = "Top GSEA Pathways") +
  theme_minimal()

ggsave(paste0(output_dir, "/top_gsea_enrichment_scores.png"),
  plot=enrichment_scores, width = 10, height = 10, 
  dpi = 300, bg='white')
# Plot a specific pathway (e.g., the top hit)
specific_pathways <-plotEnrichment(
  pathway = msigdb_sets[[top_pathways$pathway[1]]],
  stats   = ranked_genes
) +
  labs(title = top_pathways$pathway[1])

ggsave("top_gsea_pathways.png",plot = specific_pathways,  width = 8, height = 6, dpi = 300, bg='white')

####
# Create a list of enrichment plots
enrichment_plots <- lapply(top_pathways$pathway, function(pathway) {
  plotEnrichment(
    pathway = msigdb_sets[[pathway]],
    stats = ranked_genes
  ) +
    labs(title = pathway) +
    theme(plot.title = element_text(size = 8))  # Smaller title for grid
})

# Arrange in 3x3 grid
combined_plot <- wrap_plots(enrichment_plots, ncol = 3) +
  plot_annotation(title = "Top 9 Enriched Pathways from GSEA",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 14)))

# Save the combined plot
ggsave(paste0(output_dir, "/top9_gsea_pathways.png"), 
       plot = combined_plot,
       width = 15,         # Wider to accommodate 3 columns
       height = 12,        # Taller for 3 rows
       bg = 'white',
       dpi = 400)
#################
