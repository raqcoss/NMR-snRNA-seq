library(Seurat)
library(tidyverse)
library(patchwork)
library(magrittr)

# Load the Seurat object
seurat_list <- lapply(list.files(path = '/home/raquelcr/seurat/filtered_dataset2', # Adjust folder
                                pattern = "*.rds", 
                                full.names = TRUE), 
                     readRDS)
                     
                     # Get filenames (without path and .rds extension)
tissues <- c("cerebral_cortex", "hippocampus", "midbrain")
# Assign names to the list
names(seurat_list) <- tissues

head(seurat_list$cerebral_cortex)
seurat_list[[1]]
# Assign project.name based on filenames
for (i in seq_along(seurat_list)) {
  seurat_list[[i]]$project.name <- names(seurat_list)[i]
}
# Merge all Seurat objects into one
all_data <- merge(seurat_list$cerebral_cortex, y = seurat_list[c('hippocampus', 'midbrain')], add.cell.ids = tissues)
all_data <- all_data %>%
  NormalizeData(normalization.method = 'LogNormalize', scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = 'vst', nfeatures= 2000) 

# Save data into h5ad file
DefaultAssay(all_data) <- "RNA"
if (length(all_data[["RNA"]]@counts) == 0) {
    all_data[["RNA"]]@counts <- all_data[["RNA"]]@data
}


remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)

SaveLoom(all_data, filename = "nmr_brain.loom")
Convert("nmr_brain.loom", dest = "h5ad", overwrite = TRUE)
all_data <- all_data %>%
  ScaleData() %>%
  RunPCA()

setwd('/home/raquelcr/seurat')
output_dir <- paste0("outs_", format(Sys.time(), "%Y%m%d_%H%M"))
dir.create(output_dir)
print(paste0('Folder ',output_dir, ' was created'))

png(paste0(output_dir,'/elbow_plot_',tissue,'.png'))
ElbowPlot(all_data, ndims = 20) +
      ggtitle(paste(tissue,'Elbow Plot')) +
      theme(plot.title = element_text(size = 10, hjust = 0.5))
dev.off()

chosen_pc = 10
all_data <- all_data %>%
  RunUMAP(dims = 1:chosen_pc) %>%
  FindNeighbors(dims = 1:chosen_pc) %>%
  FindClusters(resolution = 0.6)

source('/home/raquelcr/seurat/ploting_functions.R')
plot_clusters_at_resolutions(all_data,0.6, output_dir,grid_cols = 1)

library(ggplot2)

# Extract UMAP coordinates and metadata
umap_df <- as.data.frame(Embeddings(all_data, "umap"))
umap_df$cluster <- as.factor(Idents(all_data))
umap_df$nCount_RNA <- all_data@meta.data$nCount_RNA

# Scale nCount_RNA to [0.2, 1] for alpha (avoid fully transparent points)
umap_df$alpha <- scales::rescale(umap_df$nCount_RNA, to = c(0.2, 1))
colnames(umap_df)
# Plot
png(paste0(output_dir,'/umap_clusters_',tissue,'_res_',chosen_res,'_alpha-nCounts.png'), width = 1000, height = 800)
plt <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cluster, alpha = alpha)) +
        geom_point(size = 1.2) +
        scale_alpha_identity() +
        theme_minimal() +
        labs(title = paste("UMAP at resolution", chosen_res),
            x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
        guides(alpha = "none")
plt
dev.off()
