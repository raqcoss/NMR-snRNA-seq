library(Seurat, tidyverse, patchwork)
library(magrittr)

# Load the Seurat object
tissue <- 'cerebral_cortex'
output_dir <- paste0("outs_", format(Sys.time(), "%Y%m%d_%H%M"))
{dir.create(output_dir)
print(paste0('Folder ',output_dir, ' was created'))}

seurat_obj <- readRDS(paste0('/home/raquelcr/seurat/outs_20250725_1854/',tissue,'_2nd_filtering.rds'))

#Normalize
seurat_obj <- NormalizeData(seurat_obj, normalization.method = 'LogNormalize', scale.factor = 10000)

#HVG's
merged_seurat <- FindVariableFeatures(seurat_obj, selection.method = 'vst', nfeatures= 2000)
HVGs <-VariableFeatures(merged_seurat)
head(HVGs, 20)
# Scale Data and perform PCA
merged_seurat <- merged_seurat %>%
    ScaleData() %>%
    RunPCA() 
merged_seurat
png(paste0(output_dir,'/elbow_plot_',tissue,'.png'))
ElbowPlot(merged_seurat, ndims = 20) +
      ggtitle(paste(tissue,'Elbow Plot')) +
      theme(plot.title = element_text(size = 10, hjust = 0.5))
dev.off()

# Testing cell scatter (temp)
CellScatter(object = merged_seurat, cell1 = colnames(merged_seurat)[1],cell2 = colnames(merged_seurat)[2] )


# Analyse with chosen PCs at given resolutions
resolutions <- c(0.2,0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7)
n_pc = 10
merged_seurat <- merged_seurat %>%
  RunUMAP(dims = 1:n_pc) %>%
  FindNeighbors(dims = 1:n_pc) %>%
  FindClusters(resolution = resolutions)
source('/home/raquelcr/seurat/ploting_functions.R')
plot_clusters_at_resolutions(merged_seurat,resolutions, output_dir)


# Add cluster column at chosen resolution
chosen_res <- 0.7
Idents(merged_seurat) <- paste0("RNA_snn_res.", chosen_res)

# Get cluster IDs
cluster_ids <- levels(Idents(merged_seurat))

# Get summary stats for each cluster
cluster_summaries <- lapply(cluster_ids, function(clust) {
  obj <- subset(merged_seurat, idents = clust)
  seurat_summary(obj)
})
names(cluster_summaries) <- cluster_ids

# Print or view summaries
cluster_summaries  # This is a named list of summary tables per cluster

library(Seurat)
library(ggplot2)

# Extract UMAP coordinates and metadata
umap_df <- as.data.frame(Embeddings(merged_seurat, "umap"))
umap_df$cluster <- as.factor(Idents(merged_seurat))
umap_df$nCount_RNA <- merged_seurat@meta.data$nCount_RNA

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

#Save object
save_seurat <- function(obj, name) {
  filename <- paste0(output_dir,"/", name, ".rds")
  saveRDS(obj, file = filename)
  return(filename)
}
save_seurat(merged_seurat, paste0(tissue,'_clusters'))

