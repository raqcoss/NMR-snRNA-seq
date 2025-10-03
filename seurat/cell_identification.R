library(Seurat)
library(ggplot2)
library(patchwork)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(proxy)
library(RcppAnnoy)
# if(!require('Vennerable',quietly = TRUE)){install.packages("Vennerable", repos="http://R-Forge.R-project.org")}
if(!require('BRETIGEA',quietly=TRUE)){
  install.packages("BRETIGEA", dependencies = TRUE)
}
library(BRETIGEA)
library(knitr) #only for visualization

#EXAMPLE
str(aba_marker_expression, list.len = 10) #input data format
str(aba_pheno_data) #input data format

ct_res = brainCells(aba_marker_expression, nMarker = 50)
kable(head(ct_res)) #output data format

cor_mic = cor.test(ct_res[, "mic"], as.numeric(aba_pheno_data$ihc_iba1_ffpe), method = "spearman")
print(cor_mic)
cor_ast = cor.test(ct_res[, "ast"], as.numeric(aba_pheno_data$ihc_gfap_ffpe), method = "spearman")
print(cor_ast)
# RISmed: para encontrar literature markers
# Analysing all samples
setwd('/home/raquelcr/seurat')

# Reload all objects
seurat_list <- lapply(list.files(path = "filtered_dataset", # Adjust folder
                                pattern = "*.rds", 
                                full.names = TRUE), 
                     readRDS)
                     
                     # Get filenames (without path and .rds extension)
tissues <- c("cerebral_cortex", "hippocampus", "midbrain")
# Assign names to the list
names(seurat_list) <- tissues
# Assign project.name based on filenames
for (i in seq_along(seurat_list)) {
  seurat_list[[i]]$project.name <- tissues[i]
}

# Summary table of objects
source('ratopin_functions.R')
lapply(seurat_list, seurat_summary)
seurat_list

### Cluster verification
source('clustering_qc.R')
head(Idents(seurat_list[[1]]))
head(GetAssayData(seurat_list[[1]],layer='scale.data'))[,1:5]

cluster_similarity_seurat <- function(seurat_obj, method = "euclidean") {
  # Get normalized data
  data <- GetAssayData(seurat_obj,layer='scale.data')
  clusters <- Idents(seurat_obj)
  source('clustering_qc.R')
  if (method == "cosine") {
    return(optimized_cosine_similarity(data, clusters))
  } else if (method == "pearson") {
    return(pearson_similarity_per_cluster(data, clusters))
  } else if (method == "euclidean") {
    return(euclidean_similarity_per_cluster(data, clusters))
  } else if (method == "jaccard") {
    return(jaccard_similarity_per_cluster(data, clusters))
  } else if (method == "annoy") {
    return(annoy_similarity(data, clusters))
  } else {
    stop("Unsupported method")
  }
}
# Get euclidean distances within clusters
library(purrr)
# Get cluster assignments
clusters <- Idents(seurat_list[[1]])
pca_embeddings <- Embeddings(seurat_list[[1]], reduction = "pca")
# Calculate Euclidean distance matrix
dist_matrix <- as.matrix(dist(pca_embeddings, method = "euclidean"))
# Calculate mean within-cluster distance
cluster_distances <- map_dbl(levels(clusters), function(cluster) {
  cells <- WhichCells(seurat_list[[1]], idents = cluster)
  if(length(cells) > 1) {
    mean(dist(pca_embeddings[cells, ], method = "euclidean"))
  } else {
    NA_real_
  }
})
names(cluster_distances) <- levels(clusters)
cluster_distances

library(pheatmap)
pheatmap(dist_matrix, 
         annotation_row = data.frame(Cluster=clusters),
         show_rownames = FALSE,
         show_colnames = FALSE)

# converting mean euclidean distance to similarity
euclidean_similarity <- 1 - (cluster_distances / max(cluster_distances, na.rm = TRUE))

# Combine into a dataframe
similarity_df <- data.frame(
  Cluster = names(cluster_distances),
  Euclidean_Similarity = euclidean_similarity,
  Cosine_Similarity = cosine_similarity,
  Annoy_Similarity = annoy_similarity)

comparison_long <- pivot_longer(comparison_df, -Cluster, 
                               names_to = "Metric", values_to = "Value")

ggplot(comparison_long, aes(x = Cluster, y = Value, color = Metric)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_line(aes(group = Metric), alpha = 0.5) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  labs(title = "Cluster Cohesion Metrics Comparison",
       y = "Normalized Similarity (0-1)",
       x = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Calculate cosine similarity
cosine_similarity_cer_cortex <- cluster_similarity_seurat(seurat_list[[1]],method='cosine')
# Calculate Annoy similarity
annoy_similarity_cer_cortex <- cluster_similarity_seurat(seurat_list[[1]],method = 'annoy')

similarity_df <- data.frame(
  Cluster = names(cosine_similarity_cer_cortex),
  Cosine_Similarity = unlist(cosine_similarity_cer_cortex),
  Annoy_Similarity = unlist(annoy_similarity_cer_cortex)
)
similarity_df
# Vizualize
output_dir <- paste0("similarity_", format(Sys.time(), "%Y%m%d_%H%M"))
{dir.create(output_dir)
print(paste0('Folder ',output_dir, ' was created'))}

comparison_plot <- ggplot(similarity_df, aes(x = Cluster)) +
  # Add points and lines for Cosine Similarity
  geom_point(aes(y = Cosine_Similarity, color = "Cosine"), size = 3) +
  geom_line(aes(y = Cosine_Similarity, group = 1, color = "Cosine"), linewidth = 1) +
  
  # Add points and lines for Annoy Similarity
  geom_point(aes(y = Annoy_Similarity, color = "Annoy"), size = 3) +
  geom_line(aes(y = Annoy_Similarity, group = 1, color = "Annoy"), linewidth = 1) +
  
  # Customize colors and labels
  scale_color_manual(name = "Metric",
                     values = c("Cosine" = "#E69F00", 
                                "Annoy" = "#56B4E9")) +
  
  # Axis labels and theme
  labs(title = "Comparison of Cosine vs Annoy Similarity by Cluster",
       y = "Similarity Score",
       x = "Cluster") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(paste0(output_dir,"/similarity_comparison_",tissue,".png"), 
       plot = comparison_plot,
       device = "png",
       width = 8,          # Width in inches
       height = 6,         # Height in inches
       dpi = 300,          # Resolution
       bg = "white")       # Background color
######

ggplot(similarity_df, aes(x = Cluster, y = Similarity)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Average Cosine Similarity per Cluster",
       y = "Mean Pairwise Cosine Similarity")
cluster_similarity_seurat(seurat_list[[1]],method='cosine')


ggplot(similarity_df, aes(x = Cluster, y = Similarity)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Average NN Similarity per Cluster",
       y = "Mean Pairwise NN Similarity")

#example

# Prepare example data (replace with your actual results)
sim_data <- data.frame(
  Cluster = factor(paste0("C", 1:10)),
  Cosine = runif(10, 0.5, 0.9),
  Pearson = runif(10, 0.4, 0.85),
  Jaccard = runif(10, 0.3, 0.8)
)

# Reshape data for plotting
plot_data <- sim_data %>%
  pivot_longer(cols = -Cluster, names_to = "Metric", values_to = "Similarity")
######## end of example
print(seurat_list[[1]])
head(colnames(seurat_list[[1]]))
head(rownames(seurat_list[[1]]))

# Join layers is needed to find markers
seurat_list <- lapply(seurat_list,JoinLayers)

all_markers <- lapply(seurat_list,function(obj){
  FindAllMarkers(
    object = obj,
    assay = "RNA",  
    only.pos = FALSE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox")
})
head(all_markers[[1]])
#Create new output dir
output_dir <- paste0("outs_", format(Sys.time(), "%Y%m%d_%H%M"))
{dir.create(output_dir)
print(paste0('Folder ',output_dir, ' was created'))}

# Save 
for(i in seq_along(all_markers)){
  write.csv(all_markers[i], paste0(output_dir, '/markers_all_clusters_',tissues[i],'.csv'))
}
# Identify cell types according to BRETIGEA
expr_matrix <- as.matrix(GetAssayData(seurat_list[[1]], assay = "RNA", slot = "data"))
unique(markers_df_brain$cell)
levels(markers_df_brain$cell)

celltype_predictions <- assign_cell_types(
  expr_mat = expr_matrix,
  markers = markers_df_brain,
  cluster_ids = Idents(seurat_list[[1]])  # Replace with your cluster column name
)
head(Idents(seurat_list[[1]]))
str(LayerData(seurat_list[[1]], assay = "RNA", layer = "counts"), list.len = 20) #input data format
str(seurat_list[[1]][["RNA"]]$data)
str(aba_pheno_data) #input data format

ct_res = brainCells(expr_matrix, nMarker = 50)
kable(head(ct_res)) #output data format

cor_mic = cor.test(ct_res[, "mic"], as.numeric(aba_pheno_data$ihc_iba1_ffpe), method = "spearman")
print(cor_mic)
cor_ast = cor.test(ct_res[, "ast"], as.numeric(aba_pheno_data$ihc_gfap_ffpe), method = "spearman")
print(cor_ast)

literature_markers <- lapply(list.files(path = "/home/raquelcr/seurat/literature_markers", # Adjust folder
                                pattern = "*.csv", 
                                full.names = TRUE), 
                     read.csv)
# Get literature markers

names(literature_markers) <- tissues
literature_markers
known_marker_list <- list()
all_known_markers <- list()
for(i in seq_along(literature_markers)){
  known_marker_list[[tissues[i]]]<- literature_markers[[i]]$marcador[literature_markers[[i]]$marcador %in% rownames(seurat_list[[i]])]
  union(all_known_marker_list,known_marker_list[[tissues[i]]])
}
known_marker_list

source('/home/raquelcr/seurat/ratopin_functions.R')
view_marker_expression(seurat_list= seurat_list, known_marker_list, output_dir = output_dir,ncol = 5, omit_violin = TRUE)



# Show top 20 HVGs in plot
top20_list = list()
for (i in seq_along(seurat_list)){
  obj <- seurat_list[[i]]
  HVGs <-VariableFeatures(obj)
  write.csv(HVGs, paste0(output_dir, '/HVG_',tissues[i],'.csv'))
  top20 <- head(HVGs,20)
  top20_list[[tissues[i]]] <- top20
  top20_plot <- VariableFeaturePlot(obj)
  top20_plot <- LabelPoints(plot=top20_plot, points =top20, repel = TRUE) 
  ggsave(filename = paste0(output_dir,'/top15HVG_', tissues[i] ,'.png'), 
  plot = top20_plot,width = 10, height = 10, dpi = 200, bg = 'white')
}
top20_list
obj
seurat_list['cerebral_cortex']
names(seurat_list)

source('/home/raquelcr/seurat/ratopin_functions.R')
view_marker_expression(seurat_list= seurat_list, known_marker_list, output_dir = output_dir,ncol = 5, omit_violin = TRUE)

