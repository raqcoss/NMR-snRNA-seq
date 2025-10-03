if (!require("Seurat", quietly = TRUE))
  install.packages("Seurat")
if (!require("harmony", quietly = TRUE))
  install.packages('harmony')
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require('celldex',quietly = TRUE)){BiocManager::install("celldex")}
if (!require("SingleR", quietly = TRUE)){BiocManager::install("SingleR")}
library(Seurat)
library(ggplot2)
library(patchwork)

# 1. Define sample paths (replace with your actual paths)
tissues <- c("cerebral_cortex", "hippocampus", "midbrain")
replicates <- c("1", "2")
num = 0
seurat_list = list()
for (tissue in tissues){
  for (replica in replicates){
    num <- num + 1
    sample_path <- paste0('/home/raquelcr/cellranger_output/out_HetGla_female_1/NMR',num,'_',tissue,'/outs/raw_feature_bc_matrix')
    print(paste0('Loading sample in ' ,sample_path))
    sample_name <- paste0(tissue,'_',replica)
    data <- Read10X(data.dir = sample_path)
    obj <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3, min.features = 200)
    obj$tissue <- tissue
    obj$replicate <- replica
    seurat_list[[sample_name]] <- obj
    print(seurat_list[[sample_name]])
  }
}

# 2. Calculate QC metrics (mitochondrial genes - adjust pattern if needed)
mito_genes <- scan("mito_genes.txt", what = character(), sep = "\n")
print(mito_genes)
for (i in seq_along(seurat_list)) {
  mito_genes_in_assay <- mito_genes[mito_genes %in% rownames(seurat_list[[i]])]
  print(mito_genes_in_assay)
  seurat_list[[i]][["percent.mt"]] <- PercentageFeatureSet(
    seurat_list[[i]], 
    features = mito_genes_in_assay,
    assay = 'RNA')}
#### Check-point

###### QC - MITO vizualization
# Create a list to store combined metric plots
thresholds <- list(
  nFeature_RNA = c(200, 5000),
  nCount_RNA = c(150, 12000),
  percent.mt = c(0, 2))
threshold_exception = list(percent.mt = list(midbrain = 2))
metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

# VIZUALIZATION
output_dir <- paste0("plots_", format(Sys.time(), "%Y%m%d_%H%M"))
{dir.create(output_dir)
print(paste0('Folder ',output_dir, ' was created'))}
source('/home/raquelcr/seurat/ploting_functions.R')
old_plot_qc_metrics(seurat_list = seurat_list, tissues= tissues, metrics = metrics, thresholds= thresholds)

plot_qc_metrics(seurat_list, file_name = paste0(output_dir,'/QC_metrics_raw.jpeg'),show_points=TRUE,
thresholds = thresholds )
plot_qc_metrics(seurat_list, file_name = paste0(output_dir,'/QC_metrics_raw_only_violin.jpeg'),
thresholds = thresholds)
######


# We check out the plots and adjust thresholds to our data
# 5. Filter low-quality cells
source('/home/raquelcr/seurat/ratopin_functions.R')
filtered_seurat_list <- filter_data(seurat_list, thresholds= thresholds)

# Vizualize filtered data
plot_qc_metrics(filtered_seurat_list, file_name = paste0(output_dir,'/QC_metrics_filtered.jpeg'),show_points=TRUE,
thresholds = thresholds , special_thresholds = NULL)
plot_qc_metrics(filtered_seurat_list, file_name = paste0(output_dir,'/QC_metrics_filtered_only_violin.jpeg'),
thresholds = thresholds , special_thresholds = NULL)

###### Checkpoint

# 6. Normalize and find variable features
filtered_seurat_list <- lapply(filtered_seurat_list, function(x) {
  x <- NormalizeData(x)
})

# Merging
# 2. Merge both replicas per tissue (and normalize again)
merged_seurat_list <- list()
for (i in seq_along(tissues)){
  obj <- merge(filtered_seurat_list[[2*i-1]], y = filtered_seurat_list[[2*i]], add.cell.ids = c("Rep1", "Rep2"), project = tissues[i])
  merged_seurat_list[tissues[i]] <-   NormalizeData(obj)
}
merged_seurat_list

# Find HVGs
processed_list <- lapply(merged_seurat_list, function(obj){
  obj <- FindVariableFeatures(obj, selection.method = 'vst', nfeatures= 2000)
})
top20_list = list()
for (i in seq_along(processed_list)){
  obj <- processed_list[[i]]
  HVGs <-VariableFeatures(obj)
  write.csv(HVGs, paste0(output_dir, '/HVG_',obj@project.name,'.csv'))
  top20 <- head(HVGs,20)
  top20_list[[obj@project.name]] <- top20
  top20_plot <- VariableFeaturePlot(obj)
  top20_plot <- LabelPoints(plot=top20_plot, points =top20, repel = TRUE) 
  ggsave(filename = paste0(output_dir,'/top15HVG_', obj@project.name ,'.png'), 
  plot = top20_plot,width = 10, height = 10, dpi = 200, bg = 'white')
}


# Process independently (without harmony batch merging)
seurat_list <- lapply(processed_list, function(obj) {
  obj %>%
    ScaleData() %>%
    RunPCA()        # PCA on each object
})
# Choose number of PC's by elbow method
plot_pca_elbows(seurat_list,file_name = paste0(output_dir,'/pca_scree_plots.jpeg'), ndims = 20)
chosen_pc = list(10,10,5) # Move to elbow in each sree plot

# Analyse with chosen PCs at given resolutions
resolutions <- c(0.2, 0.3, 0.35, 0.4, 0.5)
for (i in seq_along(seurat_list)) {
  n_pc <- chosen_pc[[i]]
  seurat_list[[i]] <- seurat_list[[i]] %>%
    RunUMAP(dims = 1:n_pc) %>%
    FindNeighbors(dims = 1:n_pc) %>%
    FindClusters(resolution = resolutions)
}
source('/home/raquelcr/seurat/ploting_functions.R')
#plot_clusters_by_res_for_all_samples(seurat_list, resolutions, output_dir)
# Create a grid of all resolutions for one sample
plot_clusters_by_sample_at_resolutions(seurat_list = seurat_list, 
resolutions = resolutions , output_dir = output_dir)

# Group cells per corresponding selected resolution 
chosen_res <- c(0.35,0.35,0.35) # same index as in sample from 1 to 3
plot_clusters_at_chosen_res(seurat_list, chosen_res, output_dir)
for (i in seq_along(seurat_list)){
  cluster_col <- paste0("RNA_snn_res.", chosen_res[[i]])
  Idents(seurat_list[[i]]) <- cluster_col
}

# For each cluster in each sample: 
# Find all positive markers for each cluster and filter for significant markers (adj. p < 0.05)
source('/home/raquelcr/seurat/ratopin_functions.R')

# Check basic object info
print(seurat_list[[1]])
head(colnames(seurat_list[[1]]))
seq_along(table(Idents(seurat_list[[1]])))-1

seurat_list <- lapply(seurat_list,function(obj){
  JoinLayers(obj)})

# 4. Find markers
markers_1 <- FindAllMarkers(
  object = seurat_list[[1]],
  assay = "RNA",  
  only.pos = FALSE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)
write.csv(markers_1, paste0(output_dir, '/markers_all_clusters_',tissues[1],".csv"))

markers_2 <- FindAllMarkers(
  object = seurat_list[[2]],
  assay = "RNA",    
  only.pos = FALSE,
  min.pct = 0.25, 
  logfc.threshold = 0.25,
  test.use = "wilcox"
)
write.csv(markers_2, paste0(output_dir, '/markers_all_clusters_',tissues[2],".csv"))

markers_3 <- FindAllMarkers(
  object = seurat_list[[3]],
  assay = "RNA",       
  only.pos = FALSE,
  min.pct = 0.25, # Try lowering if no genes found
  logfc.threshold = 0.25,
  test.use = "wilcox"
)
write.csv(markers_3, paste0(output_dir, '/markers_all_clusters_',tissues[3],".csv"))

all_markers <- list(markers_1, markers_2, markers_3)
names(all_markers) <- tissues
lapply(all_markers,head)

  # Get top 5 markers per cluster
top_markers <- lapply(all_markers, function(markers){
  top_markers <- markers %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>% # sorted by the most DE genes
    slice_head(n = 5) %>%
    ungroup()
  return(unique(top_markers$gene))
})
head(top_markers)
# For specific markers between 2 clusters
# markers01 <- FindMarkers(seurat_list[[1]], ident.1 = 0, ident.2 = 1)
#head(markers01,10)


source('/home/raquelcr/seurat/ratopin_functions.R')
view_marker_expression(seurat_list= seurat_list, top_markers, output_dir = output_dir,ncol = 5, omit_violin = TRUE)

# view_different_marker_expression(processed_list, top_found_marker_genes,output_dir,ncol=8,omit_violin = TRUE)

### CELL ANNOTATION
#CELL ANNOTATION BASED ON MARKER DATABASES

##### Test
# Load required packages
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(dplyr)
library(Seurat)

# 1. Prepare cluster markers ---------------------------------------------------
lapply(all_markers, head)
# 2. Prepare gene sets --------------------------------------------------------
# Get MSigDB gene sets (C8 = cell type signatures, H = hallmarks, C5 = GO terms)
msigdb_sets <- msigdbr(species = "Homo sapiens", collection = "C8") %>% 
  split(x = .$gene_symbol, f = .$gs_name)

# 3. Cluster-specific GSEA ----------------------------------------------------

# Process each marker set with proper error handling
results <- pblapply(names(all_markers), function(tissue_name) {
  markers <- all_markers[[tissue_name]]
  
  # 1. Run GSEA per cluster with enhanced parameters
  cluster_gsea <- lapply(unique(markers$cluster), function(clust) {
    cluster_markers <- markers %>% 
      filter(cluster == clust) %>%
      arrange(desc(avg_log2FC)) %>%
      {setNames(.$avg_log2FC, .$gene)} %>%
      .[!duplicated(names(.))]
    
    # Run GSEA with improved settings
    fgsea_res <- fgsea(
      pathways = msigdb_sets,
      stats = cluster_markers,
      minSize = 15,
      maxSize = 500,
      eps = 0.0,
      nPermSimple = 10000  # Increased from default to handle unbalanced stats
    )
    
    # Convert to clean data frame
    fgsea_res %>%
      as_tibble() %>%
      select(-leadingEdge) %>%  # Remove list column causing write.csv issues
      mutate(cluster = clust, tissue = tissue_name)
  })
  
  # 2. Combine results safely
  gsea_results <- bind_rows(cluster_gsea)
  
  # 3. Save full results after type conversion
  write.csv(
    as.data.frame(gsea_results),  # Ensure pure data frame
    file.path(output_dir, paste0("gsea_results_", tissue_name, ".csv")),
    row.names = FALSE
  )
  
  # 4. Get top pathways with NA handling
  top_pathways <- gsea_results %>%
    group_by(cluster) %>%
    filter(!is.na(padj), padj < 0.05) %>%  # Exclude NA results
    slice_min(pval, n = 5, with_ties = FALSE) %>%
    ungroup()
  
  # 5. Visualization with existence check
  if (nrow(top_pathways) > 0) {
    p <- ggplot(top_pathways, aes(x = factor(cluster), y = reorder(pathway, NES))) +
      geom_point(aes(size = -log10(padj + 1e-10), color = NES), alpha = 0.8) +
      scale_size_continuous(range = c(3, 8)) +
      scale_color_gradient2(
        low = "blue", mid = "white", high = "red", 
        midpoint = 0,
        na.value = "grey50"  # Handle NA values
      ) +
      labs(title = paste("Top Pathways in", tissue_name),
           x = "Cluster", y = "Pathway") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(
      file.path(output_dir, paste0("enriched_pathways_", tissue_name, ".png")),
      plot = p, 
      width = 10, 
      height = 6 + nrow(top_pathways)/4, 
      dpi = 300
    )
  } else {
    message("No significant pathways found for ", tissue_name)
  }
  
  return(gsea_results)
})

# Name the results
names(results) <- names(all_markers)

################################ Hasta aqui
### 1. GSEA (MSigDB) for Functional Enrichment
# Install required packages
if (!requireNamespace("msigdbr", quietly = TRUE)){
   install.packages("msigdbr")
   install.packages('msigdbdf', repos = 'https://igordot.r-universe.dev')}
if (!requireNamespace("fgsea", quietly = TRUE)) install.packages("fgsea")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")

library(msigdbr)
library(fgsea)
library(clusterProfiler)

# Example: Get marker genes from Seurat clusters
source('/home/raquelcr/seurat/ploting_functions.R')

# Get MSigDB gene sets (e.g., "C8" for cell type signatures)
msigdb_sets <- msigdbr(species = "Homo sapiens", collection = "C8") %>% 
  split(x = .$gene_symbol, f = .$gs_name)

# Aggregate duplicate genes (average logFC) and rank genes by logFC 
ranked_genes_list <- lapply(pos_marker_list, function(markers){
  markers %>% 
  group_by(gene) %>%
  summarise(avg_log2FC = mean(avg_log2FC)) %>%
  arrange(desc(avg_log2FC)) %>%
  {setNames(.$avg_log2FC, .$gene)}
  })

ranked_genes_list
# Verify no duplicates
for (i in 1:6){
  stopifnot(!any(duplicated(names(ranked_genes_list[[i]]))))
}
# Run GSEA
gsea_results_list <- lapply(ranked_genes_list,function(ranked_genes){
  fgsea(
  pathways = msigdb_sets,
  stats    = ranked_genes,
  minSize  = 15,
  maxSize  = 500,
  eps      = 0.0  # Disable zero p-value approximation
)})

# Top enriched cell types
head(gsea_results_list[[1]][order(pval), ], 5)
gsea_results_list[[1]][order(pval), ]
# Save full results
for(i in seq_along(gsea_results_list)) {
  sample_name <- processed_list[[i]]@project.name
  
  # Convert to data frame and remove list columns
  ordered_results <- as.data.frame(gsea_results_list[[i]]) %>%
    select(-any_of("leadingEdge")) %>%  # Remove the problematic column
    arrange(pval)  # Order by p-value
  
  write.csv(ordered_results,
           file.path(output_dir, paste0("gsea_results_", sample_name, ".csv")),
           row.names = FALSE)
  
  message("Saved GSEA results for: ", sample_name)
}

# Top 5 enriched pathways
top_pathways_list <- lapply(gsea_results_list,function(gsea_results){
  gsea_results %>%
  arrange(padj) %>%
  head(9)
})

plot_pathways_grid(
  seurat_list = processed_list,
  chosen_res = chosen_res,
  output_dir = output_dir,
  gsea_results_list = gsea_results_list,
  n_top_pathways = 5
)

# Save top pathways
for(i in seq_along(top_pathways_list)) {
  sample_name <- processed_list[[i]]@project.name
  
  # Convert to data frame and remove list columns
  ordered_results <- as.data.frame(top_pathways_list[[i]]) %>%
    select(-any_of("leadingEdge")) %>%  # Remove the problematic column
    arrange(pval)  # Order by p-value
  
  write.csv(ordered_results,
           file.path(output_dir, paste0("top_pathways_", sample_name, ".csv")),
           row.names = FALSE)
  
  message("Saved top pathways results for: ", sample_name)
}

# Save individual enrichment plots
### Aqui me quede antes de salir de viaje ##################\



source('/home/raquelcr/seurat/ploting_functions.R')
plot_annotated_clusters(processed_list,chosen_res, output_dir, gsea_results_list,ranked_genes_list, n_top_pathways = 5)

###########


### 2. Tabula Muris for Reference Mapping

# Installation
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment")
if (!requireNamespace("celldex", quietly = TRUE)) BiocManager::install("celldex")

library(celldex)
library(SingleCellExperiment)

# Load Tabula Muris reference
tm_ref <- TabulaMurisData()  # Mouse data (adjust if human)

# Convert reference to Seurat
tm_seurat <- as.Seurat(tm_ref)

# Run PCA and harmonize with your data
anchors <- FindTransferAnchors(
  reference = tm_seurat,
  query = seurat_list[[1]],  # Your dataset
  dims = 1:30
)

# Transfer labels
predictions <- TransferData(
  anchorset = anchors,
  refdata = tm_seurat$cell_ontology_class,
  dims = 1:30
)

# Add predictions to metadata
seurat_list[[1]]$predicted_celltype <- predictions$predicted.id

### 3. CellMarker 2.0 for Human/Mouse Markers
if (!requireNamespace("CellMarker", quietly = TRUE)) remotes::install_github("ZJUFanLab/CellMarker")
library(CellMarker)

# Get CellMarker database
cell_markers <- CellMarker::cellmarker_data(species = "Human")  # or "Mouse"

# Extract marker genes per cell type
marker_list <- split(cell_markers$GeneSymbol, cell_markers$CellType)

if (!requireNamespace("AUCell", quietly = TRUE)) BiocManager::install("AUCell")
library(AUCell)

# Get expression matrix (log-normalized)
expr_matrix <- GetAssayData(seurat_list[[1]], layer = "data")

# Calculate AUC scores
auc_scores <- AUCell_run(expr_matrix, marker_list)
celltype_scores <- t(assay(auc_scores))

# Assign the best match
seurat_list[[1]]$celltype_auc <- colnames(celltype_scores)[max.col(celltype_scores)]


### EXTRA: Using celldex
library(celldex)
library(SingleR)

# Load reference (e.g., Human Primary Cell Atlas)
ref <- celldex::HumanPrimaryCellAtlasData()

# Convert Seurat to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_list[[1]])

# Annotate
pred <- SingleR(sce, ref, labels = ref$label.main)
seurat_list[[1]]$singler_labels <- pred$labels


#### vizualization
DimPlot(seurat_list[[1]], group.by = "predicted_celltype", label = TRUE)

####
# Define your markers of interest
marker_genes <- list(
  Astrocytes = c('Gja1', 'Mfge8'),
  Microglia = c('Laptm5','Abca9'),
  Endothelial = c('Cldn5'),
  Glutamatergic_neurons = c('Tbr1','Gda'),
  Oligodendrocytes = c('Iptr2'),
  B_cells = c('Cd79a','H2-Eb1', 'H2-Ab1', 'Mef2c','Cr2', 'Fcer2a', 'Jchain'),
  Memory_B_cells = c('Ly6d', 'Mzb1'),
  Macrophages = c('Klf9'),
  Mast_cells = c('Ctsd', 'Fcer1a'),
  Monocytes = c('Fos','Cd44', 'Fcer1g'),
  Neutrophils = c('S100a8', 'S100a9', 'S100a6', 'Retn'),
  NK_cell = c('Slc40a1', 'Vcam1'),
  Red_pulp_macrophage = c('Slc40a1', 'Vcam1'),
  T_cells = c('Ccl5','Cd27', 'Cd8b1', 'Cd3g', 'Lat'),
  Timeless_DC = c('Slc46a3')
)
all_markers <- c(
  'GFAP','Gja1', 'Mfge8',          # Astrocytes
  'CCL4','Laptm5', 'Abca9', # Microglia
  'Cldn5',                  # Endothelial
  'Tbr1', 'Gda',            # Glutamatergic neurons
  'Iptr2',                  # Oligodendrocytes
  'Cd79a', 'H2-Eb1', 'H2-Ab1', 'Mef2c', 'Cr2', 'Fcer2a', 'Jchain',  # B-cells
  'Ly6d', 'Mzb1',           # Memory-B-cells
  'Klf9',                   # Macrophages
  'Ctsd', 'Fcer1a',         # Mast-cells
  'Fos', 'Cd44', 'Fcer1g',  # Monocytes
  'S100a8', 'S100a9', 'S100a6', 'Retn',  # Neutrophils
  'Slc40a1', 'Vcam1',       # NK-cell & Red-pulp-macrophage
  'Ccl5', 'Cd27', 'Cd8b1', 'Cd3g', 'Lat',  # T-cells
  'Slc46a3'                 # Timeless-high-dendritic-cell
)
# Visualize expression of markers in UMAP clusters
source('/home/raquelcr/seurat/ratopin_functions.R')
view_marker_expression(seurat_list, all_markers, output_dir = output_dir)

###########

# Save each processed object with timestamp 
saving_dir <- paste0("filtered_datasets_", format(Sys.time(), "%Y%m%d_%H%M"))
dir.create(saving_dir)
save_seurat <- function(obj, name) {
  filename <- paste0(saving_dir,"/seurat_", name, "_", format(Sys.time(), "%Y%m%d"), ".rds")
  saveRDS(obj, file = filename)
  return(filename)
} 

# Apply to all objects in seurat_list
saved_files <- mapply(save_seurat, seurat_list, names(seurat_list))

###############################

# 10. Find marker genes
markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25)
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Save results
write.csv(top_markers, "cluster_markers.csv", row.names = FALSE)
saveRDS(merged, "integrated_seurat_object.rds")