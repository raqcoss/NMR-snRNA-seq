# Summary of object
seurat_summary <- function(obj){
  metadata <- obj@meta.data
  summary_table <- data.frame(
    Metric = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'TF_counts_sum'),
    Min = c(min(metadata$nFeature_RNA),
      min(metadata$nCount_RNA),
      min(na.omit(metadata$percent.mt)),
      min(metadata$TF_counts_sum)
      ),
    Max = c(max(metadata$nFeature_RNA),
      max(metadata$nCount_RNA),
      max(na.omit(metadata$percent.mt)),
      max(metadata$TF_counts_sum)
      ),
    Mean = c(
      mean(metadata$nFeature_RNA),
      mean(metadata$nCount_RNA),
      mean(na.omit(metadata$percent.mt)),
      mean(metadata$TF_counts_sum)
    ),
    Median = c(
      median(metadata$nFeature_RNA),
      median(metadata$nCount_RNA),
      median(na.omit(metadata$percent.mt)),
      median(metadata$TF_counts_sum)
    ),
    SD = c(
      sd(metadata$nFeature_RNA),
      sd(metadata$nCount_RNA),
      sd(na.omit(metadata$percent.mt)),
      sd(metadata$TF_counts_sum)
    )
    
  )
  return(summary_table)
}



filter_data <- function(
    seurat_list,
    thresholds = list(
      nFeature_RNA = c(200, 5000),
      nCount_RNA = c(150, 12000),
      percent.mt = c(0, 0.20)
    ),
    special_thresholds = list(
      percent.mt = list(midbrain = 0.6)
    )
) {
  # Calculate percent.mt if not already present
  if (!"percent.mt" %in% names(seurat_list[[1]]@meta.data)) {
    mito_genes = c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB")
    seurat_list <- lapply(seurat_list, function(x) {
      x[["percent.mt"]] <- PercentageFeatureSet(x, features = mito_genes)
      return(x)
    })
  }
  
  # Apply filtering with special thresholds
  filtered_list <- lapply(seurat_list, function(x) {
    # Get tissue from metadata (assuming it's stored in x$tissue)
    tissue <- x$tissue[1]  # [1] gets the first value (all should be same)
    
    # Determine the appropriate percent.mt threshold
    mt_threshold <- if (tissue %in% names(special_thresholds$percent.mt)) {
      special_thresholds$percent.mt[[tissue]]
    } else {
      thresholds$percent.mt[2]
    }
    
    # Apply filtering
    subset(x, subset = 
             nFeature_RNA > thresholds[['nFeature_RNA']][1] & 
             nFeature_RNA < thresholds[['nFeature_RNA']][2] & 
             nCount_RNA > thresholds[['nCount_RNA']][1] &
             nCount_RNA < thresholds[['nCount_RNA']][2] &
             percent.mt > thresholds[['percent.mt']][1] &
             percent.mt < mt_threshold)
  })
  
  # Print filtering summary
  original_cells <- sum(sapply(seurat_list, ncol))
  filtered_cells <- sum(sapply(filtered_list, ncol))
  message(sprintf(
    "Filtering results:\n  Original cells: %d\n  Filtered cells: %d\n  Kept: %.1f%%",
    original_cells, filtered_cells, (filtered_cells/original_cells)*100
  ))
  
  return(filtered_list)
}

find_markers_in_all_objects <- function(seurat_list, output_dir, min_pct = 0.25,
    logfc_threshold = 0.25, p_value_threshold = 0.05){
  require(pbapply)
  require(Seurat)
  pos_marker_list <- pblapply(seq_along(seurat_list), function(i) {
    obj <- seurat_list[[i]]
    sample_name <- obj@project.name
    
    # Find all positive markers for each cluster
    all_markers <- FindAllMarkers(
      object = obj,
      only.pos = TRUE,
      min.pct = min_pct,
      logfc.threshold = logfc_threshold,
      test.use = "wilcox",
      assay = 'RNA')
    # Filter for significant markers (adj. p < 0.05)
    significant_markers <- all_markers[all_markers$p_val_adj < p_value_threshold, ]
    # Save to CSV
    write.csv(
      significant_markers, 
      paste0(output_dir, '/positive_markers_all_clusters_', sample_name, ".csv")
    )
    return(significant_markers)
})
return(pos_marker_list)
}

view_marker_expression <- function(seurat_list, marker_gene_list, output_dir, file_name='MarkerExpressions',
      ncol = 5, colors = c("lightgrey","#0000ff", "#ff0000"), omit_violin = FALSE , samples = c('cerebral_cortex','hippocampus', 'midbrain')){
lapply(marker_gene_list,toupper)  # Replace with your genes
for (i in seq_along(seurat_list)){
  obj <- seurat_list[[i]]
  sample_name <- samples[i]
  marker_genes <- marker_gene_list[[i]]
  markers <- marker_genes[marker_genes %in% rownames(obj)]
  # Feature plots
  fp <- FeaturePlot(obj, features = markers, ncol = ncol, cols = colors)
  ggsave(
    filename = paste0(output_dir,'/', file_name, 'InUMAP_', sample_name, ".png"),
    plot = fp, width = 3.5*ncol, height = 3.2*(ceiling(length(markers)/ncol)), dpi = 300)
  message(paste0('Saved Marker Expression Plot of ', sample_name))
  # Violin plots
  if (!omit_violin){
  vp <- VlnPlot(obj, features = markers, ncol = ncol)
  ggsave(
    filename = paste0(output_dir,'/',file_name, 'InViolin_', sample_name, ".png"),
    plot = vp, width = 17, height = 16, dpi = 300)
  message(paste0('Saved Marker Expression Violin Plot of ', sample_name))
  }
}
}

view_marker_expression_fixed_list <- function(seurat_list, marker_genes, output_dir,file_name='MarkerExpressions', 
      ncol = 5, colors = c("lightgrey","#0000ff", "#ff0000"), omit_violin = FALSE ){
marker_genes <- toupper(marker_genes) 
for (i in seq_along(seurat_list)){
  obj <- seurat_list[[i]]
  sample_name <- obj@project.name
  markers <- marker_genes[marker_genes %in% rownames(obj)]
  # Feature plots
  fp <- FeaturePlot(obj, features = markers, ncol = ncol, cols = colors)
  ggsave(
    filename = paste0(output_dir,'/02_IndMarkerExpressionsInUMAP_', sample_name, ".png"),
    plot = fp, width = 3.5*ncol, height = 3.2*(ceiling(length(markers)/ncol)), dpi = 300)
  message(paste0('Saved Marker Expression Plot of ', sample_name))
  # Violin plots
  if (!omit_violin){
  vp <- VlnPlot(obj, features = markers, ncol = ncol)
  ggsave(
    filename = paste0(output_dir,'/',file_name,'_', sample_name, ".png"),
    plot = vp, width = 17, height = 16, dpi = 300)
  message(paste0('Saved Marker Expression Violin Plot of ', sample_name))
  }
}
}
view_different_marker_expression <- function(seurat_list, marker_gene_list, output_dir, 
                                           ncol = 5, colors = c("lightgrey", "#0000ff", "#ff0000"), 
                                           omit_violin = FALSE) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Convert all genes to uppercase and combine markers by tissue type
  marker_gene_list <- lapply(marker_gene_list, toupper)
  
  # Combine markers for each tissue type (pairs of replicates)
  combined_markers <- list()
  for (i in seq(1, length(marker_gene_list), by = 2)) {
    tissue_markers <- union(marker_gene_list[[i]], marker_gene_list[[i+1]])
    combined_markers[[ceiling(i/2)]] <- tissue_markers
  }
  
  # Process each Seurat object
  for (i in seq_along(seurat_list)) {
    obj <- seurat_list[[i]]
    sample_name <- obj@project.name
    tissue_idx <- ceiling(i/2)  # Determine tissue index (1, 2, or 3)
    
    # Get markers present in this object
    markers <- combined_markers[[tissue_idx]]
    markers <- markers[markers %in% rownames(obj)]
    
    # Skip if no markers found
    if (length(markers) == 0) {
      warning(paste("No matching markers found for sample:", sample_name))
      next
    }
    
    message(paste("Plotting", length(markers), "markers for sample:", sample_name))
    
    # Feature plots
    fp <- FeaturePlot(obj, features = markers, ncol = ncol, cols = colors)
    ggsave(
      filename = file.path(output_dir, paste0('02_IndMarkerExpressionsInUMAP_', sample_name, ".png")),
      plot = fp, 
      width = 3.5 * min(ncol, length(markers)),  # Adjust width for fewer markers
      height = 3.2 * ceiling(length(markers)/ncol), 
      dpi = 300
    )
    
    # Violin plots (if not omitted)
    if (!omit_violin && length(markers) > 0) {
      vp <- VlnPlot(obj, features = markers, ncol = ncol, pt.size = 0.1)
      ggsave(
        filename = file.path(output_dir, paste0('03_IndMarkerExpressionsInViolin_', sample_name, ".png")),
        plot = vp, 
        width = min(17, 3 * length(markers)),  # Adjust width based on marker count
        height = 16, 
        dpi = 300
      )
    }
  }
}