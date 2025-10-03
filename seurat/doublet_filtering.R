# scRNA doublet finder pipeline
# run `conda activate scdblfinder_env` to activate the environment with required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(SingleCellExperiment) # required for scDblFinder
library(BiocManager)
library(scDblFinder)
packageVersion("scDblFinder")  
# Load the Seurat object
seurat_obj <- readRDS("/home/raquelcr/seurat/outs_20250725_1854/cerebral_cortex_prefiltered.rds")

layer_names <- Layers(seurat_obj, assay = "RNA")
print(layer_names)
# Initialize combined matrix using the union of all features
all_genes <- unique(unlist(lapply(layer_names, function(layer) {
  rownames(LayerData(seurat_obj, assay = "RNA", layer = layer))
})))
all_cells <- unique(unlist(lapply(layer_names, function(layer) {
  colnames(LayerData(seurat_obj, assay = "RNA", layer = layer))
})))

combined_counts <- matrix(0,
                         nrow = length(all_genes),
                         ncol = length(all_cells),
                         dimnames = list(all_genes, all_cells))

# Safe layer merging
for (layer in layer_names) {
  layer_counts <- LayerData(seurat_obj, assay = "RNA", layer = layer)
  
  # Find overlapping genes and cells
  common_genes <- intersect(rownames(layer_counts), all_genes)
  common_cells <- intersect(colnames(layer_counts), all_cells)
  
  if (length(common_genes) > 0 && length(common_cells) > 0) {
    # Convert to dense only for the overlapping subset
    layer_subset <- as.matrix(layer_counts[common_genes, common_cells, drop = FALSE])
    combined_counts[common_genes, common_cells] <- 
      combined_counts[common_genes, common_cells] + layer_subset
  }
}

# Create SingleCellExperiment (filter empty genes/cells if needed)
library(SingleCellExperiment)
sce <- SingleCellExperiment(
  assays = list(counts = combined_counts[rowSums(combined_counts) > 0, 
                                         colSums(combined_counts) > 0]),
  colData = seurat_obj[[]][colnames(combined_counts)[colSums(combined_counts) > 0], ]
)
# Add metadata from Seurat object
sce$orig.ident <- seurat_obj$orig.ident[colnames(sce)]

# Add reduced dimensions (e.g., PCA, UMAP)
reducedDims(sce) <- seurat_obj@reductions
sce <- scDblFinder(
  sce,
  clusters = TRUE,               # Auto-cluster (recommended)
  dbr = 0.02,                   # Expected doublet rate (adjust for 10x)
  nfeatures = 2000,             # Use top highly variable genes
  returnType = "sce",           # Returns annotated SCE object
  verbose = TRUE
)

# Filter doublets

singlets <- which(sce$scDblFinder.class == "singlet")
seurat_filtered <- seurat_obj[, singlets]
seurat_filtered

# Transfer doublet scores/class to original Seurat object before filtering
seurat_filtered$scDblFinder.score <- sce$scDblFinder.score
seurat_filtered$scDblFinder.class <- sce$scDblFinder.class
seurat_filtered

# Now filter
seurat_filtered <- subset(seurat_obj, subset = scDblFinder.class == "singlet")

# Check the number of cells before and after filtering
cat("Number of cells before filtering:", ncol(seurat_obj), "\n")
cat("Number of cells after filtering:", ncol(seurat_filtered), "\n")
# Save the filtered Seurat object
output_dir <- "/home/raquelcr/seurat/outs_20250725_1854"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
saveRDS(seurat_filtered, file = paste0(output_dir, "/cerebral_cortex_without_dbl.rds"))


# Here is the same but into a function for reusability

run_doublet_filtering <- function(seurat_obj, output_dir, doublet_rate = 0.02, save = FALSE) {
  ### Function to run doublet filtering on a Seurat object. 
  ### Usually scdblfinder is used for seurat (v3/v4), but this is an adaptation for Seurat v5.
  ### Args:
  ###   seurat_obj: Seurat object (v5)
  ###   output_dir: Directory to save the filtered Seurat object
  ### Returns:
  ###   Filtered Seurat object with doublets removed
  require(Seurat)
  require(scDblFinder)
  require(SingleCellExperiment)
    tissue <- seurat_obj$tissue[1]  # Assuming 'tissue' is a metadata field in the Seurat object
  if (is.null(tissue)) {tissue <- ""}
  layer_names <- Layers(seurat_obj, assay = "RNA")
  print(layer_names)
  # Initialize combined matrix using the union of all features
  all_genes <- unique(unlist(lapply(layer_names, function(layer) {
    rownames(LayerData(seurat_obj, assay = "RNA", layer = layer))
  })))
  all_cells <- unique(unlist(lapply(layer_names, function(layer) {
    colnames(LayerData(seurat_obj, assay = "RNA", layer = layer))
  })))

  combined_counts <- matrix(0,
                          nrow = length(all_genes),
                          ncol = length(all_cells),
                          dimnames = list(all_genes, all_cells))

  # Safe layer merging
  for (layer in layer_names) {
    layer_counts <- LayerData(seurat_obj, assay = "RNA", layer = layer)
    
    # Find overlapping genes and cells
    common_genes <- intersect(rownames(layer_counts), all_genes)
    common_cells <- intersect(colnames(layer_counts), all_cells)
    
    if (length(common_genes) > 0 && length(common_cells) > 0) {
      # Convert to dense only for the overlapping subset
      layer_subset <- as.matrix(layer_counts[common_genes, common_cells, drop = FALSE])
      combined_counts[common_genes, common_cells] <- 
        combined_counts[common_genes, common_cells] + layer_subset
    }
  }

  # Create SingleCellExperiment (filter empty genes/cells if needed)
  library(SingleCellExperiment)
  sce <- SingleCellExperiment(
    assays = list(counts = combined_counts[rowSums(combined_counts) > 0, 
                                          colSums(combined_counts) > 0]),
    colData = seurat_obj[[]][colnames(combined_counts)[colSums(combined_counts) > 0], ]
  )
  # Add metadata from Seurat object
  sce$orig.ident <- seurat_obj$orig.ident[colnames(sce)]

  # Add reduced dimensions (e.g., PCA, UMAP)
  reducedDims(sce) <- seurat_obj@reductions
  sce <- scDblFinder(
    sce,
    clusters = TRUE,               # Auto-cluster (recommended)
    dbr = doublet_rate,                   # Expected doublet rate (adjust for 10x)
    nfeatures = 2000,             # Use top highly variable genes
    returnType = "sce",           # Returns annotated SCE object
    verbose = TRUE
  )

  # Filter doublets

  singlets <- which(sce$scDblFinder.class == "singlet")
  seurat_filtered <- seurat_obj[, singlets]
  seurat_filtered

  # Transfer doublet scores/class to original Seurat object before filtering
  seurat_filtered$scDblFinder.score <- sce$scDblFinder.score
  seurat_filtered$scDblFinder.class <- sce$scDblFinder.class
  seurat_filtered

  # Now filter
  seurat_filtered <- subset(seurat_obj, subset = scDblFinder.class == "singlet")

  # Check the number of cells before and after filtering
  cat("Number of cells before filtering:", ncol(seurat_obj), "\n")
  cat("Number of cells after filtering:", ncol(seurat_filtered), "\n")
  # Save the filtered Seurat object
  if(save) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    saveRDS(seurat_filtered, file = paste0(output_dir, "/", tissue, "_without_dbl.rds"))
  }
  return(seurat_filtered)
}