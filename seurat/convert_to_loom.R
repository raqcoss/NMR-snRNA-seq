library(Seurat)
library(SeuratDisk)

# Input/Output dirs
input_dir  <- "/home/raquelcr/seurat/filtered_dataset2"
output_dir <- "/home/raquelcr/scanpy/nmr_data"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Find .rds files
rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
samples   <- c('cerebral_cortex', 'hippocampus', 'midbrain')

# Loop over all samples
for (i in seq_along(rds_files)) {
  message("Processing sample: ", samples[i])
  seurat_obj <- readRDS(rds_files[i])

  # Get all RNA layers
  assay <- seurat_obj[["RNA"]]
  layers <- LayerNames(assay)  # e.g. counts.sample_1, counts.sample_2, data.sample_1, data.sample_2

  # Identify counts/data pairs by suffix (_1, _2, etc.)
  suffixes <- unique(gsub("^(counts|data)\\.", "", layers))

  for (suf in suffixes) {
    message("  Exporting layer: ", suf)

    # Extract counts + data for this layer
    counts <- LayerData(assay, layer = paste0("counts.", suf))
    data   <- LayerData(assay, layer = paste0("data.", suf))

    # Build a v4-style assay
    rna_assay <- CreateAssayObject(counts = counts)
    rna_assay@data <- data

    # New Seurat object with metadata preserved
    obj <- seurat_obj
    obj[["RNA"]] <- rna_assay
    obj$layer_id <- suf  # track which layer it came from

    # Save as .loom
    SaveLoom(
      object   = obj,
      filename = file.path(output_dir, paste0(samples[i], "_", suf, ".loom")),
      overwrite = TRUE
    )
  }
    # Convert to .h5ad
    loom_files <- list.files(output_dir, pattern = "\\.loom$", full.names = TRUE)
    for (loom in loom_files) {
        h5ad_file <- sub("\\.loom$", ".h5ad", loom)
        message("  Converting to h5ad: ", h5ad_file)
        Convert(loom, dest = "h5ad", overwrite = TRUE)
    }
    }
