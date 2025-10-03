library(Seurat)
library(proxy)
library(RcppAnnoy)
library(qlcMatrix) # For ultra-fast cosine
library(parallel)
setwd('/home/raquelcr/seurat')

# Load the Seurat object
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
# Select the first Seurat object for demonstration
seurat_obj <- seurat_list[[1]]
# Ensure PCA embeddings are available
if (!"pca" %in% names(seurat_obj@reductions)) {
    stop("PCA embeddings are not available in the Seurat object.")
}
# Get PCA embeddings (using first 20 PCs for balance)
pca_emb <- Embeddings(seurat_obj, "pca")[,1:20]
clusters <- Idents(seurat_obj)

calculate_cluster_metrics <- function(cluster_id) {
  cells <- which(clusters == cluster_id)
  if(length(cells) < 2) return(c(Euclidean=NA, Cosine=NA, Annoy=NA))
  
  cluster_emb <- pca_emb[cells, ]
  
  # 1. Euclidean Similarity
  euc_sim <- 1/(1 + mean(dist(cluster_emb)))
  
  # 2. Cosine Similarity (alternative robust calculation)
  norm_emb <- cluster_emb/sqrt(rowSums(cluster_emb^2))
  cos_sim <- mean(tcrossprod(norm_emb)[upper.tri(tcrossprod(norm_emb))])
  
  # 3. Annoy Similarity (fixed retrieval)
  annoy <- new(AnnoyAngular, ncol(cluster_emb))
  for(i in 1:nrow(cluster_emb)) annoy$addItem(i-1, cluster_emb[i,])
  annoy$build(50)
  
  annoy_sims <- sapply(1:nrow(cluster_emb), function(i) {
    neighbors <- annoy$getNNsByItem(i-1, 50)  # Returns indices only
    distances <- sapply(neighbors, function(j) annoy$getDistance(i-1, j))
    1 - mean(distances)/2  # Convert angular distance to cosine similarity
  })
  
  return(c(Euclidean=euc_sim, Cosine=cos_sim, Annoy=mean(annoy_sims)))
}

# Safe execution with error handling
metrics <- do.call(rbind, lapply(levels(clusters), function(cl) {
  tryCatch({
    calculate_cluster_metrics(cl)
  }, error = function(e) {
    message(sprintf("Error in cluster %s: %s", cl, e$message))
    c(Euclidean=NA, Cosine=NA, Annoy=NA)
  })
}))
rownames(metrics) <- levels(clusters)



# Calculate for all clusters
metrics <- t(sapply(levels(clusters), calculate_cluster_metrics))

metrics
