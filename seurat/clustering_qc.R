if(!require(proxy,quietly=TRUE)){install.packages("proxy", dependencies=TRUE)}
library(proxy)

# Use after clustering and cluster asignment

# Cosine Similarity
cosine_similarity_per_cluster <- function(data, clusters) {
  cluster_ids <- unique(clusters)
  results <- list()
  
  for (cid in cluster_ids) {
    cluster_cells <- data[, clusters == cid]
    if (ncol(cluster_cells) > 1) {
      sim_matrix <- as.matrix(proxy::dist(t(cluster_cells), method = "cosine"))
      diag(sim_matrix) <- NA  # remove self-similarity
      results[[as.character(cid)]] <- mean(sim_matrix, na.rm = TRUE)
    } else {
      results[[as.character(cid)]] <- NA
    }
  }
  return(results)
}

library(Matrix)
library(irlba)
library(parallel)

optimized_cosine_similarity <- function(data, clusters, n_pcs = 10, n_cores = 4) {
  # Convert to sparse matrix if not already
  if (!is(data, "sparseMatrix")) data <- as(data, "CsparseMatrix")
  
  # PCA reduction (memory-efficient)
  pca <- irlba::prcomp_irlba(t(data), n = n_pcs, center = TRUE, scale. = FALSE)
  reduced_data <- t(pca$x)  # Cells are now columns in PC space
  
  # Parallel computation by cluster
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("reduced_data"), envir = environment())
  
  results <- parLapply(cl, unique(clusters), function(cid) {
    cluster_cells <- reduced_data[, clusters == cid, drop = FALSE]
    n <- ncol(cluster_cells)
    if (n < 2) return(NA_real_)
    
    # Fast cosine similarity in reduced space
    norm_vec <- sqrt(colSums(cluster_cells^2))
    sim_matrix <- crossprod(cluster_cells) / (tcrossprod(norm_vec))
    mean(sim_matrix[upper.tri(sim_matrix)])  # Only compute upper triangle
  })
  
  stopCluster(cl)
  names(results) <- unique(clusters)
  return(results)
}
# Pearson Correlation
pearson_similarity_per_cluster <- function(data, clusters) {
  cluster_ids <- unique(clusters)
  results <- list()
  
  for (cid in cluster_ids) {
    cluster_cells <- data[, clusters == cid]
    if (ncol(cluster_cells) > 1) {
      cor_matrix <- cor(cluster_cells, method = "pearson")
      diag(cor_matrix) <- NA
      results[[as.character(cid)]] <- mean(cor_matrix, na.rm = TRUE)
    } else {
      results[[as.character(cid)]] <- NA
    }
  }
  
  return(results)
}

# Euclidean Similarity
euclidean_similarity_per_cluster <- function(data, clusters) {
  cluster_ids <- unique(clusters)
  results <- list()
  
  for (cid in cluster_ids) {
    cluster_cells <- data[, clusters == cid]
    if (ncol(cluster_cells) > 1) {
      dist_matrix <- as.matrix(dist(t(cluster_cells), method = "euclidean"))
      # Convert distance to similarity (inverse relationship)
      sim_matrix <- 1/(1 + dist_matrix)
      diag(sim_matrix) <- NA
      results[[as.character(cid)]] <- mean(sim_matrix, na.rm = TRUE)
    } else {
      results[[as.character(cid)]] <- NA
    }
  }
  
  return(results)
}

# Jaccard Similarity

jaccard_similarity_per_cluster <- function(data, clusters) {
  cluster_ids <- unique(clusters)
  results <- list()
  
  # Binarize expression (1 if expressed, 0 otherwise)
  binary_data <- data > 0
  
  for (cid in cluster_ids) {
    cluster_cells <- binary_data[, clusters == cid]
    if (ncol(cluster_cells) > 1) {
      sim_matrix <- proxy::simil(t(cluster_cells), method = "Jaccard")
      sim_matrix <- as.matrix(sim_matrix)
      diag(sim_matrix) <- NA
      results[[as.character(cid)]] <- mean(sim_matrix, na.rm = TRUE)
    } else {
      results[[as.character(cid)]] <- NA
    }
  }
  
  return(results)
}


# Approximate Nearest Neighbors Approach (Fastest for Very Large n)
library(RcppAnnoy)

annoy_similarity <- function(data, clusters, n_trees = 50, n_neighbors = 100) {
  f <- nrow(data)
  annoy <- new(AnnoyAngular, f)
  
  # Add items (ensure no NA/Inf values)
  for (i in 1:ncol(data)) {
    annoy$addItem(i-1, data[, i]) 
  }
  annoy$build(n_trees)
  
  results <- sapply(unique(clusters), function(cid) {
    cell_idx <- which(clusters == cid)
    if (length(cell_idx) < 2) return(NA_real_)
    
    sims <- sapply(cell_idx, function(i) {
      neighbors <- annoy$getNNsByItem(i-1, n_neighbors) + 1 # 1-based indexing
      neighbors <- intersect(neighbors, cell_idx) # Only within-cluster
      if (length(neighbors) < 2) return(NA_real_)
      
      # PROPER conversion to similarity
      distances <- sapply(neighbors-1, function(j) annoy$getDistance(i-1, j))
      mean(1 - (distances/2)) # Ensures ∈ [0,1]
    })
    
    mean(sims, na.rm = TRUE)
  })
  
  return(results)
}
