plot_qc_metrics <- function(
    seurat_list,
    metrics = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    thresholds = list(
      nFeature_RNA = c(200, 5000),
      nCount_RNA = c(150, 12000),
      percent.mt = c(0, 2)
    ),
    special_thresholds = NULL,
    file_name = 'qc_plot.jpeg',
    plot_width = NULL,
    plot_height = NULL,
    dpi = 500,
    point_alpha = 0.1,
    point_size = 0.5,
    show_points = FALSE,
    color_palette = c("#f1c40f", "#9b59b6"),
    theme_style = "classic") {
  
  # Load required packages quietly
  suppressPackageStartupMessages({
    require(ggplot2)
    require(patchwork)
    require(dplyr)
    require(tidyr)
    require(purrr)
  })
  
  ## Input validation
  if (!all(metrics %in% unique(unlist(map(seurat_list, ~names(.x@meta.data)))))) {
    stop("Some specified metrics not found in Seurat objects' metadata")
  }
  
  if (!all(metrics %in% names(thresholds))) {
    stop("Thresholds not provided for all metrics")
  }
  
  ## Create combined tidy data frame
  plot_data <- imap_dfr(seurat_list, function(obj, obj_name) {
    # Extract tissue and replicate from object name
    name_parts <- strsplit(obj_name, "_")[[1]]
    tissue <- paste(name_parts[-length(name_parts)], collapse = "_")
    replicate <- paste0("Rep", name_parts[length(name_parts)])
    
    # Get metadata for specified metrics
    obj@meta.data %>%
      select(all_of(metrics)) %>%
      mutate(
        Tissue = tissue,
        Replicate = replicate,
        Cell = rownames(obj@meta.data),
        .before = 1
      )
  }) %>%
    pivot_longer(
      cols = all_of(metrics),
      names_to = "Metric",
      values_to = "Value"
    )
  
  ## Create thresholds data frame
  thresholds_df <- map_dfr(metrics, function(m) {
    tibble(
      Metric = m,
      lower = thresholds[[m]][1],
      upper = if (m %in% names(special_thresholds)) {
        # Apply special thresholds if they exist
        map_dfr(unique(plot_data$Tissue), function(tissue) {
          if (tissue %in% names(special_thresholds[[m]])) {
            tibble(Tissue = tissue, upper = special_thresholds[[m]][[tissue]])
          } else {
            tibble(Tissue = tissue, upper = thresholds[[m]][2])
          }
        })
      } else {
        tibble(Tissue = unique(plot_data$Tissue), upper = thresholds[[m]][2])
      }
    ) %>% 
      unnest(upper)
  })
  
  ## Generate plots
  generate_plot <- function(metric) {
    p <- ggplot(
      plot_data %>% filter(Metric == metric),
      aes(x = Tissue, y = Value, fill = Replicate)
    ) +
      geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +
      scale_fill_manual(values = color_palette) +
      labs(y = metric, title = metric) +
      theme(
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12)
      )
    
    # Add points if requested
    if (show_points) {
      p <- p + geom_point(
        position = position_jitterdodge(
          jitter.width = 0.2,
          jitter.height = 0,
          dodge.width = 0.9  # Must match violin dodge
        ),
        alpha = point_alpha,
        size = point_size,
        shape = 16
      )
    }
    
    # Add threshold lines
    p <- p + 
      geom_hline(
        data = thresholds_df %>% filter(Metric == metric),
        aes(yintercept = lower),
        color = "blue",
        linetype = "dashed"
      ) +
      geom_hline(
        data = thresholds_df %>% filter(Metric == metric),
        aes(yintercept = upper),
        color = "red",
        linetype = "dashed"
      )
    
    # Apply selected theme
    switch(theme_style,
           "classic" = p + theme_classic(),
           "bw" = p + theme_bw(),
           "minimal" = p + theme_minimal(),
           p
    )
  }
  
  metric_plots <- map(metrics, generate_plot)
  
  ## Combine plots
  final_plot <- wrap_plots(metric_plots, ncol = 1) + 
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ## Set default plot dimensions if not specified
  if (is.null(plot_width)) {
    plot_width <- 5 * length(unique(plot_data$Tissue))
  }
  if (is.null(plot_height)) {
    plot_height <- 3 * length(metrics)
  }
  
  ## Save plot
  ggsave(
    filename = file_name,
    plot = final_plot,
    width = plot_width,
    height = plot_height,
    dpi = dpi,
    bg = "white"
  )
  
  return(invisible(final_plot))
}
##########################################################################

old_plot_qc_metrics <- function(file_name = 'qc_plot.jpeg', seurat_list, tissues, metrics, thresholds) {
  # Load required packages
  require(ggplot2)
  require(patchwork)
  
  metric_plots <- list()

  for (metric in metrics) {
    
    # Create plot for each tissue with both replicates
    tissue_plots <- lapply(tissues, function(tissue) {
      
      # Get both replicates for current tissue
      rep_data <- lapply(1:2, function(i) {
        idx <- which(grepl(paste0(tissue, "_"), names(seurat_list)))[i]
        seurat_list[[idx]]@meta.data[, metric, drop = FALSE]
      })
      
      # Combine into dataframe for plotting
      plot_data <- data.frame(
        Value = c(rep_data[[1]][,1], rep_data[[2]][,1]),
        Replicate = rep(c("Rep1", "Rep2"), 
                      times = c(nrow(rep_data[[1]]), nrow(rep_data[[2]]))),
        Tissue = tissue
      )
      
      # Determine upper threshold - special case for midbrain
      upper_threshold <- if(tissue == "midbrain" && metric == "percent.mt") {
        0.6
      } else {
        thresholds[[metric]][2]
      }
      
      # Create violin plot
      ggplot(plot_data, aes(x = Tissue, y = Value, fill = Replicate)) +
        geom_violin(trim = TRUE) +
        scale_fill_manual(values = c("#f1c40f", "#9b59b6")) + # Gold/Purple for reps
        theme_classic() +
        theme(
          axis.title.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 12)
        ) +
        ylab(metric) +
        ggtitle(tissue) +
        geom_hline(yintercept = thresholds[[metric]][1], color = "blue", linetype = "dashed") +
        geom_hline(yintercept = upper_threshold, color = "red", linetype = "dashed")
    })
    
    # Combine tissue plots for current metric
    metric_plots[[metric]] <- wrap_plots(tissue_plots, nrow = 1) + 
      plot_annotation(title = metric) &
      theme(plot.title = element_text(hjust = 0.5, size = 14))
  }

  # Combine all metrics vertically
  final_plot <- wrap_plots(
    metric_plots,
    ncol = 1,
    heights = c(1, 1, 1.1) # Slightly taller last plot for %MT
  ) + 
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  # Save with large canvas
  ggsave(
    file_name,
    final_plot,
    width = 16,  # Wider to accommodate side-by-side reps
    height = 12, # Taller for 3 metrics
    dpi = 500,
    bg = "white"
  )
  
  # Return the plot object in case user wants to modify or view it
  return(invisible(final_plot))
}

plot_pca_elbows <- function(seurat_list, file_name = "pca_elbow_plots.jpeg", ndims = 20, ncol = 2) {
  # Load required packages
  require(patchwork)
  require(ggplot2)
  
  # Input validation
  if (!inherits(seurat_list, "list") || !all(sapply(seurat_list, function(x) inherits(x, "Seurat")))) {
    stop("seurat_list must be a list of Seurat objects")
  }
  
  # Generate all plots
  plot_list <- lapply(names(seurat_list), function(obj_name) {
    # Verify the object has PCA computed
    if (!"pca" %in% names(seurat_list[[obj_name]]@reductions)) {
      warning(paste("PCA not computed for object:", obj_name, "- Skipping"))
      return(NULL)
    }
    
    ElbowPlot(seurat_list[[obj_name]], ndims = ndims) +
      ggtitle(obj_name) +
      theme(plot.title = element_text(size = 10, hjust = 0.5))
  })
  
  # Remove NULL plots (from skipped objects)
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  
  # Check if we have any plots left
  if (length(plot_list) == 0) {
    stop("No valid plots to combine - check if PCA was computed for any objects")
  }
  
  # Combine plots
  combined_plot <- wrap_plots(plot_list, ncol = ncol) + 
    plot_annotation(title = "Elbow Plots - PCA Dimensionality Reduction",
                   theme = theme(plot.title = element_text(hjust = 0.5, size = 14)))
  
  # Calculate dynamic height based on number of rows needed
  n_plots <- length(plot_list)
  nrow <- ceiling(n_plots / ncol)
  plot_height <- max(4, 3 * nrow)  # Minimum height 4, then 3 per row
  
  # Save plot
  ggsave(
    filename = file_name,
    plot = combined_plot,
    width = 8,
    height = plot_height,
    dpi = 300
  )
  
  message("Successfully saved combined elbow plot to: ", normalizePath(file_name))
  
  # Return the plot object invisibly for further inspection if needed
  invisible(combined_plot)
}

#### Umap of PCA clusters
plot_clusters_by_sample_at_resolutions <- function(seurat_list, resolutions, output_dir, reduction = 'umap', grid_cols = 3){
    for (i in seq_along(seurat_list)){
      obj <- seurat_list[[i]]
sample_umaps <- lapply(resolutions, function(res) {
  cluster_col <- paste0("RNA_snn_res.", res)
  DimPlot(obj, 
          reduction = reduction, 
          label = TRUE, 
          group.by = cluster_col) +
    ggtitle(paste("Resolution", res)) +
    theme(legend.position = "none")
})

combined_sample_umaps <- wrap_plots(sample_umaps, ncol = grid_cols) +
  plot_annotation(title = paste("Sample:", obj@project.name),
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))

ggsave(filename = paste0(output_dir, "/umap_all_res_",obj@project.name,".jpeg"),
       combined_sample_umaps, 
       width = 18, 
       height = 10, 
       dpi = 500)
}}

plot_clusters_at_resolutions <- function(obj, resolutions, output_dir, reduction = 'umap', grid_cols = 3){
require(Seurat)
sample_umaps <- lapply(resolutions, function(res){
  cluster_col <- paste0("RNA_snn_res.", res)
  DimPlot(obj, 
          reduction = reduction, 
          label = TRUE, 
          group.by = cluster_col) +
    ggtitle(paste("Resolution", res)) +
    theme(legend.position = "none")})

combined_sample_umaps <- wrap_plots(sample_umaps, ncol = grid_cols) +
  plot_annotation(title = paste("Sample:", obj@project.name),
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))

ggsave(filename = paste0(output_dir, "/umap_all_res_",obj@project.name,".jpeg"),
       combined_sample_umaps, 
       width = 18, 
       height = 10, 
       dpi = 500)
}



# Generate plots for all samples at this resolution
plot_clusters_by_res_for_all_samples <- function(seurat_list, resolutions, output_dir, reduction = 'umap'){
    for (res in resolutions) {
  umap_plots <- lapply(seurat_list, function(obj) {
    cluster_col <- paste0("RNA_snn_res.", res)
    DimPlot(obj, 
            reduction = reduction, 
            label = TRUE, 
            group.by = cluster_col) +
      ggtitle(paste0(obj@project.name, " (res=", res, ")")) +
      theme(legend.position = "none")
  })
  
  # Combine and save
  combined_umap <- wrap_plots(umap_plots, ncol = 2) + 
    plot_annotation(title = paste("Clustering at resolution", res),
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))
  
  ggsave(filename = paste0(output_dir, "/umap_res_", res, ".jpeg"),
         combined_umap, 
         width = 12, 
         height = 16, 
         dpi = 500)
  
  message("Saved UMAP plot for resolution ", res)
}
}
plot_clusters_at_chosen_res <- function(seurat_list, chosen_res, output_dir, reduction = 'umap') {
  # Input validation
  if (length(seurat_list) != length(chosen_res)) {
    stop("seurat_list and chosen_res must be the same length")
  }
  
  umap_plots <- list()
  
  for (i in seq_along(seurat_list)) {  # Fixed seq_along usage
    res <- chosen_res[[i]]
    cluster_col <- paste0("RNA_snn_res.", res)
    
    # Check if the cluster column exists
    if (!cluster_col %in% colnames(seurat_list[[i]]@meta.data)) {
      warning(paste("Cluster column", cluster_col, "not found in object", i))
      next
    }
    
    p <- DimPlot(seurat_list[[i]], 
                reduction = reduction, 
                label = TRUE, 
                group.by = cluster_col) +
         ggtitle(paste0(seurat_list[[i]]@project.name, " (res=", res, ")")) +
         theme(legend.position = "none")
    
    umap_plots[[i]] <- p  # Fixed list appending
  }
  
  # Remove NULL elements if any warnings occurred
  umap_plots <- umap_plots[!sapply(umap_plots, is.null)]
  
  if (length(umap_plots) == 0) {
    stop("No valid plots could be generated")
  }
  
  # Combine and save
  combined_umap <- wrap_plots(umap_plots, ncol = 2) +
    plot_annotation(title = "UMAPs at chosen resolutions",
                   theme = theme(plot.title = element_text(hjust = 0.5, size = 14)))
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ggsave(filename = file.path(output_dir, "umap_chosen_res.jpeg"),
         plot = combined_umap, 
         width = 5*2, 
         height = ceiling(length(umap_plots)/2) * 4,  # Dynamic height
         dpi = 500)
  
  message("Saved UMAP plot to: ", file.path(output_dir, "umap_chosen_res.jpeg"))
  
  return(invisible(combined_umap))
}

plot_annotated_clusters <- function(seurat_list, chosen_res, output_dir, 
                                   gsea_results_list,markers_list, reduction = 'umap',
                                   n_top_pathways = 3) {
  # Input validation
  if (length(seurat_list) != length(chosen_res)) {
    stop("seurat_list and chosen_res must be the same length")
  }
  umap_plots <- list()
  for (i in seq_along(seurat_list)) {
    res <- chosen_res[[i]]
    cluster_col <- paste0("RNA_snn_res.", res)
    obj <- seurat_list[[i]]
    # Get top enriched cell type pathways from GSEA
    top_celltypes <- gsea_results_list[[i]] %>%
    arrange(padj) %>%
    head(n_top_pathways) %>%
    mutate(celltype = gsub("^[^_]+_", "", pathway))  # Extract cell type from pathway name
  
    # Check if cluster column exists
    if (!cluster_col %in% colnames(obj@meta.data)) {
      warning(paste("Cluster column", cluster_col, "not found in object", i))
      next
    }
  
    # 2. Annotate clusters with top GSEA results
    cluster_annotations <- markers_list[[i]] %>%
      group_by(cluster) %>%
      arrange(desc(avg_log2FC)) %>%
      do({
        cluster_genes <- .$gene
        # Find best matching cell type
        matches <- top_celltypes %>%
          mutate(
            overlap = sapply(pathway, function(p) sum(cluster_genes %in% msigdb_sets[[p]])),
            coverage = overlap / length(msigdb_sets[[p]])
          ) %>%
          arrange(desc(coverage))
        
        data.frame(
          predicted_celltype = ifelse(nrow(matches) > 0, 
                                   paste0(matches$celltype[1], " (", matches$pathway[1], ")"), 
                                   "Unknown"),
          stringsAsFactors = FALSE
        )
      })
    
    # Add annotations to metadata
    obj$annotated_clusters <- plyr::mapvalues(
      obj@meta.data[[cluster_col]],
      from = cluster_annotations$cluster,
      to = cluster_annotations$predicted_celltype
    )
    
    # 3. Create annotated UMAP plot
    p <- DimPlot(obj, 
                reduction = reduction,
                group.by = "annotated_clusters",
                label = TRUE,
                repel = TRUE) +
        ggtitle(paste0(obj@project.name, 
                      "\nRes: ", res,
                      "\nTop cell types: ", 
                      paste(top_celltypes$celltype, collapse = ", "))) +
        theme(legend.position = "none",
              plot.title = element_text(size = 8))
    
    umap_plots[[i]] <- p
    seurat_list[[i]] <- obj  # Update object with annotations
  }
  
  # Combine plots
  combined_umap <- wrap_plots(umap_plots, ncol = 2) +
    plot_annotation(
      title = "UMAPs with GSEA-Annotated Cell Types",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 14))
    )
  
  # Save results
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ggsave(
    filename = file.path(output_dir, "annotated_umaps.jpeg"),
    plot = combined_umap,
    width = 12,
    height = ceiling(length(umap_plots)/2) * 5,
    dpi = 500
  )
  
  # Save metadata with annotations
  for (i in seq_along(seurat_list)) {
    write.csv(
      seurat_list[[i]]@meta.data,
      file.path(output_dir, paste0("metadata_", names(seurat_list)[i], ".csv"))
    )
  }
  
  message("Saved annotated UMAPs to: ", file.path(output_dir, "annotated_umaps.jpeg"))
  return(list(plots = combined_umap, seurat_list = seurat_list))
}

plot_pathways_grid <- function(seurat_list, chosen_res, output_dir, gsea_results_list, 
                             n_top_pathways = 3,width=20, height_per_row = 6) {
  
  # Load required packages
  require(ggplot2)
  require(patchwork)
  require(dplyr)
  require(Seurat)
  # Input validation
  stopifnot(
    length(seurat_list) == length(chosen_res),
    length(seurat_list) == length(gsea_results_list)
  )
  
  # Create annotation plots for all objects
  plot_list <- lapply(seq_along(seurat_list), function(i) {
    res <- chosen_res[[i]]
    obj <- seurat_list[[i]]
    cluster_col <- paste0("RNA_snn_res.", res)
    
    # Get top pathways for this object
    top_pathways <- gsea_results_list[[i]] %>%
      arrange(padj) %>%
      head(n_top_pathways) %>%
      mutate(celltype = gsub(".*_", "", pathway))  # Clean pathway names
    
    # Create annotation plot
    DimPlot(obj, 
           reduction = "umap",
           group.by = cluster_col,
           label = TRUE,
           repel = TRUE) +
      ggtitle(paste0(obj@project.name, 
                "\nResolution: ", res,
                "\nTop cell types: ", 
                paste(top_pathways$celltype, collapse = ", "))) +
      theme(plot.title = element_text(size = 9, hjust = 0.5),
            legend.position = "none")
  })
  
  # Create enrichment score plots for each sample
  enrichment_plots <- lapply(seq_along(gsea_results_list), function(i) {
    gsea_results_list[[i]] %>%
      arrange(padj) %>%
      head(n_top_pathways) %>% 
      ggplot(aes(x = reorder(pathway, NES), y = NES)) +
      geom_col(aes(fill = padj < 0.05)) +
      coord_flip() +
      labs(x = "Pathway", y = "NES", 
           title = paste("Top Pathways:", names(seurat_list)[i])) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 10))
  })
  
  # Combine all plots in a structured grid
n_rows <- length(seurat_list)
  combined_plot <- wrap_plots(
    c(plot_list, enrichment_plots),
    ncol = 2,
    heights = rep(height_per_row, n_rows)
  )  + 
    #plot_layout(heights = c(1, 1, 1, 1)) +  # Adjust heights as needed
    plot_annotation(title = "Integrated Single-Cell Analysis with Pathway Enrichment",
                   theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))
  # Save the combined plot
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ggsave(
    filename = file.path(output_dir, "combined_annotation_grid.png"),
    plot = combined_plot,
    width = width,
    height = height_per_row * n_rows,
    dpi = 300,
    bg = "white",
    limitsize = FALSE  # Override size checks
  )
  
  message("Saved combined plot to: ", file.path(output_dir, "combined_annotation_grid.png"))
  return(combined_plot)
}
