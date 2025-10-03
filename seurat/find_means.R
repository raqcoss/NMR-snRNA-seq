if (!require("signal", quietly = TRUE))
  install.packages("signal")
if (!require("pastecs", quietly = TRUE))
  install.packages("pastecs")

library(signal)   # For Savitzky-Golay smoothing
library(pastecs)  # For turnpoints()


find_maxs_mins_in_data <- function(y, x, polynomial_order = 2, window_size = 21, thresholds=NULL, num_min = NA, num_max = NA, saving_plot_path=NULL, metadata = NA, log_transform = FALSE) {
    # 1. Smooth the data (Savitzky-Golay filter)
    y_smooth <- sgolayfilt(y, p = polynomial_order, n = window_size)

    # 2. Find local extrema using turnpoints()
    tp <- turnpoints(y_smooth)
    extrema_indices <- which(tp$peaks | tp$pits)  # Indices of maxima (peaks) and minima (pits)

    # Classify as maxima or minima
    is_maxima <- tp$peaks[extrema_indices]
    is_minima <- tp$pits[extrema_indices]

    maxima_x <- x[extrema_indices[is_maxima]]
    maxima_y <- y_smooth[extrema_indices[is_maxima]]
    minima_x <- x[extrema_indices[is_minima]]
    minima_y <- y_smooth[extrema_indices[is_minima]]

    # 4. Sort and select the two lowest minima
    sorted_mins <- order(minima_y, decreasing = FALSE)
    if(identical(num_min, NA)){num_min <- length(sorted_mins)}
    num_min <- min(length(sorted_mins),num_min)
    top_mins <- sorted_mins[1:num_min]

    lowest_minima <- data.frame(
    x = minima_x[top_mins],
    y = minima_y[top_mins],
    type = rep('minimum',num_min)
    )
    # Sort and select highest maxima
    sorted_maxs <- order(maxima_y, decreasing = TRUE)
    if(identical(num_max, NA)){num_max <- length(sorted_maxs)}
    num_max <- min(num_max,length(sorted_maxs))
    top_maxs <- sorted_maxs[1:num_max]

    highest_maxima <- data.frame(
    x = maxima_x[top_maxs],
    y = maxima_y[top_maxs],
    type = rep('maximum',num_max)
    )
    maxs_n_mins <- rbind(highest_maxima, lowest_minima)

    # Plot the results
    if (!identical(saving_plot_path, NULL)){
      png(file=saving_plot_path,width=600, height=600)
      print(paste("Saving plot to:", saving_plot_path)) 
    }
    plot(x, y, col = "gray", main = paste(metadata, if(log_transform){'(Log2+1 transformed)'},"Smooth Maxima and Minima"), 
            xlab = paste0(if(log_transform){'Log_'},metadata), ylab = "nCells")
    lines(x, y_smooth, col = "blue", lwd = 2)
    points(maxima_x, maxima_y, col = "red", pch = 19, cex = 1.2)
    points(minima_x, minima_y, col = "green", pch = 19, cex = 1.2)
    points(lowest_minima$x, lowest_minima$y, col = "black", pch = 4, lwd = 3, cex = 2)
    points(highest_maxima$x, highest_maxima$y, col = '#4e4e4e', pch = 5, lwd = 3, cex = 2)

    legend_labels <- c("Frequency Data", "Smoothed", "Maxima", "Minima", "Top Minima", "Top Maxima")
    legend_cols   <- c("gray", "blue", "red", "green", "black", '#4e4e4e')
    legend_lty    <- c(NA, 1, NA, NA, NA, NA)
    legend_pch    <- c(1, NA, 19, 19, 4, 5)

    if (!is.null(thresholds)) {
      limits <- thresholds[,metadata]
      if(log_transform){limits <- log2(limits + 1)}
      print(paste("Adding vertical lines at:", paste(limits, collapse=", ")))  # Debug line
      abline(v=limits[1], col='blue', lty=2)
      abline(v=limits[2], col='gray', lty=4)
      abline(v=limits[3], col='red', lty=2)
      legend_labels <- c(legend_labels, "Lower cut", "No-tolerance Upper Cut", "Adjusted Upper cut")
      legend_cols   <- c(legend_cols, "blue", "gray", "red")
      legend_lty    <- c(legend_lty, 2, 4, 2)
      legend_pch    <- c(legend_pch, NA, NA, NA)
    }

    legend("topright", 
        legend = legend_labels,
        col = legend_cols,
        lty = legend_lty, 
        pch = legend_pch)

    if (!identical(saving_plot_path, NULL)){
      dev.off()
    }
    return(maxs_n_mins)
}


# Generate noisy data
#set.seed(40)
#x <- seq(0, 10, length.out = 200)
#y <- sin(x) + 0.5 * rnorm(length(x))  # Noisy sine wave
#mins <- find_maxs_mins_in_data(y, x, window_size = 101,polynomial_order = 4,num_max = 2)
#mins
max_n_mins_seurat <- function(seurat_obj, metadata = 'nCount_RNA', log_transform=FALSE, save_path = NULL, bins = 300,
        thresholds = NULL,
        polynomial_order = 2,
        window_size = 5,
        num_min = 5,
        num_max = 3){
  if(log_transform){
    data <- log2(seurat_obj@meta.data[,metadata] + 1)  # Adding 1 to avoid log(0)
  }
  else {
     data <- seurat_obj@meta.data[,metadata]
  }
  # Create histogram
  h <- hist(data, breaks = bins, plot = FALSE)
  # Use midpoints of bins as x values
  x <- h$mids
  y <- h$counts
  print(paste('x length:',length(x)))
  print(paste('y length:',length(y)))
  find_maxs_mins_in_data(y,
                        x, 
                        thresholds = thresholds,
                        polynomial_order = polynomial_order,
                        window_size = window_size,
                        num_min = num_min,
                        num_max = num_max,
                        saving_plot_path = save_path,
                        metadata = metadata,
                        log_transform = log_transform)
}