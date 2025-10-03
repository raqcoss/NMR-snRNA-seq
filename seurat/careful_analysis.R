if (!require("Seurat", quietly = TRUE))
  install.packages("Seurat")
if (!require("harmony", quietly = TRUE))
  install.packages('harmony')
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager") 
install.packages('ggplot2')
if (!require('celldex',quietly = TRUE)){BiocManager::install("celldex")}
if (!require("SingleR", quietly = TRUE)){BiocManager::install("SingleR")}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
# update.packages(ask = FALSE, checkBuilt = TRUE)

# Analysing only for one sample
setwd('/home/raquelcr/seurat')
# Run one by one. Ex. Cerebral Cortex
tissues <- c("cerebral_cortex")#, 'hippocampus', 'midbrain') 
tissues <- c('hippocampus')#, 'midbrain') 

replicates <- c("1", "2")

# Load raw data

num <-2 # 0 starting from cerebral cortex, 2 from hippocampus and 4 from midbrain
seurat_list = list()
for (tissue in tissues){for (replica in replicates){
    num <- num + 1
    sample_path <- paste0('/home/raquelcr/cellranger_output/out_HetGla_female_1/NMR',num,'_',tissue,'/outs/raw_feature_bc_matrix')
    print(paste0('Loading sample in ' ,sample_path))
    sample_name <- paste0(tissue,'_',replica)
    data <- Read10X(data.dir = sample_path)
    obj <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3, min.features = 20)
    obj$tissue <- tissue
    obj$replicate <- replica
    seurat_list[[sample_name]] <- obj
    print(seurat_list[[sample_name]])
  }
}
# Drop any other column that is not data and counts
#seurat_list[[1]] <- DietSeurat(seurat_list[[1]], layers = c("counts", "data"))
#seurat_list[[2]] <- DietSeurat(seurat_list[[2]], layers = c("counts", "data"))

# Look out the count matrix
data[c("ADARB2", "RELN", "LIMS2"), 1:30] # Counts of some genes
rownames(seurat_list[[1]])[1:30] # gene list
colnames(seurat_list[[1]])[1:30] # cell barcodes
names(seurat_list) # sample names

# Merge data
merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[[2]],
  add.cell.ids = c("1", "2"),  # Prefix cell IDs to track origin
  merge.data = TRUE  # Ensure all data (counts, metadata) is merged
)

# Import mito genes from file
mito_genes <- scan("mito_genes.txt", what = character(), sep = "\n")
print(mito_genes)

# Calculate mitochondrial counts
mito_genes_in_assay <- mito_genes[mito_genes %in% rownames(merged_seurat)]
print(mito_genes_in_assay)
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, features = mito_genes_in_assay,assay = 'RNA')
colnames(merged_seurat@meta.data)
merged_seurat@meta.data$percent.mt[1:10]

#Calculate counts from Transcription Factors
tf_genes <- scan("mixed_tf_genes.txt", what = character(), sep = "\n")
tf_genes_in_assay <- tf_genes[tf_genes %in% rownames(merged_seurat)]
print(paste('TF genes in assay:', length(tf_genes_in_assay)))
layer_names <- Layers(merged_seurat, assay = "RNA")
print(layer_names)
# Initialize a named vector with all cell IDs
all_cells <- colnames(merged_seurat)
tf_counts <- setNames(numeric(length(all_cells)), all_cells)

# Sum TF counts across all layers
for (layer in layer_names) {
  # Get counts for current layer
  layer_counts <- LayerData(merged_seurat, assay = "RNA", layer = layer)
  layer_cells <- colnames(layer_counts)
  
  # Subset for TF genes that exist in this layer
  genes_in_layer <- tf_genes_in_assay[tf_genes_in_assay %in% rownames(layer_counts)]
  print(paste('TF genes in layer', layer, ':', length(genes_in_layer)))
  
  # Calculate TF sums only if genes are found
  if (length(genes_in_layer) > 0) {
    # Calculate TF sums for this layer
    layer_tf_sums <- Matrix::colSums(layer_counts[genes_in_layer, , drop = FALSE])
    
    # Add to the appropriate cells in the main vector
    tf_counts[layer_cells] <- tf_counts[layer_cells] + layer_tf_sums
  }
}

# Verify dimensions match
stopifnot(length(tf_counts) == ncol(merged_seurat))

# Add to metadata
merged_seurat$TF_counts_sum <- tf_counts[colnames(merged_seurat)]  # Ensure proper ordering

# Check results
head(merged_seurat$TF_counts_sum)
summary(merged_seurat$TF_counts_sum)


merged_seurat[["log10_nCount_RNA"]] <- log10(merged_seurat$nCount_RNA)
mean(na.omit(merged_seurat@meta.data$nCount_RNA))
source('ratopin_functions.R')
print(seurat_summary(merged_seurat))



# First annotation of outliers

# Fraction of desired cells = Expected cells / total of droplets
20000/length(merged_seurat$nCount_RNA)

# Vizualize 
#  1. Number of total counts per cell.
output_dir <- paste0("outs_", format(Sys.time(), "%Y%m%d_%H%M"))
{dir.create(output_dir)
print(paste0('Folder ',output_dir, ' was created'))}

png(file=paste0(output_dir,'/hist_raw_',tissue,'_nCount_RNA.png'),width=600, height=800)
hist(merged_seurat$nCount_RNA, breaks = 5000)
dev.off()

png(file=paste0(output_dir,'/hist_raw_',tissue,'_nCount_TF.png'),width=600, height=800)
hist(merged_seurat$TF_counts_sum, breaks = 5000)
dev.off()

# Visualize nConts and nFeatures scatter plot.
png(file=paste0(output_dir,'/scatter_raw_nCount_vs_nFeature_',tissue,'.png'),width=600, height=800)
plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
dev.off()

png(file=paste0(output_dir,'/scatter_raw_nCount_vs_mito_',tissue,'.png'),width=600, height=800)
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2
dev.off()

png(file=paste0(output_dir,'/scatter_raw_nCount_vs_TF_counts_',tissue,'.png'),width=600, height=800)
plotx <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "TF_counts_sum")
plotx
dev.off()
png(file=paste0(output_dir,'/scatter_raw_nFeature_vs_TF_counts_',tissue,'.png'),width=600, height=800)
plotx <- FeatureScatter(merged_seurat, feature1 = "nFeature_RNA", feature2 = "TF_counts_sum")
plotx
dev.off()
source('find_means.R')
max_n_mins_seurat(merged_seurat, metadata = 'TF_counts_sum', log_transform = TRUE,save_path = paste0(output_dir,'/TF_counts_sum_max_min_plot.png'), bins = 400, window_size = 91, num_min = 3, num_max = 4)
max_n_mins_seurat(merged_seurat, metadata = 'nCount_RNA', log_transform = TRUE,save_path = paste0(output_dir,'/nCounts_min_plot.png'), bins = 800, window_size = 201, num_min = 15, num_max = 6)
step = 0.001
quantiles = data.frame( quantile = seq(from=0, to=1, by=step),
  nCount_RNA = quantile(merged_seurat$nCount_RNA, probs=seq(from=0, to=1, by=step)),
  nFeature_RNA = quantile(merged_seurat$nFeature_RNA, probs=seq(from=0, to=1, by=step)),
  percent.mt = quantile(na.omit(merged_seurat$percent.mt), probs=seq(from=0, to=1, by=step))
)
expected_cells = 20000
window = ceiling(expected_cells/length(merged_seurat$nCount_RNA)/step)
window
### Please adjust the shift and tolerance prior to subsetting the data.
shift = 0.1 # 0 means no-left-shift, 1 means one window lef-shift
tolerance = min(0.1,shift) # 0 means no-tolerance, 1 means + 100% of window extension (`tolerance` must be <= `shift`)
# if `tolerance` == `shift`, doublets and some of most expressed cells are preserved
# if `tolerance` < `shift`, doublets and some of most expressed cells are removed
###

left_cut = 1/step-window*(1+shift)
orig_right_cut = 1/step-(window*shift)+1
right_cut = 1/step+(window*(tolerance-shift))+1
# Define prefiltering thresholds
prefilter_thesholds <- data.frame(
  nCount_RNA = c(quantiles$nCount_RNA[left_cut],quantiles$nCount_RNA[orig_right_cut], quantiles$nCount_RNA[right_cut]),
  nFeature_RNA = c(quantiles$nFeature[left_cut], quantiles$nFeature[orig_right_cut], quantiles$nFeature[right_cut]),
  percent.mt = c(quantiles$percent.mt[left_cut], quantiles$percent.mt[orig_right_cut], quantiles$percent.mt[right_cut])
)
prefilter_thesholds
max_n_mins_seurat(merged_seurat, metadata = 'nCount_RNA',thresholds = prefilter_thesholds, log_transform = TRUE,save_path = paste0(output_dir,'/cut_nCounts_max_min_plot.png'), bins = 800, window_size = 201, num_min = 5, num_max = 6)
# Visualize quantiles
png(file=paste0(output_dir,'/quantiles_raw_',tissue,'_nCount_RNA.png'),width=600, height=800)
plot(quantiles$quantile, quantiles$nCount, main = paste(tissue,"nCount per quantile"),
     xlab = "Quantile", ylab = "nCounts",
     pch = 19, frame = FALSE)
abline(v=quantiles$quantile[left_cut], lty=2,col="blue")
abline(v=quantiles$quantile[orig_right_cut], lty = 4, col="gray")
abline(v=quantiles$quantile[right_cut] , lty=2, col='red')
dev.off()
# Do not cut by upper threshold yet, because it is not clear if it is a doublet or not.
merged_seurat <- subset(merged_seurat, subset = 
             nCount_RNA > prefilter_thesholds$nCount_RNA[1] &
             # nCount_RNA < prefilter_thesholds$nCount_RNA[3] &
             nFeature_RNA > prefilter_thesholds$nFeature_RNA[1] &
             # nFeature_RNA < prefilter_thesholds$nFeature_RNA[3] &
             percent.mt < 5 # Default for scRNA-seq data
             )

max_n_mins_seurat(merged_seurat, metadata = 'nCount_RNA',thresholds= prefilter_thesholds,  log_transform = TRUE,save_path = paste0(output_dir,'/prefiltered2_nCounts_max_min_plot.png'), bins = 400, window_size = 101, num_min = 10, num_max = 6)
png(file=paste0(output_dir,'/prefiltered_QC_Violin plot_',tissue,'_.png'),width=800, height=800)
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha=0.2, pt.size = 0.1) +
  theme(plot.title = element_text(size = 10, hjust = 0.5))
dev.off()
# Check how many cells are left
length(merged_seurat$nCount_RNA)

# Save prefiltered object (Checkpoint)

save_seurat <- function(obj, name) {
  filename <- paste0(output_dir,"/", name, ".rds")
  saveRDS(obj, file = filename)
  return(filename)
}
save_seurat(merged_seurat, paste0(tissue,'_prefiltered'))


# Find doublets using scDblFinder (`scdblfinder_env` must be active during this step)
source('doublet_filtering.R') # go here and do the doublet filtering or run function
# Option 1. Run function
seurat_filtered <- run_doublet_filtering(merged_seurat, output_dir, save = TRUE)
# Option 2. Load prefiltered data without doublets
{
seurat_filtered = readRDS(paste0(output_dir,'/',tissue,'_without_dbl.rds'))
seurat_filtered$tissue <- tissue
seurat_filtered
}
# Vizualize prefiltered data
png(file=paste0(output_dir,'/dbl_filtered_QC_Violin plot_',tissue,'_.png'),width=800, height=800)
VlnPlot(seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha=0.2, pt.size = 0.1) +
  theme(plot.title = element_text(size = 10, hjust = 0.5))
dev.off()
source('find_means.R')
max_n_mins_seurat(seurat_filtered, metadata = 'nCount_RNA', log_transform = TRUE,save_path = paste0(output_dir,'/dbl_filtered_nCounts_max_min_plot.png'), bins = 400, window_size = 101, num_min = 10, num_max = 6)

# Perform second filtering
step = 0.005
quantiles = data.frame( quantile = seq(from=0, to=1, by=step),
  nCount_RNA = quantile(seurat_filtered$nCount_RNA, probs=seq(from=0, to=1, by=step)),
  nFeature_RNA = quantile(seurat_filtered$nFeature_RNA, probs=seq(from=0, to=1, by=step)),
  percent.mt = quantile(na.omit(seurat_filtered$percent.mt), probs=seq(from=0, to=1, by=step))
)
expected_cells = 20000
window = ceiling(expected_cells/length(seurat_filtered$nCount_RNA)/step)
window

### Please adjust the shift and tolerance prior to subsetting the data.
shift = 0.001 # 0 means no-left-shift, 1 means one window left-shift
tolerance = 0.001 # 0 means no-tolerance, 1 means + 100% of window extension
###

left_cut = 1/step-window*(1+shift)
orig_right_cut = 1/step-(window*shift)+1
right_cut = 1/step+(window*(tolerance-shift))+1
# Define prefiltering thresholds
filter_thesholds <- data.frame(
  nCount_RNA = c(quantiles$nCount_RNA[left_cut],quantiles$nCount_RNA[orig_right_cut], quantiles$nCount_RNA[right_cut]),
  nFeature_RNA = c(quantiles$nFeature[left_cut], quantiles$nFeature[orig_right_cut], quantiles$nFeature[right_cut]),
  percent.mt = c(quantiles$percent.mt[left_cut], quantiles$percent.mt[orig_right_cut], quantiles$percent.mt[right_cut])
)
filter_thesholds
max_n_mins_seurat(seurat_filtered, metadata = 'nCount_RNA',thresholds = filter_thesholds, log_transform = TRUE,save_path = paste0(output_dir,'/2nd_cut_nCounts_max_min_plot.png'), bins = 800, window_size = 201, num_min = 5, num_max = 6)
# Visualize quantiles
png(file=paste0(output_dir,'/quantiles_2nd_cut_',tissue,'_nCount_RNA.png'),width=600, height=800)
plot(quantiles$quantile, quantiles$nCount, main = paste(tissue,"nCount per quantile"),
     xlab = "Quantile", ylab = "nCounts",
     pch = 19, frame = FALSE)
abline(v=quantiles$quantile[left_cut], lty=2,col="blue")
abline(v=quantiles$quantile[orig_right_cut], lty = 4, col="gray")
abline(v=quantiles$quantile[right_cut] , lty=2, col='red')
dev.off()
# Do not cut by upper threshold if doublet filtering was done.
seurat_filtered <- subset(seurat_filtered, subset = 
             nCount_RNA > filter_thesholds$nCount_RNA[1] &
             # nCount_RNA < filter_thesholds$nCount_RNA[3] &
             nFeature_RNA > filter_thesholds$nFeature_RNA[1]
             # nFeature_RNA < filter_thesholds$nFeature_RNA[3]
             # & percent.mt < 5 # Default for scRNA-seq data
             )

max_n_mins_seurat(seurat_filtered, metadata = 'nCount_RNA',thresholds= filter_thesholds,  log_transform = TRUE,save_path = paste0(output_dir,'/filtered2_nCounts_max_min_plot.png'), bins = 400, window_size = 101, num_min = 10, num_max = 6)
png(file=paste0(output_dir,'/filtered2_QC_Violin plot_',tissue,'_.png'),width=800, height=800)
VlnPlot(seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha=0.2, pt.size = 0.1) +
  theme(plot.title = element_text(size = 10, hjust = 0.5))
dev.off()
# Check how many cells are left
length(seurat_filtered$nCount_RNA)
# Save prefiltered object (Checkpoint)

save_seurat <- function(obj, name) {
  filename <- paste0(output_dir,"/", name, ".rds")
  saveRDS(obj, file = filename)
  return(filename)
}
save_seurat(seurat_filtered, paste0(tissue,'_2nd_filtering'))

### You can continue with `clustering_data.R` script to perform clustering and further analysis.




###################################################################################################


# Visualize nConts and nFeatures scatter plot.
png(file=paste0(output_dir,'/scatter_prefiltered_nCount_vs_nFeature_',tissue,'.png'),width=600, height=800)
plot3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3
dev.off()

# Visualize nConts and nFeatures scatter plot.
png(file=paste0(output_dir,'/scatter_prefiltered_nCount_vs_scDblFinder.score_',tissue,'.png'),width=600, height=800)
plot4 <- FeatureScatter(seurat_obj, feature1 = "scDblFinder.score", feature2 = "nCount_RNA")
plot4
dev.off()

merged_seurat = seurat_obj
# Calculate quantiles and IQRs
step = 0.005
quantiles = data.frame( quantile = seq(from=0, to=1, by=step),
  nCount = quantile(merged_seurat$nCount_RNA, probs=seq(from=0, to=1, by=step)),
  nFeature = quantile(merged_seurat$nFeature_RNA, probs=seq(from=0, to=1, by=step)),
  percent.mt = quantile(na.omit(merged_seurat$percent.mt), probs=seq(from=0, to=1, by=step))
)
window = ceiling(20000/length(merged_seurat$nCount_RNA)/step)
tolerance = 1.05 # 1 means no-tolerance, 1.5 means + 50% of window extension
# Calculating Interquartile range
IQRs <- data.frame(metadata = c('nFeature', 'nCount', 'percent.mt'),
  Q1 = c(quantiles$nFeature[1/step/4+1], quantiles$nCount[1/step/4+1], quantiles$percent.mt[1/step/4+1]),
  Q3 = c(quantiles$nFeature[1/step/4*3+1], quantiles$nCount[1/step/4*3+1], quantiles$percent.mt[1/step/4*3+1])) 
IQRs$iqr = IQRs$Q3-IQRs$Q1
IQRs$upper_fence = IQRs$Q3+1.5*IQRs$iqr
log2(IQRs$upper_fence[1])
png(file=paste0(output_dir,'/scatter2_prefiltered_log2_nCount_vs_log2_nFeature_',tissue,'.png'),width=600, height=800)
p3 <- ggplot(merged_seurat@meta.data, aes(x = log2(nCount_RNA), y = log2(nFeature_RNA))) +
  geom_point() +
  geom_smooth(method = lm, color = "red", fill = "#69b3a2", se = TRUE) +
  labs(x = "log2(nCount_RNA)", y = "log2(nFeature_RNA)") +
  geom_hline(
    yintercept = log2(IQRs$upper_fence[1]),  # nFeature upper fence
    color = "darkgray", linetype = "dashed") +
  geom_vline(
    xintercept = log2(IQRs$upper_fence[2]),  # nCount upper fence
    color = "darkgray", linetype = "dashed")
  p3
dev.off()

png(file=paste0(output_dir,'/quantiles_prefiltered_',tissue,'_nCount_RNA.png'),width=600, height=800)
plot(quantiles$quantile, quantiles$nCount, main = paste(tissue,"nCount per quantile"),
     xlab = "Quantile", ylab = "nCounts",
     pch = 19, frame = FALSE)
abline(v=quantiles$quantile[left_cut], lty=2,col="blue")
abline(v=quantiles$quantile[orig_right_cut], col="gray")
abline(v=quantiles$quantile[right_cut], lty=2,col="red")
abline(v=quantiles$quantile[1/step+1], col='red')
abline(h=median(quantiles$nCount), col='black', lty = 3)
abline(h=Q1 - 1.5*IQR, col = 'darkgray', lty = 3)
abline(h=Q3 + 1.5*IQR, col = 'darkgray', lty = 3)
dev.off()

mean(quantiles$nCount)
quantiles
# Preliminary filtering
merged_seurat <- subset(merged_seurat, subset = 
             nCount_RNA > 7 &
             nFeature_RNA > 7 &
             percent.mt < 5
             #nFeature_RNA < thresholds[['nFeature_RNA']][2] & 
             #percent.mt > thresholds[['percent.mt']][1] 
             )
20000/length(merged_seurat$nCount_RNA)
step = 0.01
quantiles = data.frame( quantile = seq(from=0, to=1, by=step),
  nCount = quantile(merged_seurat$nCount_RNA, probs=seq(from=0, to=1, by=step)),
  nFeature = quantile(merged_seurat$nFeature_RNA, probs=seq(from=0, to=1, by=step)),
  percent.mt = quantile(na.omit(merged_seurat$percent.mt), probs=seq(from=0, to=1, by=step))
)
20000/length(merged_seurat$nCount_RNA)/step
# Choosing the range of cells that are between 95.5 to 99.5
valid_range <- tail(quantiles,20000/length(merged_seurat$nCount_RNA)/step+1)
valid_range
length(valid_range$quantile)
len <- length(valid_range$quantile)
valid_range[len-1,'nCount']
thresholds <- data.frame(Thresholds = c('Lower', 'Higher'),
  nCount = c(2^((log2(valid_range[1,'nCount'])+log2(valid_range[2,'nCount']))/2),2^((log2(valid_range[len,'nCount'])+log2(valid_range[len-1,'nCount']))/2)), 
  nFeature = c(2^((log2(valid_range[1,'nFeature'])+log2(valid_range[2,'nFeature']))/2),2^((log2(valid_range[len,'nFeature'])+log2(valid_range[len-1,'nFeature']))/2)),
  percent.mt = c(2^((log2(valid_range[1,'percent.mt'])+log2(valid_range[2,'percent.mt']))/2),2^((log2(valid_range[len,'percent.mt'])+log2(valid_range[len-1,'percent.mt']))/2))
)
thresholds

png(file = paste0(output_dir, '/scatter3_prefiltered_log2_nCount_vs_log2_nFeature_', tissue, '.png'), width = 600, height = 800)
p3 <- ggplot(merged_seurat@meta.data, aes(x = log2(nCount_RNA), y = log2(nFeature_RNA))) +
  geom_point() +
  labs(x = "log2(nCount_RNA)", y = "log2(nFeature_RNA)") +
  
  # Add lines and direct labels (no complex legend)
  geom_hline(yintercept = log2(IQRs$upper_fence[1]), 
             color = "darkgray", linetype = "dashed") +
  geom_vline(xintercept = log2(IQRs$upper_fence[2]), 
             color = "darkgray", linetype = "dashed") +
  geom_hline(yintercept = log2(thresholds$nFeature[2]), 
             color = "red", linetype = "dashed") +
  geom_vline(xintercept = log2(thresholds$nCount[2]), 
             color = "red", linetype = "dashed") +
  
  # Add text annotations instead of legend
  annotate("text", x = Inf, y = log2(IQRs$upper_fence[1]), 
           label = "IQR cutoff", color = "darkgray", hjust = 1.1, vjust = -0.5) +
  annotate("text", x = log2(IQRs$upper_fence[2]), y = Inf, 
           label = "IQR cutoff", color = "darkgray", hjust = 1.1, vjust = 1.5) +
  annotate("text", x = Inf, y = log2(thresholds$nFeature[2]), 
           label = "Quantile cutoff", color = "red", hjust = 1.1, vjust = -0.5) +
  annotate("text", x = log2(thresholds$nCount[2]), y = Inf, 
           label = "Quantile cutoff", color = "red", hjust = 1.1, vjust = 1.5)

print(p3)
dev.off()



png(file=paste0(output_dir,'/quantiles_',tissue,'_nCount_RNA.png'),width=600, height=800)
plot(quantiles$quantile, quantiles$nCount, main = paste(tissue,"nCount per quantile"),
     xlab = "Quantile", ylab = "nCounts",
     pch = 19, frame = FALSE)
abline(h=thresholds$nCount[1], col="blue")
abline(h=thresholds$nCount[2], col='red')
abline(v=valid_range$quantile[1]+step/2, col="blue")
abline(v=tail(valid_range,1)$quantile-step/2, col='red')
dev.off()

png(file=paste0(output_dir,'/hist_',tissue,'_nCount_RNA.png'),width=600, height=800)
hist(merged_seurat$nCount_RNA, breaks = 5000)
dev.off()
# Annotate outliers
merged_seurat <- subset(...)
dropout_mask <- merged_seurat$RNA_nCounts < thresholds$nCount[1]
doublet_mask <- merged_seurat$RNA_nCounts > thresholds$nCount[1]


merged_seurat[['dropout']] <- dropout_mask
merged_seurat[['doublet']] <- doublet_mask

# Filter cells according to found thresholds
merged_seurat <- subset(merged_seurat, subset = 
             nCount_RNA > thresholds$nCount[1] &
             nCount_RNA < thresholds$nCount[2] &
             nFeature_RNA > thresholds$nFeature[1] &
             nFeature_RNA < thresholds$nFeature[2]&
             percent.mt < thresholds$percent.mt[2]
             #nFeature_RNA < thresholds[['nFeature_RNA']][2] & 
             #percent.mt > thresholds[['percent.mt']][1] 
             )

#Look up data from selected cells
length(merged_seurat$nCount_RNA)
seurat_summary(merged_seurat)
write.csv(seurat_summary(merged_seurat),paste0(output_dir,'/',tissue,'_filtered_stats.csv'))

# Log2 transform
png(file=paste0(output_dir,'/hist_',tissue,'_log2_nCount_RNA.png'),width=600, height=800)
freqs = hist(log2(merged_seurat$nCount_RNA), breaks = 72)
freqs
dev.off()

source('find_means.R')
proposed_points <- find_maxs_mins_in_data(x=freqs$breaks[1:length(freqs$breaks)-1], y= freqs$counts, window_size = 5, saving_plot_path = paste0(output_dir,'/',tissue,'_max_min__nCounts_plot.png'))
proposed_points$original_x<- 2^proposed_points$x
proposed_points


#nFeatures
png(file=paste0(output_dir,'/hist_',tissue,'_log2_nFeatures_RNA.png'),width=600, height=800)
freqs = hist(log2(merged_seurat$nFeature_RNA), breaks = 32)
freqs
dev.off()

source('find_means.R')
proposed_points <- find_maxs_mins_in_data(x=freqs$breaks[1:length(freqs$breaks)-1], y= freqs$counts, window_size = 5, saving_plot_path = paste0(output_dir,'/',tissue,'_max_min_nFeatures_plot.png'))
proposed_points$original_x<- 2^proposed_points$x-1
proposed_points

# Calculate summary statistics
source('ratopin_functions.R')
summary_table <- seurat_summary(merged_seurat)
# Print the table
print(summary_table)

# Vizualization
source('/home/raquelcr/seurat/ploting_functions.R')
# Violin Plot
png(paste0(output_dir,'/QC_',tissue,'_violin_plot.png'),height = 1000, width = 600)
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha=0.5,)
dev.off()
# Scatter plot
plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(paste0(output_dir,'/QC_scatter_plot_',tissue,'.png'),height = 600, width = 1000)
plot1 + plot2
dev.off()

#Normalize
merged_seurat <- NormalizeData(merged_seurat)

#HVG's
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = 'vst', nfeatures= 2000)
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
resolutions <- c(0.2,0.25, 0.3, 0.35, 0.4, 0.5)
n_pc = 5
all_data <- merged_seurat %>%
  RunUMAP(dims = 1:n_pc) %>%
  FindNeighbors(dims = 1:n_pc) %>%
  FindClusters(resolution = resolutions)
source('/home/raquelcr/seurat/ploting_functions.R')
plot_clusters_at_resolutions(all_data,resolutions, output_dir)
# Add cluster column at chosen resolution
chosen_res <- 0.4
Idents(all_data) <- paste0("RNA_snn_res.", chosen_res)

#Print clusters
levels(Idents(all_data))
all_data

#Save object
save_seurat <- function(obj, name) {
  filename <- paste0(output_dir,"/", name, ".rds")
  saveRDS(obj, file = filename)
  return(filename)
}
save_seurat(all_data, paste0(tissue,'_clusters'))

