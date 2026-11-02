##=============================================================================
## PHASE 5: SCTransform, PCA, clustering, UMAP
##
## Input:  03_scRNA/03_filtering/VvProt_filtered.rds
## Output: 03_scRNA/04_analysis/VvProt_clustered.rds + plots
##
## Run on TeraStat2: sbatch 05_clustering.slurm
##=============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)

set.seed(42)

# Custom theme: black axes, no grid
theme_pub <- function() {
  theme_classic() %+replace%
    theme(axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.4))
}

# --- Paths ----
proj   <- "/lustre/home/dcapozzi/projects/scRNA_protoplasts"
indir  <- file.path(proj, "03_scRNA/03_filtering")
outdir <- file.path(proj, "03_scRNA/04_analysis")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dpi <- 300

# Stress gene list (from QC step)
stress_file <- file.path(proj, "03_scRNA/02_qc/gene_list_stress.txt")


##=============================================================================
## 1. LOAD FILTERED OBJECT
##=============================================================================

cat("=== 1. Loading filtered object ===\n")
obj <- readRDS(file.path(indir, "VvProt_filtered.rds"))
cat(sprintf("Cells: %d\n", ncol(obj)))
cat(sprintf("Genes: %d\n", nrow(obj)))


##=============================================================================
## 2. SCTRANSFORM NORMALIZATION
##=============================================================================

cat("\n=== 2. SCTransform normalization ===\n")
obj <- SCTransform(obj, verbose = TRUE, seed.use = 42)
cat("SCTransform done.\n")


##=============================================================================
## 3. CHECK VARIABLE FEATURES — stress gene dominance?
##=============================================================================

cat("\n=== 3. Checking variable features ===\n")

hvg <- VariableFeatures(obj)
cat(sprintf("Total variable features: %d\n", length(hvg)))

# Load stress gene list and check overlap
stress_genes <- c()
if (file.exists(stress_file)) {
  stress_genes <- readLines(stress_file)
  stress_genes <- stress_genes[stress_genes %in% rownames(obj)]
  stress_in_hvg <- intersect(hvg, stress_genes)
  cat(sprintf("Stress genes in object: %d\n", length(stress_genes)))
  cat(sprintf("Stress genes in HVG: %d / %d (%.1f%%)\n",
              length(stress_in_hvg), length(hvg),
              100 * length(stress_in_hvg) / length(hvg)))
  if (length(stress_in_hvg) > 0) {
    cat("Stress genes in HVG:\n")
    for (g in head(stress_in_hvg, 20)) cat(sprintf("  %s\n", g))
    if (length(stress_in_hvg) > 20) cat(sprintf("  ... and %d more\n", length(stress_in_hvg) - 20))
  }
} else {
  cat("WARNING: stress gene list not found, skipping check.\n")
}

# Variable feature plot
p_hvg <- VariableFeaturePlot(obj) +
  labs(title = sprintf("Variable features (n = %d)", length(hvg))) +
  theme_pub()
ggsave(file.path(outdir, "01_variable_features.png"), p_hvg,
       width = 8, height = 6, dpi = dpi)


##=============================================================================
## 4. PCA
##=============================================================================

cat("\n=== 4. Running PCA ===\n")
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

# Print top genes per PC
cat("\nTop genes per PC (1-5):\n")
print(obj[["pca"]], dims = 1:5, nfeatures = 5)

# Check stress gene contribution to PCs
if (length(stress_genes) > 0) {
  loadings <- Loadings(obj, reduction = "pca")
  cat("\nStress genes in top 50 loadings per PC:\n")
  for (pc in 1:10) {
    top_genes <- names(sort(abs(loadings[, pc]), decreasing = TRUE))[1:50]
    n_stress <- sum(top_genes %in% stress_genes)
    if (n_stress > 0) {
      cat(sprintf("  PC%d: %d/50 stress genes\n", pc, n_stress))
    }
  }
}

# Elbow plot
p_elbow <- ElbowPlot(obj, ndims = 50) +
  labs(title = "Elbow plot — standard deviation per PC") +
  theme_pub()
ggsave(file.path(outdir, "02_elbow_plot.png"), p_elbow,
       width = 8, height = 5, dpi = dpi)

# Variance explained
stdev <- Stdev(obj, reduction = "pca")
var_pct <- (stdev^2) / sum(stdev^2) * 100
cum_var <- cumsum(var_pct)

cat("\nVariance explained:\n")
for (i in c(10, 15, 20, 25, 30, 40, 50)) {
  cat(sprintf("  PC1-%d: %.1f%% cumulative variance\n", i, cum_var[i]))
}

# Select number of PCs — find elbow automatically
# Use the point where adding a PC contributes < 0.1% extra or < 5% of first PC
diff_var <- diff(var_pct)
n_pcs <- which(var_pct < 0.1 * var_pct[1] | var_pct < 0.1)[1] - 1
n_pcs <- max(n_pcs, 15)  # minimum 15
n_pcs <- min(n_pcs, 40)  # maximum 40
cat(sprintf("\nSelected PCs for clustering: %d\n", n_pcs))

# PCA scatter — PC1 vs PC2 colored by QC metrics
meta <- obj@meta.data
embed <- Embeddings(obj, "pca")[, 1:2]
pca_df <- cbind(as.data.frame(embed), meta)

pca_plots <- list()
for (feat in c("nCount_RNA", "nFeature_RNA", "percent.stress", "percent.chloro")) {
  pca_plots[[feat]] <- ggplot(pca_df, aes(x = PC_1, y = PC_2, color = .data[[feat]])) +
    geom_point(size = 0.3, alpha = 0.4) +
    scale_color_viridis_c() +
    labs(title = feat, x = "PC1", y = "PC2") +
    theme_pub() +
    theme(legend.key.width = unit(0.5, "cm"))
}
p_pca <- wrap_plots(pca_plots, ncol = 2)
ggsave(file.path(outdir, "03_PCA_qc_metrics.png"), p_pca,
       width = 12, height = 10, dpi = dpi)


##=============================================================================
## 5. KNN GRAPH + CLUSTERING (multiple resolutions)
##=============================================================================

cat("\n=== 5. Clustering ===\n")

obj <- FindNeighbors(obj, dims = 1:n_pcs, verbose = FALSE)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (res in resolutions) {
  obj <- FindClusters(obj, resolution = res, verbose = FALSE)
  n_clust <- length(levels(obj@meta.data[[paste0("SCT_snn_res.", res)]]))
  cat(sprintf("  Resolution %.1f: %d clusters\n", res, n_clust))
}

# Set default to 0.6 (can adjust later)
Idents(obj) <- "SCT_snn_res.0.6"
cat(sprintf("\nDefault clustering: resolution 0.6 (%d clusters)\n",
            length(levels(Idents(obj)))))


##=============================================================================
## 6. UMAP
##=============================================================================

cat("\n=== 6. Running UMAP ===\n")
obj <- RunUMAP(obj, dims = 1:n_pcs, verbose = FALSE, seed.use = 42)


##=============================================================================
## 7. PLOTS
##=============================================================================

cat("\n=== 7. Generating plots ===\n")

# 7.1 UMAP at multiple resolutions
umap_res_plots <- list()
for (res in resolutions) {
  col_name <- paste0("SCT_snn_res.", res)
  umap_res_plots[[as.character(res)]] <- DimPlot(obj, group.by = col_name, label = TRUE,
                                                  label.size = 3, pt.size = 0.3) +
    labs(title = sprintf("Resolution %.1f (%d clusters)",
                         res, length(unique(obj@meta.data[[col_name]])))) +
    theme_pub() + NoLegend()
}
p_res <- wrap_plots(umap_res_plots, ncol = 3)
ggsave(file.path(outdir, "04_UMAP_resolutions.png"), p_res,
       width = 18, height = 10, dpi = dpi)

# 7.2 UMAP default resolution (larger)
p_umap <- DimPlot(obj, label = TRUE, label.size = 4, pt.size = 0.5) +
  labs(title = sprintf("UMAP — %d clusters (resolution 0.6, %d PCs)",
                       length(levels(Idents(obj))), n_pcs)) +
  theme_pub()
ggsave(file.path(outdir, "05_UMAP_clusters.png"), p_umap,
       width = 9, height = 7, dpi = dpi)

# 7.3 UMAP colored by QC metrics
qc_feats <- c("nCount_RNA", "nFeature_RNA", "percent.stress",
               "percent.chloro", "percent.ribo", "log10GenesPerUMI")
qc_plots <- list()
for (feat in qc_feats) {
  qc_plots[[feat]] <- FeaturePlot(obj, features = feat, pt.size = 0.3) +
    scale_color_viridis_c() +
    labs(title = feat) +
    theme_pub()
}
p_qc <- wrap_plots(qc_plots, ncol = 3)
ggsave(file.path(outdir, "06_UMAP_qc_metrics.png"), p_qc,
       width = 15, height = 10, dpi = dpi)

# 7.4 UMAP colored by doublet score
p_dbl <- FeaturePlot(obj, features = "scDblFinder.score", pt.size = 0.3) +
  scale_color_viridis_c(option = "inferno") +
  labs(title = "Doublet score (scDblFinder)") +
  theme_pub()
ggsave(file.path(outdir, "07_UMAP_doublet_score.png"), p_dbl,
       width = 8, height = 7, dpi = dpi)

# 7.5 Cluster composition — cells per cluster
cluster_counts <- as.data.frame(table(Idents(obj)))
colnames(cluster_counts) <- c("Cluster", "Cells")
cluster_counts$Pct <- 100 * cluster_counts$Cells / sum(cluster_counts$Cells)

p_bar <- ggplot(cluster_counts, aes(x = Cluster, y = Cells, fill = Cluster)) +
  geom_col() +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Cells, Pct)),
            vjust = -0.2, size = 2.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = sprintf("Cells per cluster (total: %d)", ncol(obj)),
       x = "Cluster", y = "Number of cells") +
  theme_pub() + NoLegend()
ggsave(file.path(outdir, "08_cells_per_cluster.png"), p_bar,
       width = 10, height = 6, dpi = dpi)

# 7.6 Violin QC metrics per cluster
vln_plots <- list()
for (feat in c("nCount_RNA", "nFeature_RNA", "percent.stress", "percent.chloro")) {
  vln_plots[[feat]] <- VlnPlot(obj, features = feat, pt.size = 0) +
    labs(title = feat) +
    theme_pub() + NoLegend()
}
p_vln <- wrap_plots(vln_plots, ncol = 2)
ggsave(file.path(outdir, "09_violin_qc_per_cluster.png"), p_vln,
       width = 14, height = 10, dpi = dpi)


##=============================================================================
## 8. CLUSTER STATISTICS
##=============================================================================

cat("\n=== 8. Cluster statistics ===\n")

cat(sprintf("\n%-10s %8s %8s %10s %10s %10s %10s\n",
            "Cluster", "Cells", "Pct", "med_UMI", "med_genes", "med_stress", "med_chloro"))
cat(paste(rep("-", 72), collapse = ""), "\n")

for (cl in levels(Idents(obj))) {
  cells <- WhichCells(obj, idents = cl)
  m <- obj@meta.data[cells, ]
  cat(sprintf("%-10s %8d %7.1f%% %10.0f %10.0f %10.1f%% %10.1f%%\n",
              cl, length(cells),
              100 * length(cells) / ncol(obj),
              median(m$nCount_RNA),
              median(m$nFeature_RNA),
              median(m$percent.stress),
              median(m$percent.chloro)))
}

# Save cluster stats as CSV
cluster_stats <- do.call(rbind, lapply(levels(Idents(obj)), function(cl) {
  cells <- WhichCells(obj, idents = cl)
  m <- obj@meta.data[cells, ]
  data.frame(
    cluster = cl,
    n_cells = length(cells),
    pct_cells = 100 * length(cells) / ncol(obj),
    median_nCount = median(m$nCount_RNA),
    median_nFeature = median(m$nFeature_RNA),
    median_percent_stress = median(m$percent.stress),
    median_percent_chloro = median(m$percent.chloro),
    median_percent_ribo = median(m$percent.ribo)
  )
}))
write.csv(cluster_stats, file.path(outdir, "cluster_statistics.csv"), row.names = FALSE)


##=============================================================================
## 9. SAVE
##=============================================================================

cat("\n=== 9. Saving ===\n")
saveRDS(obj, file.path(outdir, "VvProt_clustered.rds"))
write.csv(obj@meta.data, file.path(outdir, "metadata_clustered.csv"), row.names = TRUE)

cat("\n=== CLUSTERING COMPLETED ===\n")
cat(sprintf("Object: VvProt_clustered.rds (%d cells, %d genes)\n", ncol(obj), nrow(obj)))
cat(sprintf("PCs used: %d\n", n_pcs))
cat(sprintf("Default resolution: 0.6 (%d clusters)\n", length(levels(Idents(obj)))))
cat("Plots in:", outdir, "\n")
cat("\nNEXT STEP: PHASE 6 — Marker genes + cell type annotation\n")
