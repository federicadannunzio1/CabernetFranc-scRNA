##=============================================================================
## PHASE 6: Marker genes per cluster
##
## Input:  03_scRNA/04_analysis/VvProt_clustered.rds
## Output: 03_scRNA/04_analysis/markers/ (tables + plots)
##
## Uses RNA assay for DE (not SCT) as recommended by Seurat developers.
## Gene labels use Grapedia equivalence table when available.
##
## Run on TeraStat2: sbatch 06_markers.slurm
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
proj    <- "/lustre/home/dcapozzi/projects/scRNA_protoplasts"
indir   <- file.path(proj, "03_scRNA/04_analysis")
outdir  <- file.path(indir, "markers")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dpi <- 300

# Grapedia equivalence table
equiv_file <- file.path(proj, "01_reference/all_on_T2T_ref_equivalences_filtered_no_copies_symbols.csv")
# Functional annotation
annot_file <- file.path(proj, "01_reference/5.1_on_T2T_ref_functional_annotation_summary_with_coordinates.tsv")


##=============================================================================
## 1. LOAD CLUSTERED OBJECT
##=============================================================================

cat("=== 1. Loading clustered object ===\n")
obj <- readRDS(file.path(indir, "VvProt_clustered.rds"))
cat(sprintf("Cells: %d\n", ncol(obj)))
cat(sprintf("Genes: %d\n", nrow(obj)))
cat(sprintf("Clusters: %d\n", length(levels(Idents(obj)))))
cat(sprintf("Active ident: %s\n", levels(Idents(obj)) |> paste(collapse = ", ")))


##=============================================================================
## 2. LOAD GENE ANNOTATION
##=============================================================================

cat("\n=== 2. Loading gene annotations ===\n")

# Grapedia equivalence table
equiv <- tryCatch(
  read.delim(equiv_file, header = TRUE, stringsAsFactors = FALSE),
  error = function(e) { cat("WARNING: equivalence table not found.\n"); NULL }
)

# Functional annotation
annot <- tryCatch({
  a <- read.delim(annot_file, header = TRUE, stringsAsFactors = FALSE)
  a$gene_seurat <- gsub("_", "-", a$gene)
  a
}, error = function(e) { cat("WARNING: annotation file not found.\n"); NULL })

# Function to get gene label
get_gene_label <- function(gene_seurat, equiv, annot) {
  gid <- gsub("-", "_", gene_seurat)
  # Grapedia name first
  if (!is.null(equiv)) {
    idx <- match(gid, equiv$target_id)
    if (!is.na(idx) && !is.na(equiv$T2T.v5.1_gene_symbol[idx]) &&
        nchar(equiv$T2T.v5.1_gene_symbol[idx]) > 0) {
      return(equiv$T2T.v5.1_gene_symbol[idx])
    }
  }
  # Fallback: eggNOG
  if (!is.null(annot)) {
    idx <- match(gid, annot$gene)
    if (!is.na(idx) && !is.na(annot$eggNOG[idx]) && nchar(annot$eggNOG[idx]) > 0) {
      return(substr(annot$eggNOG[idx], 1, 40))
    }
  }
  return("")
}

cat("Annotations loaded.\n")


##=============================================================================
## 3. FIND ALL MARKERS
##=============================================================================

cat("\n=== 3. Finding markers (this may take a while) ===\n")

# Use RNA assay for DE
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, verbose = FALSE)

markers <- FindAllMarkers(obj,
                          only.pos = FALSE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25,
                          verbose = TRUE)

cat(sprintf("Total markers found: %d\n", nrow(markers)))
cat(sprintf("Markers per cluster:\n"))
marker_summary <- table(markers$cluster)
for (cl in names(marker_summary)) {
  cat(sprintf("  Cluster %s: %d markers\n", cl, marker_summary[cl]))
}


##=============================================================================
## 4. ADD GENE ANNOTATIONS TO MARKERS TABLE
##=============================================================================

cat("\n=== 4. Annotating markers ===\n")

markers$gene_label <- sapply(markers$gene, function(g) {
  get_gene_label(g, equiv, annot)
})
markers$display_name <- ifelse(nchar(markers$gene_label) > 0,
                               paste0(markers$gene, " (", markers$gene_label, ")"),
                               markers$gene)


##=============================================================================
## 5. TOP MARKERS TABLE
##=============================================================================

cat("\n=== 5. Top markers per cluster ===\n")

top10 <- markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(n = 10, order_by = avg_log2FC) |>
  dplyr::ungroup()

top5 <- markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(n = 5, order_by = avg_log2FC) |>
  dplyr::ungroup()

# Print top 5 per cluster
for (cl in levels(Idents(obj))) {
  cat(sprintf("\n--- Cluster %s ---\n", cl))
  cl_markers <- top5[top5$cluster == cl, ]
  for (i in seq_len(nrow(cl_markers))) {
    cat(sprintf("  %s  avg_log2FC=%.2f  pct.1=%.0f%%  pct.2=%.0f%%  p_adj=%.1e\n",
                cl_markers$display_name[i],
                cl_markers$avg_log2FC[i],
                cl_markers$pct.1[i] * 100,
                cl_markers$pct.2[i] * 100,
                cl_markers$p_val_adj[i]))
  }
}


##=============================================================================
## 6. PLOTS
##=============================================================================

cat("\n=== 6. Generating plots ===\n")

# Switch back to SCT for visualization
DefaultAssay(obj) <- "SCT"

# 6.1 Heatmap top 5 markers per cluster
top5_genes <- top5$gene
# Ensure unique genes
top5_genes <- unique(top5_genes)

p_heat <- DoHeatmap(obj, features = top5_genes, size = 3, angle = 90) +
  theme(axis.text.y = element_text(size = 5)) +
  NoLegend()
ggsave(file.path(outdir, "01_heatmap_top5_markers.png"), p_heat,
       width = 16, height = 14, dpi = dpi)

# 6.2 Dot plot top 3 markers per cluster
top3 <- markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(n = 3, order_by = avg_log2FC) |>
  dplyr::ungroup()

p_dot <- DotPlot(obj, features = unique(top3$gene)) +
  RotatedAxis() +
  labs(title = "Top 3 markers per cluster") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5))
ggsave(file.path(outdir, "02_dotplot_top3_markers.png"), p_dot,
       width = 18, height = 7, dpi = dpi)

# 6.3 Feature plots — top 2 markers per cluster (on UMAP)
for (cl in levels(Idents(obj))) {
  cl_top <- top5[top5$cluster == cl, ]
  genes_to_plot <- head(cl_top$gene, 2)
  if (length(genes_to_plot) > 0) {
    plots <- FeaturePlot(obj, features = genes_to_plot, pt.size = 0.3,
                         ncol = 2, combine = TRUE) &
      scale_color_viridis_c() & theme_pub()
    ggsave(file.path(outdir, sprintf("03_feature_cluster_%s.png", cl)), plots,
           width = 12, height = 5, dpi = dpi)
  }
}

# 6.4 Violin plots — top 2 markers per cluster
for (cl in levels(Idents(obj))) {
  cl_top <- top5[top5$cluster == cl, ]
  genes_to_plot <- head(cl_top$gene, 2)
  if (length(genes_to_plot) > 0) {
    plots <- VlnPlot(obj, features = genes_to_plot, pt.size = 0, ncol = 2) &
      theme_pub()
    ggsave(file.path(outdir, sprintf("04_violin_cluster_%s.png", cl)), plots,
           width = 14, height = 5, dpi = dpi)
  }
}


##=============================================================================
## 7. SAVE RESULTS
##=============================================================================

cat("\n=== 7. Saving results ===\n")

# CSV — all markers
write.csv(markers, file.path(outdir, "all_markers.csv"), row.names = FALSE)

# XLSX
if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  wb <- createWorkbook()
  headerStyle <- createStyle(textDecoration = "bold", border = "bottom")

  # Sheet 1: All markers
  addWorksheet(wb, "All markers")
  writeData(wb, "All markers", markers)
  addStyle(wb, "All markers", headerStyle, rows = 1, cols = 1:ncol(markers))
  setColWidths(wb, "All markers", cols = 1:ncol(markers), widths = "auto")

  # Sheet 2: Top 10 per cluster
  addWorksheet(wb, "Top 10 per cluster")
  writeData(wb, "Top 10 per cluster", top10)
  addStyle(wb, "Top 10 per cluster", headerStyle, rows = 1, cols = 1:ncol(top10))
  setColWidths(wb, "Top 10 per cluster", cols = 1:ncol(top10), widths = "auto")

  # Sheet 3: Summary (markers per cluster)
  summary_df <- data.frame(
    Cluster = names(marker_summary),
    N_markers = as.integer(marker_summary),
    Top_marker = sapply(names(marker_summary), function(cl) {
      m <- markers[markers$cluster == cl, ]
      if (nrow(m) > 0) m$display_name[which.max(m$avg_log2FC)] else ""
    }),
    Top_marker_log2FC = sapply(names(marker_summary), function(cl) {
      m <- markers[markers$cluster == cl, ]
      if (nrow(m) > 0) round(max(m$avg_log2FC), 2) else NA
    })
  )
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", summary_df)
  addStyle(wb, "Summary", headerStyle, rows = 1, cols = 1:ncol(summary_df))
  setColWidths(wb, "Summary", cols = 1:ncol(summary_df), widths = "auto")

  saveWorkbook(wb, file.path(outdir, "cluster_markers.xlsx"), overwrite = TRUE)
  cat("XLSX saved: cluster_markers.xlsx\n")
} else {
  cat("NOTE: openxlsx not available.\n")
}

cat("\n=== MARKER ANALYSIS COMPLETED ===\n")
cat(sprintf("Total markers: %d\n", nrow(markers)))
cat(sprintf("Clusters: %d\n", length(levels(Idents(obj)))))
cat("Results in:", outdir, "\n")
cat("\nNEXT STEP: Cell type annotation based on marker genes\n")
