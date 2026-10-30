##=============================================================================
## PHASE 3+4 (part 2): Filtering — apply QC thresholds and remove doublets
##
## Input:  03_scRNA/02_qc/VvProt_prefilter.rds (SoupX corrected, with doublet score)
## Output: 03_scRNA/03_filtering/VvProt_filtered.rds (ready for normalization)
##
## Thresholds determined after QC plot evaluation:
##   nCount_RNA:       [1313 — 57924]   (3 MAD su log10)
##   nFeature_RNA:     [748 — 8568]     (3 MAD su log10)
##   log10GenesPerUMI: > 0.80
##   scDblFinder:      class == "singlet"
##
## Run on TeraStat2: sbatch 04_filter.slurm
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
indir   <- file.path(proj, "03_scRNA/02_qc")
outdir  <- file.path(proj, "03_scRNA/03_filtering")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- Thresholds ----
NCOUNT_MIN   <- 1313
NCOUNT_MAX   <- 57924
NFEATURE_MIN <- 748
NFEATURE_MAX <- 8568
COMPLEXITY   <- 0.80


##=============================================================================
## 1. LOAD PRE-FILTERING OBJECT
##=============================================================================

cat("=== 1. Loading pre-filtering object ===\n")
obj <- readRDS(file.path(indir, "VvProt_prefilter.rds"))
cat(sprintf("Cells pre-filtering: %d\n", ncol(obj)))
cat(sprintf("Genes: %d\n", nrow(obj)))


##=============================================================================
## 2. APPLY FILTERS
##=============================================================================

cat("\n=== 2. Applying filters ===\n")

# Label each cell with exclusion reason
obj[["pass_nCount"]]   <- obj$nCount_RNA >= NCOUNT_MIN & obj$nCount_RNA <= NCOUNT_MAX
obj[["pass_nFeature"]] <- obj$nFeature_RNA >= NFEATURE_MIN & obj$nFeature_RNA <= NFEATURE_MAX
obj[["pass_complex"]]  <- obj$log10GenesPerUMI > COMPLEXITY
obj[["pass_singlet"]]  <- obj$scDblFinder.class == "singlet"
obj[["pass_all"]]      <- obj$pass_nCount & obj$pass_nFeature & obj$pass_complex & obj$pass_singlet

# Detailed report
cat(sprintf("  Excluded for nCount < %d:        %d\n", NCOUNT_MIN, sum(obj$nCount_RNA < NCOUNT_MIN)))
cat(sprintf("  Excluded for nCount > %d:       %d\n", NCOUNT_MAX, sum(obj$nCount_RNA > NCOUNT_MAX)))
cat(sprintf("  Excluded for nFeature < %d:       %d\n", NFEATURE_MIN, sum(obj$nFeature_RNA < NFEATURE_MIN)))
cat(sprintf("  Excluded for nFeature > %d:      %d\n", NFEATURE_MAX, sum(obj$nFeature_RNA > NFEATURE_MAX)))
cat(sprintf("  Excluded for complexity < %.2f:    %d\n", COMPLEXITY, sum(obj$log10GenesPerUMI <= COMPLEXITY)))
cat(sprintf("  Excluded as doublets:              %d\n", sum(obj$scDblFinder.class == "doublet")))
cat(sprintf("  Cells passing ALL filters:         %d / %d (%.1f%%)\n",
            sum(obj$pass_all), ncol(obj), 100 * mean(obj$pass_all)))

# Filter
obj_filtered <- subset(obj, pass_all == TRUE)
cat(sprintf("\nCells after filtering: %d\n", ncol(obj_filtered)))
cat(sprintf("Genes: %d\n", nrow(obj_filtered)))

# Remove genes expressed in 0 cells after filtering
counts_filt <- GetAssayData(obj_filtered, layer = "counts")
genes_keep <- rownames(counts_filt)[Matrix::rowSums(counts_filt > 0) >= 3]
obj_filtered <- subset(obj_filtered, features = genes_keep)
cat(sprintf("Genes after removal (min 3 cells): %d\n", nrow(obj_filtered)))


##=============================================================================
## 3. PRE/POST FILTERING PLOTS
##=============================================================================

cat("\n=== 3. Pre vs post filtering plots ===\n")

dpi <- 300

# 3.1 nCount density pre/post
pre_df  <- data.frame(nCount = obj$nCount_RNA, stato = "Pre-filtering")
post_df <- data.frame(nCount = obj_filtered$nCount_RNA, stato = "Post-filtering")
df_count <- rbind(pre_df, post_df)

p1 <- ggplot(df_count, aes(x = nCount, fill = stato)) +
  geom_density(alpha = 0.4) +
  scale_x_log10() +
  scale_fill_manual(values = c("Pre-filtering" = "grey60", "Post-filtering" = "steelblue")) +
  geom_vline(xintercept = c(NCOUNT_MIN, NCOUNT_MAX), linetype = "dashed", color = "red") +
  labs(title = "UMI distribution — pre vs post filtering",
       x = "nCount_RNA (log10)", y = "Density", fill = "") +
  theme_pub()
ggsave(file.path(outdir, "01_density_nCount_pre_post.png"), p1,
       width = 8, height = 5, dpi = dpi)

# 3.2 nFeature density pre/post
pre_df2  <- data.frame(nFeature = obj$nFeature_RNA, stato = "Pre-filtering")
post_df2 <- data.frame(nFeature = obj_filtered$nFeature_RNA, stato = "Post-filtering")
df_feat <- rbind(pre_df2, post_df2)

p2 <- ggplot(df_feat, aes(x = nFeature, fill = stato)) +
  geom_density(alpha = 0.4) +
  scale_x_log10() +
  scale_fill_manual(values = c("Pre-filtering" = "grey60", "Post-filtering" = "coral")) +
  geom_vline(xintercept = c(NFEATURE_MIN, NFEATURE_MAX), linetype = "dashed", color = "red") +
  labs(title = "Gene distribution — pre vs post filtering",
       x = "nFeature_RNA (log10)", y = "Density", fill = "") +
  theme_pub()
ggsave(file.path(outdir, "02_density_nFeature_pre_post.png"), p2,
       width = 8, height = 5, dpi = dpi)

# 3.3 Scatter post-filtering
p3 <- ggplot(obj_filtered@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(size = 0.2, alpha = 0.2, color = "steelblue") +
  scale_x_log10() + scale_y_log10() +
  labs(title = sprintf("nCount vs nFeature post-filtering (%d cells)", ncol(obj_filtered)),
       x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
  theme_pub()
ggsave(file.path(outdir, "03_scatter_post_filter.png"), p3,
       width = 7, height = 6, dpi = dpi)

# 3.4 Violin post-filtering
p4 <- VlnPlot(obj_filtered, features = c("nFeature_RNA", "nCount_RNA"),
              ncol = 2, pt.size = 0) & theme_pub() & NoLegend()
ggsave(file.path(outdir, "04_violin_post_filter.png"), p4,
       width = 10, height = 5, dpi = dpi)


##=============================================================================
## 4. POST-FILTERING STATISTICS
##=============================================================================

cat("\n=== POST-FILTERING STATISTICS ===\n")
cat(sprintf("Cells: %d (from %d, removed %.1f%%)\n",
            ncol(obj_filtered), ncol(obj),
            100 * (1 - ncol(obj_filtered) / ncol(obj))))
cat(sprintf("Genes: %d\n", nrow(obj_filtered)))

# Helper
has_meta <- function(obj, col) col %in% colnames(obj@meta.data)

# --- Summary table pre vs post filtering ---
calc_stats <- function(x) {
  c(min = min(x, na.rm = TRUE), Q1 = unname(quantile(x, 0.25, na.rm = TRUE)),
    median = median(x, na.rm = TRUE), Q3 = unname(quantile(x, 0.75, na.rm = TRUE)),
    max = max(x, na.rm = TRUE))
}

metrics <- c("nCount_RNA", "nFeature_RNA", "log10GenesPerUMI")
if (has_meta(obj, "percent.ribo"))   metrics <- c(metrics, "percent.ribo")
if (has_meta(obj, "percent.chloro")) metrics <- c(metrics, "percent.chloro")
if (has_meta(obj, "percent.stress")) metrics <- c(metrics, "percent.stress")

pre_rows  <- list()
post_rows <- list()
for (m in metrics) {
  pre_rows[[m]]  <- calc_stats(obj@meta.data[[m]])
  post_rows[[m]] <- calc_stats(obj_filtered@meta.data[[m]])
}

summary_table <- data.frame(
  metrica     = rep(metrics, each = 1),
  pre_min     = sapply(pre_rows, `[`, "min"),
  pre_Q1      = sapply(pre_rows, `[`, "Q1"),
  pre_median  = sapply(pre_rows, `[`, "median"),
  pre_Q3      = sapply(pre_rows, `[`, "Q3"),
  pre_max     = sapply(pre_rows, `[`, "max"),
  post_min    = sapply(post_rows, `[`, "min"),
  post_Q1     = sapply(post_rows, `[`, "Q1"),
  post_median = sapply(post_rows, `[`, "median"),
  post_Q3     = sapply(post_rows, `[`, "Q3"),
  post_max    = sapply(post_rows, `[`, "max"),
  row.names   = NULL
)

# Add row with cell and gene counts
counts_row <- data.frame(
  metrica = c("Cells", "Genes"),
  pre_min = c(ncol(obj), nrow(obj)),
  pre_Q1 = NA, pre_median = NA, pre_Q3 = NA, pre_max = NA,
  post_min = c(ncol(obj_filtered), nrow(obj_filtered)),
  post_Q1 = NA, post_median = NA, post_Q3 = NA, post_max = NA
)

summary_table <- rbind(counts_row, summary_table)

# Print to console
cat("\n=== SUMMARY TABLE PRE vs POST FILTERING ===\n\n")
cat(sprintf("  %-20s | %10s %10s %10s %10s %10s | %10s %10s %10s %10s %10s\n",
            "", "PRE_min", "PRE_Q1", "PRE_med", "PRE_Q3", "PRE_max",
            "POST_min", "POST_Q1", "POST_med", "POST_Q3", "POST_max"))
cat(paste(rep("-", 135), collapse = ""), "\n")
for (i in seq_len(nrow(summary_table))) {
  row <- summary_table[i, ]
  if (row$metrica %in% c("Cells", "Genes")) {
    cat(sprintf("  %-20s | %10.0f %10s %10s %10s %10s | %10.0f %10s %10s %10s %10s\n",
                row$metrica, row$pre_min, "", "", "", "",
                row$post_min, "", "", "", ""))
  } else {
    cat(sprintf("  %-20s | %10.1f %10.1f %10.1f %10.1f %10.1f | %10.1f %10.1f %10.1f %10.1f %10.1f\n",
                row$metrica,
                row$pre_min, row$pre_Q1, row$pre_median, row$pre_Q3, row$pre_max,
                row$post_min, row$post_Q1, row$post_median, row$post_Q3, row$post_max))
  }
}
cat("\n")

# Save table as CSV and XLSX
write.csv(summary_table, file.path(outdir, "summary_pre_vs_post.csv"), row.names = FALSE)

if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  wb <- createWorkbook()

  # Sheet 1: Pre vs post summary
  addWorksheet(wb, "Pre vs Post")
  writeData(wb, "Pre vs Post", summary_table)
  headerStyle <- createStyle(textDecoration = "bold", border = "bottom")
  addStyle(wb, "Pre vs Post", headerStyle, rows = 1, cols = 1:ncol(summary_table))
  setColWidths(wb, "Pre vs Post", cols = 1:ncol(summary_table), widths = "auto")

  # Sheet 2: Applied thresholds
  soglie_df <- data.frame(
    Filter = c("nCount_RNA", "nFeature_RNA", "log10GenesPerUMI", "scDblFinder", "SoupX"),
    Threshold_min = c(NCOUNT_MIN, NFEATURE_MIN, COMPLEXITY, NA, NA),
    Threshold_max = c(NCOUNT_MAX, NFEATURE_MAX, NA, NA, NA),
    Method = c("3 MAD (log10)", "3 MAD (log10)", "Standard", "scDblFinder", "SoupX v1.6.2"),
    Notes = c("UMI per cell", "Genes detected per cell",
              "Transcriptional complexity", "class == singlet, threshold 0.414",
              "Estimated contamination 2.9%")
  )
  addWorksheet(wb, "Applied thresholds")
  writeData(wb, "Applied thresholds", soglie_df)
  addStyle(wb, "Applied thresholds", headerStyle, rows = 1, cols = 1:ncol(soglie_df))
  setColWidths(wb, "Applied thresholds", cols = 1:ncol(soglie_df), widths = "auto")

  # Sheet 3: Filtering counts
  filter_counts <- data.frame(
    Description = c("Cells pre-filtering",
                    "Excluded for low nCount", "Excluded for high nCount",
                    "Excluded for low nFeature", "Excluded for high nFeature",
                    "Excluded for low complexity", "Excluded as doublets",
                    "Cells post-filtering", "Genes post-filtering"),
    Value = c(ncol(obj),
               sum(obj$nCount_RNA < NCOUNT_MIN), sum(obj$nCount_RNA > NCOUNT_MAX),
               sum(obj$nFeature_RNA < NFEATURE_MIN), sum(obj$nFeature_RNA > NFEATURE_MAX),
               sum(obj$log10GenesPerUMI <= COMPLEXITY),
               sum(obj$scDblFinder.class == "doublet"),
               ncol(obj_filtered), nrow(obj_filtered))
  )
  addWorksheet(wb, "Filtering counts")
  writeData(wb, "Filtering counts", filter_counts)
  addStyle(wb, "Filtering counts", headerStyle, rows = 1, cols = 1:ncol(filter_counts))
  setColWidths(wb, "Filtering counts", cols = 1:ncol(filter_counts), widths = "auto")

  saveWorkbook(wb, file.path(outdir, "QC_filtering_summary.xlsx"), overwrite = TRUE)
  cat("XLSX table saved: QC_filtering_summary.xlsx\n")
} else {
  cat("NOTE: openxlsx package not available, XLSX not generated.\n")
  cat("To install: install.packages('openxlsx')\n")
}


##=============================================================================
## 5. SAVE FILTERED OBJECT
##=============================================================================

# Remove temporary columns
obj_filtered$pass_nCount   <- NULL
obj_filtered$pass_nFeature <- NULL
obj_filtered$pass_complex  <- NULL
obj_filtered$pass_singlet  <- NULL
obj_filtered$pass_all      <- NULL

saveRDS(obj_filtered, file.path(outdir, "VvProt_filtered.rds"))
write.csv(obj_filtered@meta.data, file.path(outdir, "metadata_post_filter.csv"), row.names = TRUE)

cat("\n=== FILTERING COMPLETED ===\n")
cat("Filtered object: VvProt_filtered.rds\n")
cat("Metadata: metadata_post_filter.csv\n")
cat("Plots in:", outdir, "\n")
cat("\nNEXT STEP: PHASE 5 — SCTransform, PCA, clustering, UMAP\n")
