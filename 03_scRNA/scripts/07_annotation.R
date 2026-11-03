##=============================================================================
## PHASE 7: Cell type annotation — module scores + FindAllMarkers integration
##
## Input:  03_scRNA/04_analysis/VvProt_clustered.rds
##         03_scRNA/04_analysis/markers/all_markers.csv
## Output: 03_scRNA/04_analysis/annotation/ (plots + tables)
##
## Strategy:
##   A) Module scores from known leaf marker gene sets (Arabidopsis orthologs)
##   B) Overlap between FindAllMarkers DE genes and known marker gene sets
##   C) Combined evidence table: module score rank + DE overlap count + top DE
##      markers matching known sets → suggested cell type annotation
##
## Run on TeraStat2: sbatch 07_annotation.slurm
##=============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)

set.seed(42)

# Custom theme: black axes, no grid
theme_pub <- function() {
  theme_classic() %+replace%
    theme(axis.line  = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.4))
}

# --- Paths ----
proj   <- "/lustre/home/dcapozzi/projects/scRNA_protoplasts"
indir  <- file.path(proj, "03_scRNA/04_analysis")
outdir <- file.path(indir, "annotation")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dpi <- 300

# Reference files
equiv_file   <- file.path(proj, "01_reference/all_on_T2T_ref_equivalences_filtered_no_copies_symbols.csv")
annot_file   <- file.path(proj, "01_reference/5.1_on_T2T_ref_functional_annotation_summary_with_coordinates.tsv")
markers_file <- file.path(indir, "markers/all_markers.csv")


##=============================================================================
## 1. LOAD DATA
##=============================================================================

cat("=== 1. Loading data ===\n")
obj <- readRDS(file.path(indir, "VvProt_clustered.rds"))
DefaultAssay(obj) <- "SCT"

equiv <- read.delim(equiv_file, header = TRUE, stringsAsFactors = FALSE)
annot <- read.delim(annot_file, header = TRUE, stringsAsFactors = FALSE)
annot$gene_seurat <- gsub("_", "-", annot$gene)

# Load FindAllMarkers results
de_markers <- read.csv(markers_file, stringsAsFactors = FALSE)
cat(sprintf("Cells: %d, Clusters: %d\n", ncol(obj), length(levels(Idents(obj)))))
cat(sprintf("DE markers loaded: %d rows across %d clusters\n",
            nrow(de_markers), length(unique(de_markers$cluster))))


##=============================================================================
## 2. DEFINE MARKER GENE SETS FOR LEAF CELL TYPES
##=============================================================================

cat("\n=== 2. Defining marker gene sets ===\n")

# Function to find grapevine genes by Grapedia symbol pattern
find_by_symbol <- function(pattern, equiv, obj) {
  idx <- grep(pattern, equiv$T2T.v5.1_gene_symbol, ignore.case = TRUE)
  if (length(idx) == 0) return(character(0))
  genes_seurat <- gsub("_", "-", equiv$target_id[idx])
  genes_seurat[genes_seurat %in% rownames(obj)]
}

# Function to find grapevine genes by functional annotation keyword
find_by_annot <- function(pattern, annot, obj,
                          cols = c("Pfam", "PANTHER", "GO_biological_process",
                                   "GO_molecular_function", "eggNOG")) {
  hits <- character(0)
  for (col in cols) {
    if (col %in% colnames(annot)) {
      idx <- grep(pattern, annot[[col]], ignore.case = TRUE)
      hits <- union(hits, annot$gene_seurat[idx])
    }
  }
  hits[hits %in% rownames(obj)]
}

# --- Define cell type markers ---

# MESOPHYLL (palisade + spongy): photosynthesis, chloroplast
mesophyll_markers <- list(
  LHCB = find_by_symbol("^LHCB", equiv, obj),
  LHCA = find_by_symbol("^LHCA", equiv, obj),
  RBCS = find_by_symbol("^RBCS", equiv, obj),
  CAB  = find_by_symbol("^CAB[0-9]", equiv, obj),
  PSBS = find_by_symbol("^PSBS", equiv, obj),
  RCA  = find_by_symbol("^RCA$", equiv, obj),
  photosystem = find_by_annot("photosystem.II|oxygen.evolving|PSII", annot, obj)
)

# EPIDERMIS: cuticle, wax, lipid transfer
epidermis_markers <- list(
  PDF  = find_by_symbol("^PDF[12]", equiv, obj),
  CER  = find_by_symbol("^CER[0-9]", equiv, obj),
  KCS  = find_by_symbol("^KCS[0-9]", equiv, obj),
  LTP  = find_by_annot("lipid.transfer.protein|nsLTP", annot, obj),
  GDSL = find_by_annot("GDSL.*lipase|GDSL.*esterase", annot, obj),
  ATML1 = find_by_symbol("^ATML1", equiv, obj),
  FDH  = find_by_symbol("^FDH$", equiv, obj)
)

# GUARD CELLS: stomata
guard_markers <- list(
  SLAC1 = find_by_symbol("^SLAC[0-9]", equiv, obj),
  KAT   = find_by_symbol("^KAT[0-9]", equiv, obj),
  MYB60 = find_by_symbol("^MYB60", equiv, obj),
  FAMA  = find_by_symbol("^FAMA", equiv, obj),
  OST1  = find_by_symbol("^OST1", equiv, obj),
  slow_anion = find_by_annot("slow.anion.channel|SLAC", annot, obj),
  stomatal   = find_by_annot("stomat", annot, obj)
)

# VASCULAR - XYLEM
xylem_markers <- list(
  CESA_secondary = find_by_symbol("^CESA[478]", equiv, obj),
  IRX = find_by_symbol("^IRX", equiv, obj),
  VND = find_by_symbol("^VND[0-9]", equiv, obj),
  LAC = find_by_symbol("^LAC[0-9]", equiv, obj),
  XCP = find_by_symbol("^XCP[0-9]", equiv, obj),
  xylem = find_by_annot("xylan|secondary.cell.wall|lignin.biosyn", annot, obj)
)

# VASCULAR - PHLOEM
phloem_markers <- list(
  SEOR  = find_by_symbol("^SEOR", equiv, obj),
  PP2   = find_by_symbol("^PP2$|^PP2-", equiv, obj),
  CLE41 = find_by_symbol("^CLE41|^CLE44", equiv, obj),
  APL   = find_by_symbol("^APL$", equiv, obj),
  sieve = find_by_annot("sieve.element|phloem.protein|callose.synthase.*phloem", annot, obj)
)

# COMPANION CELLS
companion_markers <- list(
  SUC = find_by_symbol("^SUC[0-9]|^SUT[0-9]", equiv, obj),
  AHA3 = find_by_symbol("^AHA3", equiv, obj),
  sucrose_transp = find_by_annot("sucrose.*transport|sucrose.*symport", annot, obj)
)

# BUNDLE SHEATH
bundle_sheath_markers <- list(
  SCR = find_by_symbol("^SCR$|^SCARECROW", equiv, obj),
  SHR = find_by_symbol("^SHR$", equiv, obj),
  SCL = find_by_symbol("^SCL[0-9]", equiv, obj)
)

# DIVIDING CELLS (S-phase / proliferation)
dividing_markers <- list(
  histone_H3 = find_by_annot("histone.H3|histone.H4", annot, obj),
  PCNA   = find_by_symbol("^PCNA", equiv, obj),
  cyclin = find_by_annot("^cyclin[- ]|cyclin.dependent", annot, obj),
  MCM    = find_by_symbol("^MCM[0-9]", equiv, obj)
)

# Collect all cell type marker lists
cell_types <- list(
  Mesophyll        = unique(unlist(mesophyll_markers)),
  Epidermis        = unique(unlist(epidermis_markers)),
  Guard_cells      = unique(unlist(guard_markers)),
  Xylem            = unique(unlist(xylem_markers)),
  Phloem           = unique(unlist(phloem_markers)),
  Companion_cells  = unique(unlist(companion_markers)),
  Bundle_sheath    = unique(unlist(bundle_sheath_markers)),
  Dividing         = unique(unlist(dividing_markers))
)

# Report
cat("\nMarker genes found per cell type:\n")
for (ct in names(cell_types)) {
  cat(sprintf("  %-20s: %d genes\n", ct, length(cell_types[[ct]])))
}


##=============================================================================
## 3. MODULE SCORES PER CELL TYPE
##=============================================================================

cat("\n=== 3. Computing module scores ===\n")

# Switch to RNA assay: normalize and join layers for Seurat v5 compatibility
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- NormalizeData(obj, verbose = FALSE)
cat("RNA assay normalized for module score computation\n")

for (ct in names(cell_types)) {
  genes <- cell_types[[ct]]
  genes <- genes[genes %in% rownames(obj)]
  if (length(genes) >= 5) {
    obj <- AddModuleScore(obj, features = list(genes),
                          name = paste0("score_", ct), seed = 42)
    colnames(obj@meta.data)[ncol(obj@meta.data)] <- paste0("score_", ct)
    cat(sprintf("  %s: module score computed (%d genes)\n", ct, length(genes)))
  } else {
    cat(sprintf("  %s: skipped (only %d genes, need >= 5)\n", ct, length(genes)))
  }
}

##=============================================================================
## 3b. CELL CYCLE SCORING
##=============================================================================

cat("\n=== 3b. Cell cycle scoring ===\n")

# Define S-phase and G2M gene sets for grapevine (orthologs of Tirosh et al. 2016)
# RNA assay is already active and normalized

s_genes_vv <- unique(c(
  find_by_symbol("^MCM[2-7]", equiv, obj),
  find_by_symbol("^PCNA", equiv, obj),
  find_by_symbol("^RPA[0-9]", equiv, obj),
  find_by_symbol("^RFC[0-9]", equiv, obj),
  find_by_symbol("^POLD[0-9]", equiv, obj),
  find_by_symbol("^PRIM[0-9]", equiv, obj),
  find_by_symbol("^FEN1", equiv, obj),
  find_by_symbol("^RRM[12]", equiv, obj),
  find_by_symbol("^RAD51$", equiv, obj),
  find_by_symbol("^NASP", equiv, obj),
  find_by_annot("DNA.replication|DNA.polymerase.*alpha|DNA.primase", annot, obj)
))

g2m_genes_vv <- unique(c(
  find_by_symbol("^CYCB[0-9]", equiv, obj),
  find_by_symbol("^CDKB[0-9]", equiv, obj),
  find_by_symbol("^CDC20", equiv, obj),
  find_by_symbol("^AURORA", equiv, obj),
  find_by_symbol("^TOP2", equiv, obj),
  find_by_symbol("^NACK[12]", equiv, obj),
  find_by_symbol("^KN$|^KNOLLE", equiv, obj),
  find_by_annot("condensin|kinesin.*mitot|septin|cytokinesis", annot, obj)
))

cat(sprintf("S-phase genes found: %d\n", length(s_genes_vv)))
cat(sprintf("G2M genes found: %d\n", length(g2m_genes_vv)))

if (length(s_genes_vv) >= 5 && length(g2m_genes_vv) >= 5) {
  obj <- CellCycleScoring(obj,
                          s.features   = s_genes_vv,
                          g2m.features = g2m_genes_vv,
                          set.ident    = FALSE)
  cat(sprintf("\nCell cycle phase distribution:\n"))
  phase_tab <- table(obj$Phase)
  for (ph in names(phase_tab)) {
    cat(sprintf("  %s: %d cells (%.1f%%)\n", ph, phase_tab[ph],
                100 * phase_tab[ph] / ncol(obj)))
  }

  # Per-cluster cell cycle summary
  cc_per_cluster <- do.call(rbind, lapply(levels(Idents(obj)), function(cl) {
    cells <- WhichCells(obj, idents = cl)
    ph <- obj$Phase[cells]
    n <- length(cells)
    data.frame(
      Cluster  = cl,
      N_cells  = n,
      Pct_G1   = round(100 * sum(ph == "G1") / n, 1),
      Pct_S    = round(100 * sum(ph == "S") / n, 1),
      Pct_G2M  = round(100 * sum(ph == "G2M") / n, 1),
      Pct_cycling = round(100 * sum(ph != "G1") / n, 1)
    )
  }))
  cat("\nCell cycle per cluster:\n")
  print(cc_per_cluster, row.names = FALSE)
} else {
  cat("WARNING: Not enough S/G2M genes found for cell cycle scoring\n")
  obj$Phase <- "Unknown"
  obj$S.Score <- 0
  obj$G2M.Score <- 0
  cc_per_cluster <- NULL
}

# Switch back to SCT for visualization
DefaultAssay(obj) <- "SCT"


##=============================================================================
## 4. DE MARKERS vs KNOWN MARKER GENE SETS — OVERLAP ANALYSIS
##=============================================================================

cat("\n=== 4. Overlap between DE markers and known gene sets ===\n")

# Keep only significant upregulated markers for overlap (adjusted p < 0.05, logFC > 0)
de_up <- de_markers[de_markers$p_val_adj < 0.05 & de_markers$avg_log2FC > 0, ]
cat(sprintf("Significant upregulated DE markers: %d\n", nrow(de_up)))

clusters <- sort(unique(de_up$cluster))

# For each cluster x cell type: count overlap and list matching genes
overlap_list <- list()

for (cl in clusters) {
  cl_genes <- de_up$gene[de_up$cluster == cl]
  row <- list(cluster = cl)
  for (ct in names(cell_types)) {
    ct_genes <- cell_types[[ct]]
    hits <- intersect(cl_genes, ct_genes)
    row[[paste0(ct, "_n")]]    <- length(hits)
    row[[paste0(ct, "_genes")]] <- paste(hits, collapse = ";")
  }
  overlap_list[[length(overlap_list) + 1]] <- row
}

# Build overlap count matrix (clusters x cell types)
overlap_mat <- do.call(rbind, lapply(overlap_list, function(r) {
  counts <- sapply(names(cell_types), function(ct) r[[paste0(ct, "_n")]])
  c(cluster = r$cluster, counts)
}))
overlap_mat <- as.data.frame(overlap_mat, stringsAsFactors = FALSE)
for (ct in names(cell_types)) overlap_mat[[ct]] <- as.integer(overlap_mat[[ct]])

cat("\nDE marker overlap with known gene sets (count):\n")
print(overlap_mat, row.names = FALSE)


##=============================================================================
## 5. COMBINED ANNOTATION TABLE
##=============================================================================

cat("\n=== 5. Building combined annotation ===\n")

# Gene label helper
get_gene_label <- function(gene_seurat, equiv) {
  gid <- gsub("-", "_", gene_seurat)
  idx <- match(gid, equiv$target_id)
  if (!is.na(idx) && !is.na(equiv$T2T.v5.1_gene_symbol[idx]) &&
      nchar(equiv$T2T.v5.1_gene_symbol[idx]) > 0) {
    return(equiv$T2T.v5.1_gene_symbol[idx])
  }
  return(gene_seurat)
}

score_cols <- grep("^score_", colnames(obj@meta.data), value = TRUE)

combined_annot <- do.call(rbind, lapply(levels(Idents(obj)), function(cl) {
  cells <- WhichCells(obj, idents = cl)
  m <- obj@meta.data[cells, ]

  # --- A: Module scores ---
  scores <- sapply(score_cols, function(sc) mean(m[[sc]], na.rm = TRUE))
  names(scores) <- gsub("score_", "", names(scores))
  score_rank <- names(sort(scores, decreasing = TRUE))
  best_score_ct <- score_rank[1]
  best_score_val <- scores[best_score_ct]

  # --- B: DE overlap counts ---
  ol_row <- overlap_list[[which(sapply(overlap_list, function(r) r$cluster == cl))]]
  ol_counts <- sapply(names(cell_types), function(ct) ol_row[[paste0(ct, "_n")]])
  best_overlap_ct <- names(which.max(ol_counts))
  best_overlap_n  <- max(ol_counts)

  # Matching gene names for best overlap
  best_overlap_genes <- ol_row[[paste0(best_overlap_ct, "_genes")]]
  # Get Grapedia labels for those genes
  if (nchar(best_overlap_genes) > 0) {
    gene_ids <- strsplit(best_overlap_genes, ";")[[1]]
    gene_labs <- sapply(gene_ids[seq_len(min(5, length(gene_ids)))],
                        function(g) get_gene_label(g, equiv))
    top_overlap_labels <- paste(gene_labs, collapse = ", ")
  } else {
    top_overlap_labels <- ""
  }

  # --- C: Cell cycle info (informational only, no correction applied) ---
  pct_cycling <- 0
  pct_s <- 0
  pct_g2m <- 0
  if (!is.null(cc_per_cluster)) {
    cc_row <- cc_per_cluster[cc_per_cluster$Cluster == cl, ]
    pct_cycling <- cc_row$Pct_cycling
    pct_s   <- cc_row$Pct_S
    pct_g2m <- cc_row$Pct_G2M
  }

  # --- D: Combined verdict (module score + DE overlap) ---
  if (best_score_ct == best_overlap_ct) {
    suggested <- best_score_ct
    confidence <- "high"
  } else if (best_overlap_n == 0) {
    suggested <- best_score_ct
    confidence <- "low (no DE overlap)"
  } else {
    suggested <- paste0(best_score_ct, " / ", best_overlap_ct)
    confidence <- "medium (methods disagree)"
  }

  # Top 5 positive DE markers for this cluster (any cell type)
  cl_top <- de_up[de_up$cluster == cl, ]
  cl_top <- cl_top[order(cl_top$avg_log2FC, decreasing = TRUE), ]
  top5 <- head(cl_top, 5)
  top5_labels <- sapply(top5$gene, function(g) get_gene_label(g, equiv))
  top5_str <- paste(top5_labels, collapse = ", ")

  data.frame(
    Cluster          = cl,
    N_cells          = length(cells),
    Pct_cells        = round(100 * length(cells) / ncol(obj), 1),
    ModuleScore_best = best_score_ct,
    ModuleScore_val  = round(best_score_val, 4),
    DEoverlap_best   = best_overlap_ct,
    DEoverlap_n      = best_overlap_n,
    DEoverlap_genes  = top_overlap_labels,
    Pct_cycling      = pct_cycling,
    Pct_S            = pct_s,
    Pct_G2M          = pct_g2m,
    Suggested_type   = suggested,
    Confidence       = confidence,
    Top5_DE_markers  = top5_str,
    Median_stress    = round(median(m$percent.stress), 1),
    Median_chloro    = round(median(m$percent.chloro), 1),
    stringsAsFactors = FALSE
  )
}))

# Print summary
cat("\n=== COMBINED ANNOTATION SUMMARY ===\n\n")
cat(sprintf("%-5s %5s %5s  %-16s %7s  %-16s %3s  %6s %4s %4s  %-24s  %-28s  %6s %6s\n",
            "Cl", "Cells", "Pct", "ModScore", "Val",
            "DEoverlap", "N", "Cycl%", "S%", "G2M%",
            "Suggested", "Confidence", "Str%", "Chl%"))
cat(paste(rep("-", 160), collapse = ""), "\n")
for (i in seq_len(nrow(combined_annot))) {
  r <- combined_annot[i, ]
  cat(sprintf("%-5s %5d %4.1f%%  %-16s %7.4f  %-16s %3d  %5.1f%% %4.1f %4.1f  %-24s  %-28s  %5.1f%% %5.1f%%\n",
              r$Cluster, r$N_cells, r$Pct_cells,
              r$ModuleScore_best, r$ModuleScore_val,
              r$DEoverlap_best, r$DEoverlap_n,
              r$Pct_cycling, r$Pct_S, r$Pct_G2M,
              r$Suggested_type, r$Confidence,
              r$Median_stress, r$Median_chloro))
}


##=============================================================================
## 6. PLOTS — MODULE SCORES ON UMAP
##=============================================================================

cat("\n=== 6. Plotting module scores on UMAP ===\n")

score_plots <- list()
for (sc in score_cols) {
  ct_name <- gsub("score_", "", sc)
  ct_name <- gsub("_", " ", ct_name)
  score_plots[[sc]] <- FeaturePlot(obj, features = sc, pt.size = 0.3) +
    scale_color_viridis_c(option = "magma") +
    labs(title = ct_name) +
    theme_pub()
}
p_scores <- wrap_plots(score_plots, ncol = 3)
ggsave(file.path(outdir, "01_UMAP_module_scores.png"), p_scores,
       width = 16, height = ceiling(length(score_plots) / 3) * 5, dpi = dpi)


##=============================================================================
## 7. PLOTS — VIOLIN MODULE SCORES PER CLUSTER
##=============================================================================

cat("\n=== 7. Violin module scores per cluster ===\n")

vln_score_plots <- list()
for (sc in score_cols) {
  ct_name <- gsub("score_", "", sc)
  ct_name <- gsub("_", " ", ct_name)
  vln_score_plots[[sc]] <- VlnPlot(obj, features = sc, pt.size = 0) +
    labs(title = ct_name) +
    theme_pub() + NoLegend()
}
p_vln_scores <- wrap_plots(vln_score_plots, ncol = 2)
ggsave(file.path(outdir, "02_violin_module_scores.png"), p_vln_scores,
       width = 14, height = ceiling(length(vln_score_plots) / 2) * 4, dpi = dpi)


##=============================================================================
## 8. PLOTS — DE OVERLAP HEATMAP
##=============================================================================

cat("\n=== 8. DE overlap heatmap ===\n")

# Build matrix for heatmap
heat_data <- overlap_mat[, names(cell_types)]
rownames(heat_data) <- overlap_mat$cluster

# Reshape for ggplot
heat_long <- data.frame(
  Cluster   = rep(overlap_mat$cluster, times = length(cell_types)),
  Cell_type = rep(names(cell_types), each = nrow(overlap_mat)),
  Count     = as.integer(unlist(heat_data)),
  stringsAsFactors = FALSE
)
heat_long$Cluster <- factor(heat_long$Cluster,
                            levels = as.character(sort(as.integer(unique(heat_long$Cluster)))))
heat_long$Cell_type <- factor(heat_long$Cell_type, levels = rev(names(cell_types)))

p_heat <- ggplot(heat_long, aes(x = Cluster, y = Cell_type, fill = Count)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(Count > 0, Count, "")), size = 3) +
  scale_fill_gradient(low = "white", high = "#B2182B", name = "N markers") +
  labs(title = "DE marker overlap with known cell type gene sets",
       x = "Cluster", y = "") +
  theme_pub() +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5))
ggsave(file.path(outdir, "03_heatmap_DE_overlap.png"), p_heat,
       width = 10, height = 6, dpi = dpi)


##=============================================================================
## 9. DOT PLOT — KEY KNOWN MARKERS FOUND IN DE RESULTS
##=============================================================================

cat("\n=== 9. Dot plot: known markers found in DE ===\n")

# For each cell type, find known marker genes that are also DE in at least 1 cluster
de_known <- list()
for (ct in names(cell_types)) {
  ct_genes <- cell_types[[ct]]
  de_ct <- de_up[de_up$gene %in% ct_genes, ]
  if (nrow(de_ct) > 0) {
    # Pick top 5 by max logFC across clusters
    gene_maxfc <- tapply(de_ct$avg_log2FC, de_ct$gene, max)
    top_genes <- names(sort(gene_maxfc, decreasing = TRUE))[seq_len(min(5, length(gene_maxfc)))]
    de_known[[ct]] <- top_genes
  }
}

# Build gene list with labels
all_dot_genes <- c()
dot_labels <- c()
for (ct in names(de_known)) {
  for (g in de_known[[ct]]) {
    if (!(g %in% all_dot_genes)) {
      all_dot_genes <- c(all_dot_genes, g)
      lab <- get_gene_label(g, equiv)
      dot_labels <- c(dot_labels, setNames(paste0(lab, " [", ct, "]"), g))
    }
  }
}

if (length(all_dot_genes) > 0) {
  p_dot <- DotPlot(obj, features = all_dot_genes) +
    scale_x_discrete(labels = function(x) {
      sapply(x, function(g) {
        if (g %in% names(dot_labels)) dot_labels[g] else g
      })
    }) +
    labs(title = "Known cell type markers found among DE genes") +
    theme_pub() +
    theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5))
  ggsave(file.path(outdir, "04_dotplot_known_markers_in_DE.png"), p_dot,
         width = max(12, length(all_dot_genes) * 0.5), height = 7, dpi = dpi)
}


##=============================================================================
## 10. FEATURE PLOTS — TOP KNOWN MARKER PER CELL TYPE (FROM DE)
##=============================================================================

cat("\n=== 10. Feature plots: top known marker per cell type ===\n")

feat_plots <- list()
for (ct in names(de_known)) {
  g <- de_known[[ct]][1]
  lab <- get_gene_label(g, equiv)
  feat_plots[[ct]] <- FeaturePlot(obj, features = g, pt.size = 0.3) +
    scale_color_viridis_c() +
    labs(title = paste0(ct, ": ", lab)) +
    theme_pub()
}
if (length(feat_plots) > 0) {
  p_feat <- wrap_plots(feat_plots, ncol = 3)
  ggsave(file.path(outdir, "05_feature_top_known_markers.png"), p_feat,
         width = 16, height = ceiling(length(feat_plots) / 3) * 5, dpi = dpi)
}


##=============================================================================
## 11. UMAP WITH SUGGESTED ANNOTATION
##=============================================================================

cat("\n=== 11. Annotated UMAP ===\n")

# Create cluster → suggested label map
label_map <- setNames(combined_annot$Suggested_type, combined_annot$Cluster)
annot_vec <- label_map[as.character(Idents(obj))]
names(annot_vec) <- colnames(obj)
obj$annotation <- annot_vec

# UMAP with cluster numbers
p_umap_cl <- DimPlot(obj, group.by = "seurat_clusters", label = TRUE, pt.size = 0.3) +
  labs(title = "Clusters") +
  theme_pub() + NoLegend()

# UMAP with suggested annotation
p_umap_an <- DimPlot(obj, group.by = "annotation", label = TRUE,
                     label.size = 3, repel = TRUE, pt.size = 0.3) +
  labs(title = "Suggested annotation") +
  theme_pub() + NoLegend()

p_umap_both <- p_umap_cl | p_umap_an
ggsave(file.path(outdir, "06_UMAP_annotation.png"), p_umap_both,
       width = 18, height = 7, dpi = dpi)


##=============================================================================
## 12. CELL CYCLE PLOTS
##=============================================================================

cat("\n=== 12. Cell cycle plots ===\n")

if ("Phase" %in% colnames(obj@meta.data) && !all(obj$Phase == "Unknown")) {
  # UMAP colored by cell cycle phase
  p_cc_umap <- DimPlot(obj, group.by = "Phase", pt.size = 0.3,
                        cols = c("G1" = "grey80", "S" = "#E69F00", "G2M" = "#D55E00")) +
    labs(title = "Cell cycle phase") +
    theme_pub()

  # UMAP colored by S and G2M scores
  p_s <- FeaturePlot(obj, features = "S.Score", pt.size = 0.3) +
    scale_color_viridis_c(option = "magma") +
    labs(title = "S-phase score") + theme_pub()
  p_g2m <- FeaturePlot(obj, features = "G2M.Score", pt.size = 0.3) +
    scale_color_viridis_c(option = "magma") +
    labs(title = "G2M score") + theme_pub()

  p_cc_all <- p_cc_umap | p_s | p_g2m
  ggsave(file.path(outdir, "07_UMAP_cell_cycle.png"), p_cc_all,
         width = 21, height = 6, dpi = dpi)

  # Bar plot: cell cycle phase per cluster
  cc_data <- data.frame(
    Cluster = as.character(Idents(obj)),
    Phase   = obj$Phase
  )
  cc_data$Cluster <- factor(cc_data$Cluster,
                            levels = as.character(sort(as.integer(unique(cc_data$Cluster)))))
  cc_data$Phase <- factor(cc_data$Phase, levels = c("G1", "S", "G2M"))

  p_cc_bar <- ggplot(cc_data, aes(x = Cluster, fill = Phase)) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = c("G1" = "grey80", "S" = "#E69F00", "G2M" = "#D55E00")) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Cell cycle phase composition per cluster",
         x = "Cluster", y = "Proportion") +
    theme_pub()
  ggsave(file.path(outdir, "08_barplot_cell_cycle_per_cluster.png"), p_cc_bar,
         width = 10, height = 5, dpi = dpi)

  # Overall phase distribution (pie chart)
  phase_counts <- as.data.frame(table(obj$Phase))
  colnames(phase_counts) <- c("Phase", "Count")
  phase_counts$Pct <- round(100 * phase_counts$Count / sum(phase_counts$Count), 1)
  phase_counts$Label <- paste0(phase_counts$Phase, "\n", phase_counts$Count,
                               " (", phase_counts$Pct, "%)")
  phase_counts$Phase <- factor(phase_counts$Phase, levels = c("G1", "S", "G2M"))

  p_pie <- ggplot(phase_counts, aes(x = "", y = Count, fill = Phase)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_manual(values = c("G1" = "grey80", "S" = "#E69F00", "G2M" = "#D55E00")) +
    labs(title = "Overall cell cycle phase distribution") +
    theme_void() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(outdir, "09_pie_cell_cycle_overall.png"), p_pie,
         width = 7, height = 6, dpi = dpi)

  # S.Score vs G2M.Score scatter, colored by phase
  p_scatter <- ggplot(data.frame(S = obj$S.Score, G2M = obj$G2M.Score, Phase = obj$Phase),
                      aes(x = S, y = G2M, color = Phase)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_manual(values = c("G1" = "grey70", "S" = "#E69F00", "G2M" = "#D55E00")) +
    labs(title = "S-phase vs G2M score per cell",
         x = "S-phase score", y = "G2M score") +
    theme_pub()
  ggsave(file.path(outdir, "10_scatter_S_vs_G2M_score.png"), p_scatter,
         width = 8, height = 6, dpi = dpi)

  # Violin of S and G2M scores per cluster
  p_vln_s <- VlnPlot(obj, features = "S.Score", pt.size = 0) +
    labs(title = "S-phase score per cluster") +
    theme_pub() + NoLegend()
  p_vln_g2m <- VlnPlot(obj, features = "G2M.Score", pt.size = 0) +
    labs(title = "G2M score per cluster") +
    theme_pub() + NoLegend()
  p_vln_cc <- p_vln_s / p_vln_g2m
  ggsave(file.path(outdir, "11_violin_cell_cycle_scores.png"), p_vln_cc,
         width = 12, height = 8, dpi = dpi)

  cat("Cell cycle plots saved\n")
}


##=============================================================================
## 13. SAVE
##=============================================================================

cat("\n=== 13. Saving ===\n")

# CSV — combined annotation
write.csv(combined_annot, file.path(outdir, "cluster_annotation.csv"), row.names = FALSE)

# CSV — marker gene lists per cell type
marker_list_df <- do.call(rbind, lapply(names(cell_types), function(ct) {
  genes <- cell_types[[ct]]
  if (length(genes) == 0) return(NULL)
  data.frame(
    cell_type   = ct,
    gene        = genes,
    gene_symbol = sapply(genes, function(g) get_gene_label(g, equiv)),
    stringsAsFactors = FALSE
  )
}))
write.csv(marker_list_df, file.path(outdir, "celltype_marker_genes.csv"), row.names = FALSE)

# CSV — full DE overlap detail (which genes match for each cluster x cell type)
overlap_detail <- do.call(rbind, lapply(overlap_list, function(r) {
  do.call(rbind, lapply(names(cell_types), function(ct) {
    genes_str <- r[[paste0(ct, "_genes")]]
    if (nchar(genes_str) == 0) return(NULL)
    gene_ids <- strsplit(genes_str, ";")[[1]]
    data.frame(
      cluster     = r$cluster,
      cell_type   = ct,
      gene        = gene_ids,
      gene_symbol = sapply(gene_ids, function(g) get_gene_label(g, equiv)),
      stringsAsFactors = FALSE
    )
  }))
}))
if (!is.null(overlap_detail) && nrow(overlap_detail) > 0) {
  write.csv(overlap_detail, file.path(outdir, "DE_overlap_detail.csv"), row.names = FALSE)
}

# XLSX
if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  wb <- createWorkbook()
  hs <- createStyle(textDecoration = "bold", border = "bottom")

  # Sheet 1: Combined annotation
  addWorksheet(wb, "Annotation")
  writeData(wb, "Annotation", combined_annot)
  addStyle(wb, "Annotation", hs, rows = 1, cols = seq_len(ncol(combined_annot)))
  setColWidths(wb, "Annotation", cols = seq_len(ncol(combined_annot)), widths = "auto")

  # Sheet 2: DE overlap counts
  addWorksheet(wb, "DE overlap counts")
  writeData(wb, "DE overlap counts", overlap_mat)
  addStyle(wb, "DE overlap counts", hs, rows = 1, cols = seq_len(ncol(overlap_mat)))
  setColWidths(wb, "DE overlap counts", cols = seq_len(ncol(overlap_mat)), widths = "auto")

  # Sheet 3: DE overlap detail (gene-level)
  if (!is.null(overlap_detail) && nrow(overlap_detail) > 0) {
    addWorksheet(wb, "DE overlap detail")
    writeData(wb, "DE overlap detail", overlap_detail)
    addStyle(wb, "DE overlap detail", hs, rows = 1, cols = seq_len(ncol(overlap_detail)))
    setColWidths(wb, "DE overlap detail", cols = seq_len(ncol(overlap_detail)), widths = "auto")
  }

  # Sheet 4: Module scores per cluster
  score_summary <- do.call(rbind, lapply(levels(Idents(obj)), function(cl) {
    cells <- WhichCells(obj, idents = cl)
    m <- obj@meta.data[cells, ]
    row <- data.frame(Cluster = cl)
    for (sc in score_cols) {
      ct_name <- gsub("score_", "", sc)
      row[[ct_name]] <- round(mean(m[[sc]], na.rm = TRUE), 4)
    }
    row
  }))
  addWorksheet(wb, "Module scores")
  writeData(wb, "Module scores", score_summary)
  addStyle(wb, "Module scores", hs, rows = 1, cols = seq_len(ncol(score_summary)))
  setColWidths(wb, "Module scores", cols = seq_len(ncol(score_summary)), widths = "auto")

  # Sheet 5: Cell cycle per cluster
  if (!is.null(cc_per_cluster)) {
    addWorksheet(wb, "Cell cycle")
    writeData(wb, "Cell cycle", cc_per_cluster)
    addStyle(wb, "Cell cycle", hs, rows = 1, cols = seq_len(ncol(cc_per_cluster)))
    setColWidths(wb, "Cell cycle", cols = seq_len(ncol(cc_per_cluster)), widths = "auto")
  }

  # Sheet 6: Known marker genes used
  addWorksheet(wb, "Marker genes")
  writeData(wb, "Marker genes", marker_list_df)
  addStyle(wb, "Marker genes", hs, rows = 1, cols = seq_len(ncol(marker_list_df)))
  setColWidths(wb, "Marker genes", cols = seq_len(ncol(marker_list_df)), widths = "auto")

  saveWorkbook(wb, file.path(outdir, "annotation_summary.xlsx"), overwrite = TRUE)
  cat("XLSX saved: annotation_summary.xlsx\n")
}

cat("\n=== ANNOTATION COMPLETED ===\n")
cat("Results in:", outdir, "\n")
