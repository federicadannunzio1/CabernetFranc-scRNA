##=============================================================================
## PHASE 7: Cell type annotation using known leaf markers
##
## Input:  03_scRNA/04_analysis/VvProt_clustered.rds
## Output: 03_scRNA/04_analysis/annotation/ (plots + tables)
##
## Strategy: use known marker genes from Arabidopsis leaf scRNA-seq atlases
## (Kim 2021, Lopez-Anido 2021, etc.) and find grapevine orthologs via
## Grapedia gene symbols and functional annotation keywords.
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
    theme(axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.4))
}

# --- Paths ----
proj   <- "/lustre/home/dcapozzi/projects/scRNA_protoplasts"
indir  <- file.path(proj, "03_scRNA/04_analysis")
outdir <- file.path(indir, "annotation")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dpi <- 300

# Reference files
equiv_file <- file.path(proj, "01_reference/all_on_T2T_ref_equivalences_filtered_no_copies_symbols.csv")
annot_file <- file.path(proj, "01_reference/5.1_on_T2T_ref_functional_annotation_summary_with_coordinates.tsv")


##=============================================================================
## 1. LOAD DATA
##=============================================================================

cat("=== 1. Loading data ===\n")
obj <- readRDS(file.path(indir, "VvProt_clustered.rds"))
DefaultAssay(obj) <- "SCT"

equiv <- read.delim(equiv_file, header = TRUE, stringsAsFactors = FALSE)
annot <- read.delim(annot_file, header = TRUE, stringsAsFactors = FALSE)
annot$gene_seurat <- gsub("_", "-", annot$gene)

cat(sprintf("Cells: %d, Clusters: %d\n", ncol(obj), length(levels(Idents(obj)))))


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
find_by_annot <- function(pattern, annot, obj, cols = c("Pfam", "PANTHER", "GO_biological_process", "GO_molecular_function", "eggNOG")) {
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
  CAB = find_by_symbol("^CAB[0-9]", equiv, obj),
  PSBS = find_by_symbol("^PSBS", equiv, obj),
  RCA = find_by_symbol("^RCA$", equiv, obj),
  photosystem = find_by_annot("photosystem.II|oxygen.evolving|PSII", annot, obj)
)

# EPIDERMIS: cuticle, wax, lipid transfer
epidermis_markers <- list(
  PDF = find_by_symbol("^PDF[12]", equiv, obj),
  CER = find_by_symbol("^CER[0-9]", equiv, obj),
  KCS = find_by_symbol("^KCS[0-9]", equiv, obj),
  LTP = find_by_annot("lipid.transfer.protein|nsLTP", annot, obj),
  GDSL = find_by_annot("GDSL.*lipase|GDSL.*esterase", annot, obj),
  ATML1 = find_by_symbol("^ATML1", equiv, obj),
  FDH = find_by_symbol("^FDH$", equiv, obj)
)

# GUARD CELLS: stomata
guard_markers <- list(
  SLAC1 = find_by_symbol("^SLAC[0-9]", equiv, obj),
  KAT = find_by_symbol("^KAT[0-9]", equiv, obj),
  MYB60 = find_by_symbol("^MYB60", equiv, obj),
  FAMA = find_by_symbol("^FAMA", equiv, obj),
  OST1 = find_by_symbol("^OST1", equiv, obj),
  slow_anion = find_by_annot("slow.anion.channel|SLAC", annot, obj),
  stomatal = find_by_annot("stomat", annot, obj)
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
  SEOR = find_by_symbol("^SEOR", equiv, obj),
  PP2 = find_by_symbol("^PP2$|^PP2-", equiv, obj),
  CLE41 = find_by_symbol("^CLE41|^CLE44", equiv, obj),
  APL = find_by_symbol("^APL$", equiv, obj),
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
  PCNA = find_by_symbol("^PCNA", equiv, obj),
  cyclin = find_by_annot("^cyclin[- ]|cyclin.dependent", annot, obj),
  MCM = find_by_symbol("^MCM[0-9]", equiv, obj)
)

# Collect all cell type marker lists
cell_types <- list(
  Mesophyll = unique(unlist(mesophyll_markers)),
  Epidermis = unique(unlist(epidermis_markers)),
  Guard_cells = unique(unlist(guard_markers)),
  Xylem = unique(unlist(xylem_markers)),
  Phloem = unique(unlist(phloem_markers)),
  Companion_cells = unique(unlist(companion_markers)),
  Bundle_sheath = unique(unlist(bundle_sheath_markers)),
  Dividing = unique(unlist(dividing_markers))
)

# Report
cat("\nMarker genes found per cell type:\n")
for (ct in names(cell_types)) {
  cat(sprintf("  %-20s: %d genes\n", ct, length(cell_types[[ct]])))
}


##=============================================================================
## 3. MODULE SCORE PER CELL TYPE
##=============================================================================

cat("\n=== 3. Computing module scores ===\n")

for (ct in names(cell_types)) {
  genes <- cell_types[[ct]]
  if (length(genes) >= 3) {
    obj <- AddModuleScore(obj, features = list(genes),
                          name = paste0("score_", ct), seed = 42)
    # AddModuleScore appends "1" to the name
    colnames(obj@meta.data)[ncol(obj@meta.data)] <- paste0("score_", ct)
    cat(sprintf("  %s: module score computed (%d genes)\n", ct, length(genes)))
  } else {
    cat(sprintf("  %s: skipped (only %d genes found)\n", ct, length(genes)))
  }
}


##=============================================================================
## 4. PLOTS — MODULE SCORES ON UMAP
##=============================================================================

cat("\n=== 4. Plotting module scores ===\n")

score_cols <- grep("^score_", colnames(obj@meta.data), value = TRUE)

# 4.1 UMAP colored by module scores
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

# 4.2 Violin of module scores per cluster
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
## 5. DOT PLOT — KEY MARKERS PER CELL TYPE
##=============================================================================

cat("\n=== 5. Key marker dot plots ===\n")

# Helper: get gene label
get_gene_label <- function(gene_seurat, equiv) {
  gid <- gsub("-", "_", gene_seurat)
  idx <- match(gid, equiv$target_id)
  if (!is.na(idx) && !is.na(equiv$T2T.v5.1_gene_symbol[idx]) &&
      nchar(equiv$T2T.v5.1_gene_symbol[idx]) > 0) {
    return(equiv$T2T.v5.1_gene_symbol[idx])
  }
  return(gene_seurat)
}

# Select top representative markers per cell type (max 5 each)
key_markers <- list()
for (ct in names(cell_types)) {
  genes <- cell_types[[ct]]
  if (length(genes) > 5) {
    # Pick the most variable ones
    var_data <- GetAssayData(obj, layer = "data")[genes, ]
    gene_var <- apply(var_data, 1, var)
    genes <- names(sort(gene_var, decreasing = TRUE))[1:5]
  }
  key_markers[[ct]] <- genes
}

# Build ordered gene list with labels
all_key_genes <- c()
gene_labels_vec <- c()
for (ct in names(key_markers)) {
  for (g in key_markers[[ct]]) {
    if (!(g %in% all_key_genes)) {
      all_key_genes <- c(all_key_genes, g)
      lab <- get_gene_label(g, equiv)
      gene_labels_vec <- c(gene_labels_vec,
                           setNames(paste0(lab, " [", ct, "]"), g))
    }
  }
}

if (length(all_key_genes) > 0) {
  p_dot_ct <- DotPlot(obj, features = all_key_genes) +
    scale_x_discrete(labels = function(x) {
      sapply(x, function(g) {
        if (g %in% names(gene_labels_vec)) gene_labels_vec[g] else g
      })
    }) +
    labs(title = "Key markers by cell type") +
    theme_pub() +
    theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5))
  ggsave(file.path(outdir, "03_dotplot_celltype_markers.png"), p_dot_ct,
         width = 22, height = 7, dpi = dpi)
}


##=============================================================================
## 6. FEATURE PLOTS — TOP MARKER PER CELL TYPE
##=============================================================================

cat("\n=== 6. Feature plots per cell type ===\n")

feat_plots <- list()
for (ct in names(key_markers)) {
  genes <- key_markers[[ct]]
  if (length(genes) >= 1) {
    g <- genes[1]  # top marker
    lab <- get_gene_label(g, equiv)
    feat_plots[[ct]] <- FeaturePlot(obj, features = g, pt.size = 0.3) +
      scale_color_viridis_c() +
      labs(title = paste0(ct, ": ", lab)) +
      theme_pub()
  }
}
if (length(feat_plots) > 0) {
  p_feat_ct <- wrap_plots(feat_plots, ncol = 3)
  ggsave(file.path(outdir, "04_feature_celltype_top_markers.png"), p_feat_ct,
         width = 16, height = ceiling(length(feat_plots) / 3) * 5, dpi = dpi)
}


##=============================================================================
## 7. ANNOTATION SUMMARY TABLE
##=============================================================================

cat("\n=== 7. Building annotation summary ===\n")

# For each cluster, find which cell type has the highest module score
score_cols_available <- grep("^score_", colnames(obj@meta.data), value = TRUE)

annotation_table <- do.call(rbind, lapply(levels(Idents(obj)), function(cl) {
  cells <- WhichCells(obj, idents = cl)
  m <- obj@meta.data[cells, ]

  # Mean score per cell type
  scores <- sapply(score_cols_available, function(sc) mean(m[[sc]], na.rm = TRUE))
  names(scores) <- gsub("score_", "", names(scores))

  best_ct <- names(which.max(scores))
  best_score <- max(scores)
  second_ct <- names(sort(scores, decreasing = TRUE))[2]
  second_score <- sort(scores, decreasing = TRUE)[2]

  data.frame(
    Cluster = cl,
    N_cells = length(cells),
    Pct_cells = round(100 * length(cells) / ncol(obj), 1),
    Best_celltype = best_ct,
    Best_score = round(best_score, 4),
    Second_celltype = second_ct,
    Second_score = round(second_score, 4),
    Median_stress = round(median(m$percent.stress), 1),
    Median_chloro = round(median(m$percent.chloro), 1)
  )
}))

# Print
cat("\n=== ANNOTATION SUMMARY ===\n\n")
cat(sprintf("%-8s %6s %6s  %-18s %8s  %-18s %8s  %7s %7s\n",
            "Cluster", "Cells", "Pct", "Best type", "Score",
            "2nd type", "Score", "Stress", "Chloro"))
cat(paste(rep("-", 100), collapse = ""), "\n")
for (i in seq_len(nrow(annotation_table))) {
  r <- annotation_table[i, ]
  cat(sprintf("%-8s %6d %5.1f%%  %-18s %8.4f  %-18s %8.4f  %6.1f%% %6.1f%%\n",
              r$Cluster, r$N_cells, r$Pct_cells,
              r$Best_celltype, r$Best_score,
              r$Second_celltype, r$Second_score,
              r$Median_stress, r$Median_chloro))
}


##=============================================================================
## 8. SAVE
##=============================================================================

cat("\n=== 8. Saving ===\n")

# CSV
write.csv(annotation_table, file.path(outdir, "cluster_annotation.csv"), row.names = FALSE)

# Marker gene lists per cell type
marker_list_df <- do.call(rbind, lapply(names(cell_types), function(ct) {
  genes <- cell_types[[ct]]
  if (length(genes) == 0) return(NULL)
  data.frame(
    cell_type = ct,
    gene = genes,
    gene_symbol = sapply(genes, function(g) get_gene_label(g, equiv)),
    stringsAsFactors = FALSE
  )
}))
write.csv(marker_list_df, file.path(outdir, "celltype_marker_genes.csv"), row.names = FALSE)

# XLSX
if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  wb <- createWorkbook()
  headerStyle <- createStyle(textDecoration = "bold", border = "bottom")

  # Sheet 1: Annotation
  addWorksheet(wb, "Cluster annotation")
  writeData(wb, "Cluster annotation", annotation_table)
  addStyle(wb, "Cluster annotation", headerStyle, rows = 1, cols = 1:ncol(annotation_table))
  setColWidths(wb, "Cluster annotation", cols = 1:ncol(annotation_table), widths = "auto")

  # Sheet 2: Marker genes used
  addWorksheet(wb, "Marker genes")
  writeData(wb, "Marker genes", marker_list_df)
  addStyle(wb, "Marker genes", headerStyle, rows = 1, cols = 1:ncol(marker_list_df))
  setColWidths(wb, "Marker genes", cols = 1:ncol(marker_list_df), widths = "auto")

  # Sheet 3: Module scores summary per cluster
  score_summary <- do.call(rbind, lapply(levels(Idents(obj)), function(cl) {
    cells <- WhichCells(obj, idents = cl)
    m <- obj@meta.data[cells, ]
    row <- data.frame(Cluster = cl)
    for (sc in score_cols_available) {
      ct_name <- gsub("score_", "", sc)
      row[[ct_name]] <- round(mean(m[[sc]], na.rm = TRUE), 4)
    }
    row
  }))
  addWorksheet(wb, "Module scores")
  writeData(wb, "Module scores", score_summary)
  addStyle(wb, "Module scores", headerStyle, rows = 1, cols = 1:ncol(score_summary))
  setColWidths(wb, "Module scores", cols = 1:ncol(score_summary), widths = "auto")

  saveWorkbook(wb, file.path(outdir, "annotation_summary.xlsx"), overwrite = TRUE)
  cat("XLSX saved: annotation_summary.xlsx\n")
}

cat("\n=== ANNOTATION COMPLETED ===\n")
cat("Results in:", outdir, "\n")
