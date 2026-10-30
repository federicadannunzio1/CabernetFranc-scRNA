##=============================================================================
## FASE 3+4: QC completo — scRNA-seq protoplasti Vitis vinifera (Cabernet Franc)
##
## Ordine operazioni:
##   1. SoupX — correzione ambient RNA (usa matrice raw + filtrata)
##   2. Metriche QC per cellula (nCount, nFeature, complessita', ribo, stress)
##   3. scDblFinder — predizione doublet
##   4. Plot diagnostici con tutte le metriche + doublet score
##   5. Soglie suggerite (MAD-based) — NON filtra, decidiamo dopo
##
## Eseguire su TeraStat2: sbatch 03_qc.slurm
##=============================================================================

library(Seurat)
library(SoupX)
library(scDblFinder)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)

set.seed(42)

# --- Paths ----
proj   <- "/lustre/home/dcapozzi/projects/scRNA_protoplasts"
outdir <- file.path(proj, "03_scRNA/02_qc")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cr_outs      <- file.path(proj, "03_scRNA/01_alignment/VvProt_scRNA/outs")
filtered_dir <- file.path(cr_outs, "filtered_feature_bc_matrix")
raw_dir      <- file.path(cr_outs, "raw_feature_bc_matrix")


##=============================================================================
## 1. SoupX — RIMOZIONE AMBIENT RNA
##=============================================================================
## SoupX stima la "zuppa" di RNA ambientale dai barcode vuoti (matrice raw)
## e la sottrae dalla matrice filtrata (cellule vere).
## Usa il clustering di Cell Ranger per stimare i profili di contaminazione.

cat("=== 1. SoupX: correzione ambient RNA ===\n")

# Carica con il metodo automatico di SoupX (legge clustering CR da analysis/)
sc <- load10X(cr_outs)

# Stima automatica della frazione di contaminazione
sc <- autoEstCont(sc)
cat(sprintf("Frazione contaminazione stimata: %.3f (%.1f%%)\n",
            sc$fit$rhoEst, 100 * sc$fit$rhoEst))

# Correggi le conte
counts_corrected <- adjustCounts(sc)

cat(sprintf("Dimensioni matrice corretta: %d geni x %d cellule\n",
            nrow(counts_corrected), ncol(counts_corrected)))

# Plot contaminazione SoupX
tryCatch({
  png(file.path(outdir, "00_soupx_contamination.png"), width = 800, height = 600, res = 150)
  plot(sc, "rhoEst")
  dev.off()
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  cat("NOTA: plot SoupX fallito (non critico):", conditionMessage(e), "\n")
})


##=============================================================================
## 2. CREA OGGETTO SEURAT CON CONTE CORRETTE
##=============================================================================

cat("\n=== 2. Creazione oggetto Seurat (conte SoupX-corrette) ===\n")

obj <- CreateSeuratObject(counts = counts_corrected, project = "VvProt",
                          min.cells = 3, min.features = 200)

cat(sprintf("Cellule: %d (dopo min.cells=3, min.features=200)\n", ncol(obj)))
cat(sprintf("Geni: %d\n", nrow(obj)))
cat(sprintf("Esempio nomi geni: %s\n", paste(head(rownames(obj), 10), collapse = ", ")))


##=============================================================================
## 3. METRICHE QC PER CELLULA
##=============================================================================

cat("\n=== 3. Calcolo metriche QC ===\n")

# 3.1 Geni ribosomali
#     I rownames possono essere ID (Vitvi05_...) o nomi genici.
ribo_by_name <- grep("(?i)(^RPL|^RPS[0-9]|^MRPL|^MRPS|ribosomal protein)",
                     rownames(obj), value = TRUE, perl = TRUE)

if (length(ribo_by_name) > 0) {
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, features = ribo_by_name)
  cat(sprintf("Geni ribosomali trovati (per nome): %d\n", length(ribo_by_name)))
} else {
  cat("NOTA: nessun gene ribosomale trovato per pattern nome.\n")
  cat("I rownames sono probabilmente ID Vitvi05_. Cercando con pattern ID...\n")
  # Se i nomi non funzionano, percent.ribo resta NA
  obj[["percent.ribo"]] <- NA
}

# 3.2 Geni mitocondriali / plastidiali
#     chr00 escluso dal reference -> le sequenze organellari NON mappano.
#     NON possiamo calcolare %mito/%chloro classico.
#     Cerchiamo geni nucleari che codificano proteine organellari.
#     ATTENZIONE: questi NON sono contaminazione, sono geni nucleari reali.

mito_genes <- grep("(?i)(^NAD[0-9]|^COX[0-9]|^ATP[0-9]|^CYTB|^COB|^SDH[0-9])",
                   rownames(obj), value = TRUE, perl = TRUE)
chloro_genes <- grep("(?i)(^PSB[A-Z]|^PSA[A-Z]|^LHCB|^LHCA|^RBCS|^RBCL|^PETC|^ATPC)",
                     rownames(obj), value = TRUE, perl = TRUE)

cat(sprintf("Geni con pattern mitocondriale nucleare: %d\n", length(mito_genes)))
if (length(mito_genes) > 0) cat("  ", paste(head(mito_genes, 15), collapse = ", "), "\n")
cat(sprintf("Geni con pattern plastidiale nucleare: %d\n", length(chloro_genes)))
if (length(chloro_genes) > 0) cat("  ", paste(head(chloro_genes, 15), collapse = ", "), "\n")

if (length(chloro_genes) > 0) {
  obj[["percent.chloro"]] <- PercentageFeatureSet(obj, features = chloro_genes)
}

# 3.3 Complessita': log10(nFeature) / log10(nCount)
#     Valori bassi = pochi geni su molti UMI -> dying cells, ambient RNA,
#     o cellule dominate da pochi trascritti altamente espressi
obj[["log10GenesPerUMI"]] <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)

# 3.4 Geni stress / protoplastizzazione
#     La digestione enzimatica della parete induce geni di stress.
stress_genes <- grep("(?i)(^HSP|^WRKY|^ERF|^JAZ|^MYC2|^PR-|^PR[0-9]|^LOX|^AOS|^AOC|^OPR|^PAL|^CHS|wound|defense|pathogen|chitinase|glucanase|peroxidase)",
                     rownames(obj), value = TRUE, perl = TRUE)

cat(sprintf("Geni stress/difesa trovati per pattern: %d\n", length(stress_genes)))
if (length(stress_genes) > 0) {
  cat("  ", paste(head(stress_genes, 20), collapse = ", "), "\n")
  obj[["percent.stress"]] <- PercentageFeatureSet(obj, features = stress_genes)
}


##=============================================================================
## 4. scDblFinder — PREDIZIONE DOUBLET
##=============================================================================

cat("\n=== 4. scDblFinder: predizione doublet ===\n")

# Converti in SingleCellExperiment per scDblFinder
sce <- as.SingleCellExperiment(obj)
sce <- scDblFinder(sce)

# Trasferisci risultati nel Seurat object
obj[["scDblFinder.score"]]  <- sce$scDblFinder.score
obj[["scDblFinder.class"]]  <- sce$scDblFinder.class

n_singlet <- sum(obj$scDblFinder.class == "singlet")
n_doublet <- sum(obj$scDblFinder.class == "doublet")
cat(sprintf("Singlet: %d (%.1f%%)\n", n_singlet, 100 * n_singlet / ncol(obj)))
cat(sprintf("Doublet: %d (%.1f%%)\n", n_doublet, 100 * n_doublet / ncol(obj)))

rm(sce)


##=============================================================================
## 5. STATISTICHE RIASSUNTIVE
##=============================================================================

cat("\n=== STATISTICHE QC PRE-FILTRAGGIO ===\n")
cat(sprintf("Cellule totali: %d\n", ncol(obj)))

qc_summary <- function(x, name) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(invisible(NULL))
  cat(sprintf("  %-20s min=%-8.1f Q1=%-8.1f median=%-8.1f Q3=%-8.1f max=%-8.1f\n",
              name, min(x), quantile(x, 0.25), median(x), quantile(x, 0.75), max(x)))
}

qc_summary(obj$nCount_RNA, "nCount_RNA")
qc_summary(obj$nFeature_RNA, "nFeature_RNA")
qc_summary(obj$log10GenesPerUMI, "log10GenesPerUMI")
qc_summary(obj$scDblFinder.score, "doublet.score")
if (!all(is.na(obj$percent.ribo)))      qc_summary(obj$percent.ribo, "percent.ribo")
if (!is.null(obj[["percent.chloro"]]))  qc_summary(obj$percent.chloro, "percent.chloro")
if (!is.null(obj[["percent.stress"]]))  qc_summary(obj$percent.stress, "percent.stress")


##=============================================================================
## 6. SOGLIE SUGGERITE (MAD)
##=============================================================================

log_nCount   <- log10(obj$nCount_RNA)
log_nFeature <- log10(obj$nFeature_RNA)

nCount_med   <- median(log_nCount)
nCount_mad   <- mad(log_nCount)
nFeature_med <- median(log_nFeature)
nFeature_mad <- mad(log_nFeature)

nCount_lo   <- 10^(nCount_med - 3 * nCount_mad)
nCount_hi   <- 10^(nCount_med + 3 * nCount_mad)
nFeature_lo <- 10^(nFeature_med - 3 * nFeature_mad)
nFeature_hi <- 10^(nFeature_med + 3 * nFeature_mad)

cat("\n=== SOGLIE SUGGERITE (3 MAD su log10) ===\n")
cat(sprintf("  nCount_RNA:       [%.0f — %.0f]\n", nCount_lo, nCount_hi))
cat(sprintf("  nFeature_RNA:     [%.0f — %.0f]\n", nFeature_lo, nFeature_hi))
cat(sprintf("  log10GenesPerUMI: > 0.80\n"))
cat(sprintf("  scDblFinder:      class == 'singlet'\n"))

# Conta cellule che passano ogni filtro
pass_count   <- obj$nCount_RNA >= nCount_lo & obj$nCount_RNA <= nCount_hi
pass_feature <- obj$nFeature_RNA >= nFeature_lo & obj$nFeature_RNA <= nFeature_hi
pass_complex <- obj$log10GenesPerUMI > 0.80
pass_singlet <- obj$scDblFinder.class == "singlet"
pass_all     <- pass_count & pass_feature & pass_complex & pass_singlet

cat(sprintf("\n  Passano nCount:       %d / %d (%.1f%%)\n",
            sum(pass_count), ncol(obj), 100 * mean(pass_count)))
cat(sprintf("  Passano nFeature:     %d / %d (%.1f%%)\n",
            sum(pass_feature), ncol(obj), 100 * mean(pass_feature)))
cat(sprintf("  Passano complessita': %d / %d (%.1f%%)\n",
            sum(pass_complex), ncol(obj), 100 * mean(pass_complex)))
cat(sprintf("  Passano singlet:      %d / %d (%.1f%%)\n",
            sum(pass_singlet), ncol(obj), 100 * mean(pass_singlet)))
cat(sprintf("  Passano TUTTI:        %d / %d (%.1f%%)\n",
            sum(pass_all), ncol(obj), 100 * mean(pass_all)))
cat("==========================================\n")


##=============================================================================
## 7. PLOT QC
##=============================================================================

cat("\n=== 7. Generazione plot ===\n")

# 7.1 Violin nFeature e nCount
p_vln <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA"),
                 ncol = 2, pt.size = 0) & theme_minimal() & NoLegend()
ggsave(file.path(outdir, "01_violin_nFeature_nCount.png"), p_vln,
       width = 10, height = 5, dpi = 200)

# 7.2 Violin metriche extra (complessita', ribo, chloro, stress)
extra_features <- "log10GenesPerUMI"
if (!all(is.na(obj$percent.ribo)))      extra_features <- c(extra_features, "percent.ribo")
if (!is.null(obj[["percent.chloro"]]))  extra_features <- c(extra_features, "percent.chloro")
if (!is.null(obj[["percent.stress"]]))  extra_features <- c(extra_features, "percent.stress")

p_vln2 <- VlnPlot(obj, features = extra_features,
                  ncol = length(extra_features), pt.size = 0) &
  theme_minimal() & NoLegend()
ggsave(file.path(outdir, "02_violin_complexity_extra.png"), p_vln2,
       width = 5 * length(extra_features), height = 5, dpi = 200)

# 7.3 Scatter nCount vs nFeature
p_scatter1 <- ggplot(obj@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(size = 0.2, alpha = 0.2, color = "grey30") +
  scale_x_log10() + scale_y_log10() +
  geom_smooth(method = "loess", se = FALSE, color = "dodgerblue", linewidth = 0.7) +
  labs(title = "nCount vs nFeature (log-log)", x = "UMI totali", y = "Geni rilevati") +
  theme_minimal()
ggsave(file.path(outdir, "03_scatter_count_vs_feature.png"), p_scatter1,
       width = 7, height = 6, dpi = 200)

# 7.4 Densita' nCount e nFeature con soglie MAD
p_dens1 <- ggplot(obj@meta.data, aes(x = nCount_RNA)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  scale_x_log10() +
  geom_vline(xintercept = c(nCount_lo, nCount_hi), linetype = "dashed", color = "red") +
  annotate("text", x = nCount_lo, y = Inf, label = sprintf("%.0f", nCount_lo),
           vjust = 2, hjust = 1.1, color = "red", size = 3) +
  annotate("text", x = nCount_hi, y = Inf, label = sprintf("%.0f", nCount_hi),
           vjust = 2, hjust = -0.1, color = "red", size = 3) +
  labs(title = "Distribuzione UMI per cellula", x = "nCount_RNA (log10)", y = "Densita'") +
  theme_minimal()

p_dens2 <- ggplot(obj@meta.data, aes(x = nFeature_RNA)) +
  geom_density(fill = "coral", alpha = 0.5) +
  scale_x_log10() +
  geom_vline(xintercept = c(nFeature_lo, nFeature_hi), linetype = "dashed", color = "red") +
  annotate("text", x = nFeature_lo, y = Inf, label = sprintf("%.0f", nFeature_lo),
           vjust = 2, hjust = 1.1, color = "red", size = 3) +
  annotate("text", x = nFeature_hi, y = Inf, label = sprintf("%.0f", nFeature_hi),
           vjust = 2, hjust = -0.1, color = "red", size = 3) +
  labs(title = "Distribuzione geni per cellula", x = "nFeature_RNA (log10)", y = "Densita'") +
  theme_minimal()

p_dens <- p_dens1 + p_dens2
ggsave(file.path(outdir, "04_density_count_feature.png"), p_dens,
       width = 12, height = 5, dpi = 200)

# 7.5 Densita' complessita'
p_complex <- ggplot(obj@meta.data, aes(x = log10GenesPerUMI)) +
  geom_density(fill = "darkgreen", alpha = 0.5) +
  geom_vline(xintercept = 0.80, linetype = "dashed", color = "red") +
  labs(title = "Complessita' trascrizionale",
       x = "log10(Genes)/log10(UMI)", y = "Densita'") +
  theme_minimal()
ggsave(file.path(outdir, "05_density_complexity.png"), p_complex,
       width = 6, height = 5, dpi = 200)

# 7.6 Scatter con soglie MAD sovrapposte
p_thresh <- ggplot(obj@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(size = 0.2, alpha = 0.2, color = "grey30") +
  scale_x_log10() + scale_y_log10() +
  geom_vline(xintercept = c(nCount_lo, nCount_hi), linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(nFeature_lo, nFeature_hi), linetype = "dashed", color = "blue") +
  labs(title = "nCount vs nFeature con soglie 3-MAD",
       subtitle = sprintf("nCount [%.0f-%.0f], nFeature [%.0f-%.0f]",
                          nCount_lo, nCount_hi, nFeature_lo, nFeature_hi),
       x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
  theme_minimal()
ggsave(file.path(outdir, "06_scatter_thresholds.png"), p_thresh,
       width = 8, height = 7, dpi = 200)

# 7.7 Scatter nCount vs complessita'
p_sc_compl <- ggplot(obj@meta.data, aes(x = nCount_RNA, y = log10GenesPerUMI)) +
  geom_point(size = 0.2, alpha = 0.2, color = "grey30") +
  scale_x_log10() +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
  labs(title = "UMI vs Complessita'",
       x = "nCount_RNA (log10)", y = "log10GenesPerUMI") +
  theme_minimal()
ggsave(file.path(outdir, "07_scatter_count_vs_complexity.png"), p_sc_compl,
       width = 7, height = 6, dpi = 200)

# 7.8 Scatter stress vs nCount
if (!is.null(obj[["percent.stress"]])) {
  p_stress <- ggplot(obj@meta.data, aes(x = nCount_RNA, y = percent.stress)) +
    geom_point(size = 0.2, alpha = 0.2, color = "grey30") +
    scale_x_log10() +
    labs(title = "UMI vs % geni stress/protoplastizzazione",
         x = "nCount_RNA (log10)", y = "% UMI geni stress") +
    theme_minimal()
  ggsave(file.path(outdir, "08_scatter_stress.png"), p_stress,
         width = 7, height = 6, dpi = 200)
}

# 7.9 Doublet score — distribuzione e scatter
p_dbl_dens <- ggplot(obj@meta.data, aes(x = scDblFinder.score, fill = scDblFinder.class)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("singlet" = "steelblue", "doublet" = "red")) +
  labs(title = "Distribuzione doublet score (scDblFinder)",
       x = "Doublet score", y = "Densita'", fill = "Classe") +
  theme_minimal()
ggsave(file.path(outdir, "09_density_doublet_score.png"), p_dbl_dens,
       width = 7, height = 5, dpi = 200)

p_dbl_scatter <- ggplot(obj@meta.data,
                        aes(x = nCount_RNA, y = nFeature_RNA, color = scDblFinder.class)) +
  geom_point(size = 0.3, alpha = 0.3) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = c("singlet" = "grey60", "doublet" = "red")) +
  labs(title = "nCount vs nFeature — doublet evidenziati",
       x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)", color = "Classe") +
  theme_minimal()
ggsave(file.path(outdir, "10_scatter_doublets.png"), p_dbl_scatter,
       width = 8, height = 7, dpi = 200)

# Doublet per posizione nello spazio: i doublet tendono ad avere nCount/nFeature alti
p_dbl_violin <- VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA"),
                        group.by = "scDblFinder.class", ncol = 2, pt.size = 0) &
  theme_minimal()
ggsave(file.path(outdir, "11_violin_singlet_vs_doublet.png"), p_dbl_violin,
       width = 10, height = 5, dpi = 200)


##=============================================================================
## 8. TOP 30 GENI PIU' ESPRESSI
##=============================================================================

counts_mat <- GetAssayData(obj, layer = "counts")
gene_totals <- sort(Matrix::rowSums(counts_mat), decreasing = TRUE)
total_umi <- sum(gene_totals)
top30 <- head(gene_totals, 30)

cat("\n=== TOP 30 GENI PIU' ESPRESSI ===\n")
for (i in seq_along(top30)) {
  cat(sprintf("  %2d. %-35s %10d UMI  (%.2f%%)\n",
              i, names(top30)[i], top30[i], 100 * top30[i] / total_umi))
}

top_pct <- cumsum(top30) / total_umi * 100
cat(sprintf("\n  Top 10 geni: %.1f%% UMI totali\n", top_pct[10]))
cat(sprintf("  Top 30 geni: %.1f%% UMI totali\n", top_pct[30]))

# Plot
top30_df <- data.frame(
  gene = factor(names(top30), levels = rev(names(top30))),
  total_umi = as.numeric(top30),
  pct = 100 * as.numeric(top30) / total_umi
)
p_top <- ggplot(top30_df, aes(x = total_umi, y = gene)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = sprintf("%.1f%%", pct)), hjust = -0.1, size = 2.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Top 30 geni piu' espressi (UMI totali)", x = "UMI totali", y = "") +
  theme_minimal()
ggsave(file.path(outdir, "12_top30_genes.png"), p_top, width = 9, height = 8, dpi = 200)

# % cellule che esprimono i top 10
cat("\n=== % CELLULE CHE ESPRIMONO I TOP 10 GENI ===\n")
for (i in 1:min(10, length(top30))) {
  gene_name <- names(top30)[i]
  pct_cells <- 100 * sum(counts_mat[gene_name, ] > 0) / ncol(counts_mat)
  cat(sprintf("  %-35s espresso in %.1f%% delle cellule\n", gene_name, pct_cells))
}


##=============================================================================
## 9. ISTOGRAMMA GENI PER CELLULA
##=============================================================================

p_hist <- ggplot(obj@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100, fill = "grey40", color = "white") +
  geom_vline(xintercept = c(nFeature_lo, nFeature_hi), linetype = "dashed", color = "red") +
  labs(title = "Istogramma geni per cellula", x = "nFeature_RNA", y = "Numero cellule") +
  theme_minimal()
ggsave(file.path(outdir, "13_histogram_nFeature.png"), p_hist,
       width = 7, height = 5, dpi = 200)


##=============================================================================
## 10. SALVA OGGETTO E METADATA
##=============================================================================

write.csv(obj@meta.data, file.path(outdir, "metadata_pre_filter.csv"), row.names = TRUE)
saveRDS(obj, file.path(outdir, "VvProt_prefilter.rds"))

# Salva anche le soglie suggerite in un file per riferimento
soglie <- data.frame(
  metrica = c("nCount_RNA_min", "nCount_RNA_max",
              "nFeature_RNA_min", "nFeature_RNA_max",
              "log10GenesPerUMI_min", "doublet_filter"),
  valore = c(round(nCount_lo), round(nCount_hi),
             round(nFeature_lo), round(nFeature_hi),
             0.80, NA),
  nota = c("3 MAD (log10)", "3 MAD (log10)",
           "3 MAD (log10)", "3 MAD (log10)",
           "standard", "scDblFinder.class == singlet")
)
write.csv(soglie, file.path(outdir, "soglie_suggerite.csv"), row.names = FALSE)

cat("\n=== QC COMPLETATO ===\n")
cat("Plot salvati in:", outdir, "\n")
cat("Oggetto Seurat (SoupX corretto, con doublet score): VvProt_prefilter.rds\n")
cat("Metadata completo: metadata_pre_filter.csv\n")
cat("Soglie suggerite: soglie_suggerite.csv\n")
cat(sprintf("Contaminazione SoupX: %.1f%%\n", 100 * sc$fit$rhoEst))
cat(sprintf("Doublet predetti: %d / %d (%.1f%%)\n", n_doublet, ncol(obj),
            100 * n_doublet / ncol(obj)))
cat("\nPROSSIMO STEP:\n")
cat("  1. Trasferire plot in locale, valutare soglie\n")
cat("  2. Applicare filtraggio definitivo (script 04_filter.R)\n")
