##=============================================================================
## FASE 3+4 (parte 2): Filtraggio — applica soglie QC e rimuove doublet
##
## Input:  03_scRNA/02_qc/VvProt_prefilter.rds (SoupX corretto, con doublet score)
## Output: 03_scRNA/03_filtering/VvProt_filtered.rds (pronto per normalizzazione)
##
## Soglie decise dopo valutazione dei plot QC:
##   nCount_RNA:       [1313 — 57924]   (3 MAD su log10)
##   nFeature_RNA:     [748 — 8568]     (3 MAD su log10)
##   log10GenesPerUMI: > 0.80
##   scDblFinder:      class == "singlet"
##
## Eseguire su TeraStat2: sbatch 04_filter.slurm
##=============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)

set.seed(42)

# --- Paths ----
proj    <- "/lustre/home/dcapozzi/projects/scRNA_protoplasts"
indir   <- file.path(proj, "03_scRNA/02_qc")
outdir  <- file.path(proj, "03_scRNA/03_filtering")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- Soglie ----
NCOUNT_MIN   <- 1313
NCOUNT_MAX   <- 57924
NFEATURE_MIN <- 748
NFEATURE_MAX <- 8568
COMPLEXITY   <- 0.80


##=============================================================================
## 1. CARICA OGGETTO PRE-FILTRAGGIO
##=============================================================================

cat("=== 1. Caricamento oggetto pre-filtraggio ===\n")
obj <- readRDS(file.path(indir, "VvProt_prefilter.rds"))
cat(sprintf("Cellule pre-filtraggio: %d\n", ncol(obj)))
cat(sprintf("Geni: %d\n", nrow(obj)))


##=============================================================================
## 2. APPLICA FILTRI
##=============================================================================

cat("\n=== 2. Applicazione filtri ===\n")

# Etichetta ogni cellula con il motivo di esclusione
obj[["pass_nCount"]]   <- obj$nCount_RNA >= NCOUNT_MIN & obj$nCount_RNA <= NCOUNT_MAX
obj[["pass_nFeature"]] <- obj$nFeature_RNA >= NFEATURE_MIN & obj$nFeature_RNA <= NFEATURE_MAX
obj[["pass_complex"]]  <- obj$log10GenesPerUMI > COMPLEXITY
obj[["pass_singlet"]]  <- obj$scDblFinder.class == "singlet"
obj[["pass_all"]]      <- obj$pass_nCount & obj$pass_nFeature & obj$pass_complex & obj$pass_singlet

# Report dettagliato
cat(sprintf("  Escluse per nCount < %d:        %d\n", NCOUNT_MIN, sum(obj$nCount_RNA < NCOUNT_MIN)))
cat(sprintf("  Escluse per nCount > %d:       %d\n", NCOUNT_MAX, sum(obj$nCount_RNA > NCOUNT_MAX)))
cat(sprintf("  Escluse per nFeature < %d:       %d\n", NFEATURE_MIN, sum(obj$nFeature_RNA < NFEATURE_MIN)))
cat(sprintf("  Escluse per nFeature > %d:      %d\n", NFEATURE_MAX, sum(obj$nFeature_RNA > NFEATURE_MAX)))
cat(sprintf("  Escluse per complessita' < %.2f:  %d\n", COMPLEXITY, sum(obj$log10GenesPerUMI <= COMPLEXITY)))
cat(sprintf("  Escluse come doublet:            %d\n", sum(obj$scDblFinder.class == "doublet")))
cat(sprintf("  Cellule che passano TUTTI:       %d / %d (%.1f%%)\n",
            sum(obj$pass_all), ncol(obj), 100 * mean(obj$pass_all)))

# Filtra
obj_filtered <- subset(obj, pass_all == TRUE)
cat(sprintf("\nCellule dopo filtraggio: %d\n", ncol(obj_filtered)))
cat(sprintf("Geni: %d\n", nrow(obj_filtered)))

# Rimuovi geni espressi in 0 cellule dopo il filtraggio
counts_filt <- GetAssayData(obj_filtered, layer = "counts")
genes_keep <- rownames(counts_filt)[Matrix::rowSums(counts_filt > 0) >= 3]
obj_filtered <- subset(obj_filtered, features = genes_keep)
cat(sprintf("Geni dopo rimozione (min 3 cellule): %d\n", nrow(obj_filtered)))


##=============================================================================
## 3. PLOT PRIMA/DOPO FILTRAGGIO
##=============================================================================

cat("\n=== 3. Plot confronto prima/dopo ===\n")

dpi <- 300

# 3.1 Densita' nCount prima/dopo
pre_df  <- data.frame(nCount = obj$nCount_RNA, stato = "Pre-filtraggio")
post_df <- data.frame(nCount = obj_filtered$nCount_RNA, stato = "Post-filtraggio")
df_count <- rbind(pre_df, post_df)

p1 <- ggplot(df_count, aes(x = nCount, fill = stato)) +
  geom_density(alpha = 0.4) +
  scale_x_log10() +
  scale_fill_manual(values = c("Pre-filtraggio" = "grey60", "Post-filtraggio" = "steelblue")) +
  geom_vline(xintercept = c(NCOUNT_MIN, NCOUNT_MAX), linetype = "dashed", color = "red") +
  labs(title = "Distribuzione UMI — prima vs dopo filtraggio",
       x = "nCount_RNA (log10)", y = "Densita'", fill = "") +
  theme_minimal()
ggsave(file.path(outdir, "01_density_nCount_pre_post.png"), p1,
       width = 8, height = 5, dpi = dpi)

# 3.2 Densita' nFeature prima/dopo
pre_df2  <- data.frame(nFeature = obj$nFeature_RNA, stato = "Pre-filtraggio")
post_df2 <- data.frame(nFeature = obj_filtered$nFeature_RNA, stato = "Post-filtraggio")
df_feat <- rbind(pre_df2, post_df2)

p2 <- ggplot(df_feat, aes(x = nFeature, fill = stato)) +
  geom_density(alpha = 0.4) +
  scale_x_log10() +
  scale_fill_manual(values = c("Pre-filtraggio" = "grey60", "Post-filtraggio" = "coral")) +
  geom_vline(xintercept = c(NFEATURE_MIN, NFEATURE_MAX), linetype = "dashed", color = "red") +
  labs(title = "Distribuzione geni — prima vs dopo filtraggio",
       x = "nFeature_RNA (log10)", y = "Densita'", fill = "") +
  theme_minimal()
ggsave(file.path(outdir, "02_density_nFeature_pre_post.png"), p2,
       width = 8, height = 5, dpi = dpi)

# 3.3 Scatter post-filtraggio
p3 <- ggplot(obj_filtered@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(size = 0.2, alpha = 0.2, color = "steelblue") +
  scale_x_log10() + scale_y_log10() +
  labs(title = sprintf("nCount vs nFeature post-filtraggio (%d cellule)", ncol(obj_filtered)),
       x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
  theme_minimal()
ggsave(file.path(outdir, "03_scatter_post_filter.png"), p3,
       width = 7, height = 6, dpi = dpi)

# 3.4 Violin post-filtraggio
p4 <- VlnPlot(obj_filtered, features = c("nFeature_RNA", "nCount_RNA"),
              ncol = 2, pt.size = 0) & theme_minimal() & NoLegend()
ggsave(file.path(outdir, "04_violin_post_filter.png"), p4,
       width = 10, height = 5, dpi = dpi)


##=============================================================================
## 4. STATISTICHE POST-FILTRAGGIO
##=============================================================================

cat("\n=== STATISTICHE POST-FILTRAGGIO ===\n")
cat(sprintf("Cellule: %d (da %d, rimosso %.1f%%)\n",
            ncol(obj_filtered), ncol(obj),
            100 * (1 - ncol(obj_filtered) / ncol(obj))))
cat(sprintf("Geni: %d\n", nrow(obj_filtered)))

qc_summary <- function(x, name) {
  cat(sprintf("  %-20s min=%-8.1f Q1=%-8.1f median=%-8.1f Q3=%-8.1f max=%-8.1f\n",
              name, min(x), quantile(x, 0.25), median(x), quantile(x, 0.75), max(x)))
}

qc_summary(obj_filtered$nCount_RNA, "nCount_RNA")
qc_summary(obj_filtered$nFeature_RNA, "nFeature_RNA")
qc_summary(obj_filtered$log10GenesPerUMI, "log10GenesPerUMI")


##=============================================================================
## 5. SALVA OGGETTO FILTRATO
##=============================================================================

# Rimuovi colonne temporanee
obj_filtered$pass_nCount   <- NULL
obj_filtered$pass_nFeature <- NULL
obj_filtered$pass_complex  <- NULL
obj_filtered$pass_singlet  <- NULL
obj_filtered$pass_all      <- NULL

saveRDS(obj_filtered, file.path(outdir, "VvProt_filtered.rds"))
write.csv(obj_filtered@meta.data, file.path(outdir, "metadata_post_filter.csv"), row.names = TRUE)

cat("\n=== FILTRAGGIO COMPLETATO ===\n")
cat("Oggetto filtrato: VvProt_filtered.rds\n")
cat("Metadata: metadata_post_filter.csv\n")
cat("Plot in:", outdir, "\n")
cat("\nPROSSIMO STEP: FASE 5 — SCTransform, PCA, clustering, UMAP\n")
