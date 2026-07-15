#!/usr/bin/env Rscript
# =============================================================================
#  deseq2_pipeline.R
#  -------------------------------------------------------------------------
#  Analisis de Expresion Diferencial (DE) gen-nivel con DESeq2 sobre las
#  12 librerias paired-end de Microplastic-transcriptomics (StringTie ->
#  prepDE_counts.py -> gene_count_matrix.csv).
#
#  Diseno experimental:
#    Ch14 (Chapala 2022, alta invierno)   x3
#    EP08 (Punta Banda 2022, alta)        x3
#    Se17 (Sesma 2022, baja)              x3
#    SES  (Sesma 2024, lote nuevo)        x3   <-- batch confundido con sitio
#
#  Decision de diseno:
#    El lote 2024 (SES) esta perfectamente confundido con la cohorte SES,
#    asi que NO se puede modelar `batch` como covariable (rango deficiente).
#    Estrategia:
#      (1) PCA QC con todas las muestras
#      (2) DE entre cohortes 2022 (Ch14/EP08/Se17)  -> diseno ~ cohort
#      (3) SES vs Se17 reportado APARTE con caveat de batch
#
#  Uso:
#    Rscript deseq2_pipeline.R [COUNTS_CSV] [META_TSV] [OUT_DIR]
# =============================================================================

suppressPackageStartupMessages({
  required <- c("DESeq2", "apeglm", "EnhancedVolcano", "pheatmap",
                "ggplot2", "dplyr", "RColorBrewer", "matrixStats")
  bioc <- c("DESeq2", "apeglm", "EnhancedVolcano")
  to_install <- setdiff(required, rownames(installed.packages()))
  if (length(to_install)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(intersect(to_install, bioc), update = FALSE, ask = FALSE)
    cran <- setdiff(to_install, bioc)
    if (length(cran)) install.packages(cran, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(required, library, character.only = TRUE))
})

# ---- 1. Configuracion --------------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
counts_csv  <- if (length(args) >= 1) args[1] else "gene_count_matrix.csv"
meta_tsv    <- if (length(args) >= 2) args[2] else "sample_metadata.tsv"
out_dir     <- if (length(args) >= 3) args[3] else "deseq2_out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(2026)
alpha_fdr   <- 0.05
lfc_thresh  <- 1.0
sample_ord  <- c("Ch14_1","Ch14_2","Ch14_3",
                 "EP08_1","EP08_2","EP08_3",
                 "Se17_1","Se17_2","Se17_3",
                 "SES1",  "SES2",  "SES3")

# ---- 2. Carga de conteos y metadatos ----------------------------------------
counts <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts[, sample_ord])
storage.mode(counts) <- "integer"

coldata <- read.delim(meta_tsv, row.names = "sample")
coldata <- coldata[sample_ord, ]
coldata$cohort <- factor(coldata$cohort,
                         levels = c("Ch14","EP08","Se17","SES"))
coldata$batch  <- factor(coldata$batch, levels = c("2022","2024"))

stopifnot(identical(colnames(counts), rownames(coldata)))
cat(sprintf("[INFO] %d genes x %d muestras\n", nrow(counts), ncol(counts)))

# =============================================================================
# (A) QC GLOBAL - VST + PCA con las 12 muestras
# =============================================================================
dds_all <- DESeqDataSetFromMatrix(counts, coldata, design = ~ cohort)
keep    <- rowSums(counts(dds_all)) >= 10
dds_all <- dds_all[keep, ]
cat(sprintf("[INFO] Genes tras pre-filtro (suma >=10): %d\n", nrow(dds_all)))

vsd_all <- vst(dds_all, blind = TRUE)

# PCA con etiqueta de cohorte + batch
pca_df <- plotPCA(vsd_all, intgroup = c("cohort","batch"), returnData = TRUE,
                  ntop = 5000)
pv <- round(100 * attr(pca_df, "percentVar"))
pal <- c(Ch14="#1f77b4", EP08="#d62728", Se17="#2ca02c", SES="#9467bd")
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = cohort, shape = batch)) +
  geom_hline(yintercept = 0, lty = 2, color = "grey80") +
  geom_vline(xintercept = 0, lty = 2, color = "grey80") +
  geom_point(size = 4, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = rownames(pca_df)),
                           size = 3, show.legend = FALSE) +
  scale_color_manual(values = pal, name = "Cohorte") +
  scale_shape_manual(values = c("2022" = 16, "2024" = 17), name = "Lote") +
  labs(title = "PCA - VST (DESeq2) | 12 librerias",
       x = sprintf("PC1 (%d%%)", pv[1]),
       y = sprintf("PC2 (%d%%)", pv[2])) +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "QC_PCA_vst.pdf"), p_pca, width = 8, height = 6)
ggsave(file.path(out_dir, "QC_PCA_vst.png"), p_pca, width = 8, height = 6, dpi = 300)

# Heatmap de distancias entre muestras
sampleDists <- dist(t(assay(vsd_all)))
distMat     <- as.matrix(sampleDists)
pheatmap(distMat,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(RColorBrewer::brewer.pal(9,"Blues")))(255),
         annotation_col = data.frame(cohort = coldata$cohort,
                                     batch  = coldata$batch,
                                     row.names = rownames(coldata)),
         filename = file.path(out_dir, "QC_sample_distances.pdf"),
         width = 8, height = 6.5)

# =============================================================================
# (B) MODELO PRINCIPAL - solo cohortes de 2022 (clean biological comparison)
# =============================================================================
keep_2022    <- coldata$batch == "2022"
counts_22    <- counts[, keep_2022]
coldata_22   <- droplevels(coldata[keep_2022, ])

dds <- DESeqDataSetFromMatrix(counts_22, coldata_22, design = ~ cohort)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
cat("[INFO] resultsNames(dds):\n"); print(resultsNames(dds))

# Plot de dispersion
pdf(file.path(out_dir, "QC_dispersion_2022.pdf"), width = 7, height = 6)
plotDispEsts(dds, main = "Dispersion estimates (cohortes 2022)")
dev.off()

# Guardar objeto fitted para iteracion downstream
saveRDS(dds, file.path(out_dir, "dds_2022_fitted.rds"))

# ---- contrastes pareados ----------------------------------------------------
contrasts_list <- list(
  EP08_vs_Ch14 = c("cohort", "EP08", "Ch14"),
  Se17_vs_Ch14 = c("cohort", "Se17", "Ch14"),
  Se17_vs_EP08 = c("cohort", "Se17", "EP08")
)

run_contrast <- function(nm, ctr) {
  res <- results(dds, contrast = ctr, alpha = alpha_fdr,
                 independentFilter = TRUE)
  # apeglm requires coef name; cuando es contrast generico usamos ashr
  res_sh <- lfcShrink(dds, contrast = ctr, type = "ashr", res = res)
  df  <- as.data.frame(res_sh)
  df  <- df[order(df$padj, na.last = TRUE), ]
  df$gene_id <- rownames(df)
  write.csv(df, file.path(out_dir, sprintf("DE_%s_all.csv", nm)),
            row.names = FALSE)
  sig <- df[!is.na(df$padj) & df$padj < alpha_fdr &
            abs(df$log2FoldChange) > lfc_thresh, ]
  write.csv(sig, file.path(out_dir, sprintf("DE_%s_significant.csv", nm)),
            row.names = FALSE)

  # Volcano
  EnhancedVolcano::EnhancedVolcano(df,
      lab = df$gene_id, x = "log2FoldChange", y = "padj",
      pCutoff = alpha_fdr, FCcutoff = lfc_thresh,
      title = nm,
      subtitle = sprintf("DESeq2 + ashr | %d significativos | |LFC|>%g, padj<%g",
                         nrow(sig), lfc_thresh, alpha_fdr),
      legendPosition = "right")
  ggsave(file.path(out_dir, sprintf("volcano_%s.pdf", nm)),
         width = 9, height = 7)

  # ranked list para GSEA
  rk <- df[!is.na(df$pvalue), ]
  rk$score <- sign(rk$log2FoldChange) * -log10(rk$pvalue + 1e-300)
  rk <- rk[order(rk$score, decreasing = TRUE), c("gene_id","score")]
  write.table(rk, file.path(out_dir, sprintf("ranked_%s.rnk", nm)),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  cat(sprintf("  %-15s -> total=%d  sig=%d  up=%d  down=%d\n",
              nm, nrow(df), nrow(sig),
              sum(sig$log2FoldChange > 0), sum(sig$log2FoldChange < 0)))
  invisible(df)
}

cat("\n[INFO] Ejecutando contrastes pareados (2022):\n")
res_list <- mapply(run_contrast, names(contrasts_list), contrasts_list,
                   SIMPLIFY = FALSE)

# Heatmap de la union de top-30 DEGs por contraste
top_genes <- unique(unlist(lapply(res_list, function(d) {
  d <- d[!is.na(d$padj), ]
  head(d$gene_id[order(d$padj)], 30)
})))
mat <- assay(vst(dds, blind = FALSE))[top_genes, ]
mat <- mat - rowMeans(mat)
ann <- data.frame(cohort = coldata_22$cohort,
                  row.names = rownames(coldata_22))
pheatmap(mat,
         annotation_col = ann,
         annotation_colors = list(cohort = pal[levels(coldata_22$cohort)]),
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, fontsize_col = 10,
         filename = file.path(out_dir, "heatmap_top_DEGs_2022.pdf"),
         width = 7, height = 9)

# =============================================================================
# (C) Contraste SES vs Se17 - mismo sitio, distinto lote (CAVEAT)
# =============================================================================
keep_sesma  <- coldata$cohort %in% c("Se17","SES")
coldata_ses <- droplevels(coldata[keep_sesma, ])
counts_ses  <- counts[, keep_sesma]
dds_ses <- DESeqDataSetFromMatrix(counts_ses, coldata_ses, design = ~ cohort)
dds_ses <- dds_ses[rowSums(counts(dds_ses)) >= 10, ]
dds_ses <- DESeq(dds_ses)
res_ses <- results(dds_ses, contrast = c("cohort","SES","Se17"),
                   alpha = alpha_fdr)
res_ses_sh <- lfcShrink(dds_ses, contrast = c("cohort","SES","Se17"),
                        type = "ashr", res = res_ses)
df_ses <- as.data.frame(res_ses_sh)
df_ses$gene_id <- rownames(df_ses)
df_ses <- df_ses[order(df_ses$padj, na.last = TRUE), ]
write.csv(df_ses, file.path(out_dir, "DE_SES_vs_Se17_BATCH_CONFOUNDED.csv"),
          row.names = FALSE)
sig_ses <- df_ses[!is.na(df_ses$padj) & df_ses$padj < alpha_fdr &
                  abs(df_ses$log2FoldChange) > lfc_thresh, ]
cat(sprintf("\n[CAVEAT] SES vs Se17 (lote 2024 vs 2022): %d DEGs - efecto biologico ininterpretable sin replicacion entre lotes.\n",
            nrow(sig_ses)))

# =============================================================================
# (D) Resumen
# =============================================================================
summary_df <- data.frame(
  contrast = c(names(contrasts_list), "SES_vs_Se17_BATCH"),
  n_DEG    = c(sapply(res_list, function(d) sum(!is.na(d$padj) &
                       d$padj < alpha_fdr & abs(d$log2FoldChange) > lfc_thresh)),
               nrow(sig_ses)),
  notas    = c(rep("clean 2022", 3), "batch confundido - solo exploratorio")
)
write.csv(summary_df, file.path(out_dir, "DE_summary.csv"), row.names = FALSE)
print(summary_df)
cat("\n[DONE]\n")
