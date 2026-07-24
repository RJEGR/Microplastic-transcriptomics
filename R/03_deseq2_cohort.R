# =============================================================================
# 03_deseq2_cohort.R
# Modelo DE "B" (limpio): ~cohort sobre 9 muestras 2022 (Ch14, EP08, Se17).
# SES (2024) excluido: confundido con batch.
#
# Se corre para AMBOS niveles: gene (LOC*) y transcript (XM_/XR_).
# Metadata derivada del sample_id (sin Manifest.tsv).
#
# Salidas por nivel:
#   - results/dds_modelB_<level>.rds
#   - results/vsd_modelB_<level>.rds
#   - results/DE_results_modelB_<level>.rds
#   - figures/pca_modelB_<level>.png
#   - figures/ma_<contrast>_<level>.png
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
})

here <- normalizePath(dirname(sys.frame(1)$ofile %||% "R/03_deseq2_cohort.R"),
                      mustWork = FALSE)
if (!file.exists(file.path(here, "config.R"))) here <- "R"
source(file.path(here, "config.R"))
source(file.path(here, "functions.R"))

log_msg("03_deseq2_cohort: iniciando")

# ---- Metadata desde sample_id (no hay Manifest.tsv) ------------------------
build_meta <- function(sample_ids) {
  tibble(sample_id = sample_ids) %>%
    mutate(cohort = CFG$cohort_from_sample(sample_id)) %>%
    filter(cohort %in% CFG$cohorts_keep) %>%
    mutate(cohort = factor(cohort, levels = CFG$cohorts_keep)) %>%
    arrange(cohort, sample_id)
}

# ---- Loop por nivel --------------------------------------------------------
run_deseq_level <- function(level = c("gene","transcript")) {
  level <- match.arg(level)
  counts_path <- switch(level,
    gene       = CFG$gene_counts_csv,
    transcript = CFG$transcript_counts_csv)

  require_file(counts_path,
    "Matriz de conteos (regenerada desde ballgown/*.gtf)")

  log_msg("  === nivel: ", level, " ===")
  log_msg("  counts: ", counts_path)

  counts_raw <- readr::read_csv(counts_path, show_col_types = FALSE)
  colnames(counts_raw)[1] <- ifelse(level == "gene", "gene_id", "transcript_id")

  sample_cols <- colnames(counts_raw)[-1]
  meta <- build_meta(sample_cols)
  log_msg("  muestras retenidas: ", nrow(meta))
  print(dplyr::count(meta, cohort))

  counts_mat <- counts_raw %>%
    column_to_rownames(colnames(counts_raw)[1]) %>%
    as.matrix()
  counts_mat <- counts_mat[, meta$sample_id, drop = FALSE]
  mode(counts_mat) <- "integer"

  log_msg("  matriz: ", nrow(counts_mat), " x ", ncol(counts_mat))

  coldata <- as.data.frame(meta)
  rownames(coldata) <- coldata$sample_id

  dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                                colData   = coldata,
                                design    = CFG$design_formula)
  keep <- rowSums(counts(dds) >= 10) >= 3
  log_msg("  prefiltro: ", sum(keep), " / ", length(keep), " retenidos")
  dds <- dds[keep, ]

  dds <- DESeq(dds)
  backup_write(dds, CFG$out$dds[[level]])

  vsd <- vst(dds, blind = FALSE)
  backup_write(vsd, CFG$out$vsd[[level]])

  # PCA
  pca_df <- plotPCA(vsd, intgroup = "cohort", returnData = TRUE)
  pv     <- round(100 * attr(pca_df, "percentVar"))
  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = cohort, label = name)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
    labs(x = paste0("PC1 (", pv[1], "%)"),
         y = paste0("PC2 (", pv[2], "%)"),
         title = sprintf("PCA - Modelo B (~cohort) - %s-level", level)) +
    theme_bw(base_size = 12)
  ggsave(p_pca, filename = sprintf("pca_modelB_%s.png", level),
         path = CFG$figures_dir, width = 6.5, height = 5, dpi = 300)

  # Contrastes
  all_res <- purrr::map_dfr(names(CFG$contrasts), function(nm) {
    ctr <- CFG$contrasts[[nm]]
    log_msg("  contraste ", nm, ": ", ctr[1], " vs ", ctr[2])
    res <- get_res(dds, ctr, alpha_cutoff = CFG$alpha) %>%
      mutate(
        contrast = nm,
        is_DE    = padj < CFG$alpha & abs(log2FoldChange) > CFG$lfc_thresh,
        direction = case_when(
          is_DE & log2FoldChange >  CFG$lfc_thresh ~ "up",
          is_DE & log2FoldChange < -CFG$lfc_thresh ~ "down",
          TRUE ~ "ns")
      )

    p_ma <- ggplot(res, aes(x = log2(baseMean + 1), y = log2FoldChange,
                            color = direction)) +
      geom_point(alpha = 0.4, size = 0.7) +
      scale_color_manual(values = c(up="#d62728", down="#1f77b4", ns="grey70")) +
      geom_hline(yintercept = c(-CFG$lfc_thresh, CFG$lfc_thresh),
                 linetype = "dashed", color = "grey40") +
      labs(title = paste(nm, "-", level), x = "log2(baseMean+1)",
           y = "log2 Fold Change") +
      theme_bw(base_size = 11)
    ggsave(p_ma, filename = sprintf("ma_%s_%s.png", nm, level),
           path = CFG$figures_dir, width = 5.5, height = 4, dpi = 300)
    res
  })

  log_msg("  DEGs (padj<", CFG$alpha, " & |LFC|>", CFG$lfc_thresh, "):")
  # Print robusto (evita pivot_wider por si direction tiene un solo nivel en algun contraste)
  try({
    smry <- all_res %>%
      dplyr::filter(is_DE) %>%
      dplyr::count(contrast, direction)
    print(as.data.frame(smry))
  }, silent = FALSE)

  backup_write(all_res, CFG$out$de_results[[level]])
}

for (lvl in CFG$levels) run_deseq_level(lvl)

log_msg("03_deseq2_cohort: OK")
