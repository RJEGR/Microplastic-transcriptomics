#!/usr/bin/env Rscript
# =============================================================================
#  build_matrix_pca.R
#  -------------------------------------------------------------------------
#  Construye la matriz de expresion gen-nivel (TPM y FPKM) a partir de los
#  archivos *.gene_abund.tab generados por StringTie (-A) para las 12
#  librerias paired-end del proyecto Microplastic-transcriptomics, y ejecuta
#  un Analisis de Componentes Principales (PCA) etiquetado por cohorte.
#
#  Cohortes biologicas:
#    - Ch14 (Chapala, alta - invierno 2022)        : Ch14_1, Ch14_2, Ch14_3
#    - EP08 (Punta Banda, alta - 2022)             : EP08_1, EP08_2, EP08_3
#    - Se17 (Sesma,    baja  - 2022)               : Se17_1, Se17_2, Se17_3
#    - SES  (Sesma,    lote 2024 sin clasificar)   : SES1,  SES2,  SES3
#
#  Autor : Bioinformatica - UABC (rgomez41@uabc.edu.mx)
#  Fecha : 2026-06-24
#  Uso   : Rscript build_matrix_pca.R [STRINGTIE_DIR] [OUT_DIR]
# =============================================================================

suppressPackageStartupMessages({
  required <- c("data.table", "ggplot2", "ggrepel", "matrixStats")
  to_install <- required[!required %in% rownames(installed.packages())]
  if (length(to_install)) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(matrixStats)
})

# ---- 1. Configuracion --------------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
stringtie_d <- if (length(args) >= 1) args[1] else "../04_stringtie/per_sample"
out_dir     <- if (length(args) >= 2) args[2] else "."
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Orden canonico de muestras y cohortes (1 fila por libreria)
sample_meta <- data.table(
  sample = c("Ch14_1","Ch14_2","Ch14_3",
             "EP08_1","EP08_2","EP08_3",
             "Se17_1","Se17_2","Se17_3",
             "SES1",  "SES2",  "SES3"),
  cohort = c(rep("Ch14", 3), rep("EP08", 3), rep("Se17", 3), rep("SES", 3)),
  origin = c(rep("Chapala (alta, invierno 2022)", 3),
             rep("Punta Banda (alta, 2022)",       3),
             rep("Sesma (baja, 2022)",             3),
             rep("Sesma (lote 2024, sin clasif.)", 3)),
  batch  = c(rep("2022", 9), rep("2024", 3))
)
fwrite(sample_meta,
       file.path(out_dir, "sample_metadata.tsv"),
       sep = "\t")

cat(sprintf("[INFO] %d librerias | %d cohortes\n",
            nrow(sample_meta), length(unique(sample_meta$cohort))))

# ---- 2. Lectura modular de gene_abund.tab ------------------------------------
read_gene_abund <- function(sample, dir) {
  f <- file.path(dir, paste0(sample, ".gene_abund.tab"))
  if (!file.exists(f)) stop(sprintf("Archivo faltante: %s", f))
  dt <- fread(f, sep = "\t", header = TRUE,
              colClasses = c("character","character","character","character",
                             "integer","integer","numeric","numeric","numeric"))
  # Algunas filas pueden repetir Gene ID (transcritos en multiples loci).
  # Colapsar por Gene ID sumando expresion y conservando 1 etiqueta de nombre/ref.
  dt <- dt[, .(`Gene Name` = `Gene Name`[1],
               Reference   = paste(unique(Reference), collapse = ","),
               Coverage    = sum(Coverage,  na.rm = TRUE),
               FPKM        = sum(FPKM,      na.rm = TRUE),
               TPM         = sum(TPM,       na.rm = TRUE)),
           by = `Gene ID`]
  setnames(dt,
           c("Coverage","FPKM","TPM"),
           paste0(sample, c("__Coverage","__FPKM","__TPM")))
  dt
}

abund_list <- lapply(sample_meta$sample,
                     read_gene_abund, dir = stringtie_d)

# ---- 3. Union outer por Gene ID ---------------------------------------------
# Mantener union completa: genes ausentes en alguna muestra -> 0
union_dt <- Reduce(function(x, y) merge(x, y,
                                        by      = "Gene ID",
                                        all     = TRUE,
                                        suffixes = c("", ".y")),
                   abund_list)

# Reconciliar columnas Gene Name / Reference (puede haber .y)
name_cols <- grep("^Gene Name", names(union_dt), value = TRUE)
ref_cols  <- grep("^Reference",  names(union_dt), value = TRUE)
union_dt[, `Gene Name` := do.call(coalesce_chr <- function(...) {
  v <- c(...); v <- v[!is.na(v) & v != ""]
  if (length(v)) v[1] else NA_character_
}, .SD), .SDcols = name_cols, by = `Gene ID`]
union_dt[, Reference  := do.call(coalesce_chr, .SD),
         .SDcols = ref_cols, by = `Gene ID`]
union_dt[, (setdiff(c(name_cols, ref_cols), c("Gene Name","Reference"))) := NULL]

# Reemplazar NA por 0 en columnas numericas
num_cols <- grep("__(Coverage|FPKM|TPM)$", names(union_dt), value = TRUE)
for (col in num_cols) set(union_dt, which(is.na(union_dt[[col]])), col, 0)

# ---- 4. Matrices finales (filas = genes, columnas = muestras) ---------------
make_matrix <- function(metric) {
  cols <- paste0(sample_meta$sample, "__", metric)
  m    <- union_dt[, c("Gene ID","Gene Name","Reference", ..cols)]
  setnames(m, cols, sample_meta$sample)
  setcolorder(m, c("Gene ID","Gene Name","Reference", sample_meta$sample))
  m
}

mat_tpm      <- make_matrix("TPM")
mat_fpkm     <- make_matrix("FPKM")
mat_coverage <- make_matrix("Coverage")

fwrite(mat_tpm,      file.path(out_dir, "gene_matrix_TPM.tsv"),      sep = "\t")
fwrite(mat_fpkm,     file.path(out_dir, "gene_matrix_FPKM.tsv"),     sep = "\t")
fwrite(mat_coverage, file.path(out_dir, "gene_matrix_Coverage.tsv"), sep = "\t")

cat(sprintf("[INFO] Matriz TPM: %d genes x %d muestras\n",
            nrow(mat_tpm), length(sample_meta$sample)))

# ---- 5. Pre-procesamiento para PCA ------------------------------------------
# Usamos TPM por estar normalizada a longitud y profundidad, log2(TPM+1)
# para estabilizar varianza, y filtramos genes de baja senal.
expr <- as.matrix(mat_tpm[, sample_meta$sample, with = FALSE])
rownames(expr) <- mat_tpm$`Gene ID`

# Filtro: mantener genes con TPM >= 1 en al menos 3 muestras (n minimo por cohorte)
keep   <- rowSums(expr >= 1) >= 3
expr_f <- log2(expr[keep, ] + 1)

cat(sprintf("[INFO] Genes tras filtro (TPM>=1 en >=3 muestras): %d\n",
            sum(keep)))

# Seleccionar top-N genes mas variables (mejora separacion biologica vs ruido)
top_n  <- min(5000, nrow(expr_f))
vars   <- matrixStats::rowVars(expr_f)
expr_f <- expr_f[order(vars, decreasing = TRUE)[seq_len(top_n)], ]

# ---- 6. PCA -----------------------------------------------------------------
pca <- prcomp(t(expr_f), center = TRUE, scale. = TRUE)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2) * 100

scores <- as.data.table(pca$x[, 1:4, drop = FALSE], keep.rownames = "sample")
scores <- merge(scores, sample_meta, by = "sample", sort = FALSE)
fwrite(scores, file.path(out_dir, "pca_scores.tsv"), sep = "\t")

loadings <- as.data.table(pca$rotation[, 1:4, drop = FALSE],
                          keep.rownames = "Gene ID")
fwrite(loadings, file.path(out_dir, "pca_loadings.tsv"), sep = "\t")

# ---- 7. Visualizacion PC1 vs PC2 --------------------------------------------
cohort_palette <- c("Ch14" = "#1f77b4",
                    "EP08" = "#d62728",
                    "Se17" = "#2ca02c",
                    "SES"  = "#9467bd")

p_main <- ggplot(scores, aes(PC1, PC2, color = cohort, shape = batch)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(aes(label = sample), size = 3.2,
                  box.padding = 0.4, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = cohort_palette,
                     name = "Cohorte",
                     labels = c("Ch14" = "Ch14 - Chapala",
                                "EP08" = "EP08 - Punta Banda",
                                "Se17" = "Se17 - Sesma 2022",
                                "SES"  = "SES  - Sesma 2024")) +
  scale_shape_manual(values = c("2022" = 16, "2024" = 17), name = "Lote") +
  labs(title = "PCA - expresion genica (StringTie TPM, log2)",
       subtitle = sprintf("Top %d genes mas variables | n = %d librerias",
                          top_n, ncol(expr_f)),
       x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey30"))

ggsave(file.path(out_dir, "PCA_PC1_PC2.png"),
       p_main, width = 8.5, height = 6, dpi = 300)
ggsave(file.path(out_dir, "PCA_PC1_PC2.pdf"),
       p_main, width = 8.5, height = 6)

# Scree plot
scree <- data.table(PC = factor(paste0("PC", seq_along(var_exp)),
                                levels = paste0("PC", seq_along(var_exp))),
                    var = var_exp)[1:min(10, length(var_exp))]
p_scree <- ggplot(scree, aes(PC, var)) +
  geom_col(fill = "#4477AA") +
  geom_text(aes(label = sprintf("%.1f%%", var)),
            vjust = -0.4, size = 3) +
  labs(title = "Varianza explicada por componente",
       x = NULL, y = "% varianza") +
  theme_bw(base_size = 12) +
  ylim(0, max(scree$var) * 1.15)

ggsave(file.path(out_dir, "PCA_scree.png"),
       p_scree, width = 7, height = 4.5, dpi = 300)

# ---- 8. Resumen en consola --------------------------------------------------
cat("\n[RESUMEN PCA]\n")
cat(sprintf("  PC1: %.2f%%   PC2: %.2f%%   PC3: %.2f%%\n",
            var_exp[1], var_exp[2], var_exp[3]))
cat("\n[OUTPUTS]\n")
cat(sprintf("  - %s\n", normalizePath(out_dir)))
cat("    gene_matrix_TPM.tsv / FPKM.tsv / Coverage.tsv\n")
cat("    sample_metadata.tsv\n")
cat("    pca_scores.tsv  pca_loadings.tsv\n")
cat("    PCA_PC1_PC2.png/.pdf  PCA_scree.png\n")
