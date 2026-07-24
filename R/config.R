# =============================================================================
# config.R
# Rutas, umbrales y parametros centralizados para el pipeline GO/DE
# Proyecto: Microplastic-transcriptomics (Magallana gigas)
# Autor: Ricardo Gomez-Reyes | 2026
# -----------------------------------------------------------------------------
# EDITA SOLO ESTE ARCHIVO cuando cambien rutas o umbrales.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ---- Operador null-coalescing (R < 4.4 compat) ------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

CFG <- list()
CFG$repo_dir  <- normalizePath(
  file.path(dirname(sys.frame(1)$ofile %||% "R/config.R"), ".."),
  mustWork = FALSE)
if (!dir.exists(CFG$repo_dir)) CFG$repo_dir <- getwd()

# ---- Directorios de trabajo -------------------------------------------------
CFG$R_dir       <- file.path(CFG$repo_dir, "R")
CFG$results_dir <- file.path(CFG$repo_dir, "results")
CFG$figures_dir <- file.path(CFG$repo_dir, "figures")
CFG$backup_dir  <- file.path(CFG$results_dir, "backup")
CFG$logs_dir    <- file.path(CFG$results_dir, "logs")

for (d in c(CFG$results_dir, CFG$figures_dir, CFG$backup_dir, CFG$logs_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Entradas ---------------------------------------------------------------

# 1) Anotacion eggNOG-mapper (proteoma expresado xbMagGiga1.1)
CFG$emapper_file <- file.path(
  CFG$repo_dir, "workdir", "10_emapper",
  "eggnog_xbMagGiga1.1.emapper.annotations")

# 2) GTF merged de StringTie (para mapa gene <-> transcript RefSeq)
#    Producido por 04b_stringtie_merge.slurm con -G xbMagGiga1.1.gtf
#    Atributos utilizados: gene_id, transcript_id, ref_gene_id
CFG$stringtie_gtf <- file.path(
  CFG$repo_dir, "workdir", "04_stringtie", "merged", "merged.gtf")

# 3) Matrices de conteos regeneradas desde ballgown/*.gtf (paso 04c)
#    gene_id = LOC*/MSTRG.*/nombre_gen ; transcript_id = XM_/XR_/NM_/NR_/MSTRG.*.*
CFG$gene_counts_csv       <- file.path(
  CFG$repo_dir, "workdir", "05_matrix_pca", "gene_count_matrix.csv")
CFG$transcript_counts_csv <- file.path(
  CFG$repo_dir, "workdir", "05_matrix_pca", "transcript_count_matrix.csv")

# 4) Metadata: NO hay Manifest.tsv. Se deriva la cohorte del sample_id via regex.
#    Ch14_1..3 -> Ch14 ; EP08_1..3 -> EP08 ; Se17_1..3 -> Se17 ; SES{1..3} -> SES
CFG$cohort_from_sample <- function(sid) {
  sub("_?[0-9]+$", "", sid)
}

# ---- Diseno experimental ----------------------------------------------------
# Modelo B "limpio": ~cohort sobre las 9 muestras 2022 (Ch14, EP08, Se17).
# SES excluido (batch confounded con 2024).
CFG$cohorts_keep   <- c("Ch14", "EP08", "Se17")
CFG$design_formula <- ~ cohort

# Contrastes pareados: log2FC = num vs den
CFG$contrasts <- list(
  EP08_vs_Ch14 = c("EP08", "Ch14"),
  Se17_vs_Ch14 = c("Se17", "Ch14"),
  Se17_vs_EP08 = c("Se17", "EP08")
)

# ---- Niveles de analisis ----------------------------------------------------
# Ambos: gene-level (LOC*) y transcript-level (XM_*). Se corren en paralelo.
CFG$levels <- c("gene", "transcript")

# ---- Umbrales de significancia ---------------------------------------------
CFG$alpha       <- 0.05    # padj
CFG$lfc_thresh  <- 1       # |log2FoldChange|

# ---- topGO -----------------------------------------------------------------
CFG$ontologies         <- c("BP", "MF", "CC")
CFG$topgo_nodes        <- 50
CFG$topgo_conservative <- TRUE

# ---- rrvgo -----------------------------------------------------------------
# NOTA: Magallana gigas no tiene orgdb en Bioconductor.
# org.Hs.eg.db se usa SOLO como fuente de la estructura semantica de GO.
CFG$rrvgo_orgdb     <- "org.Hs.eg.db"
CFG$rrvgo_threshold <- 0.9

# ---- Backup ----------------------------------------------------------------
CFG$stamp <- format(Sys.time(), "%Y-%m-%d")

# ---- Ficheros de salida (por nivel) ----------------------------------------
CFG$out <- list(
  emapper_tidy = file.path(CFG$results_dir, "emapper_annotations_tidy.rds"),
  # gene2GO por xm_id (transcrito RefSeq) - key = XM_...
  gene2GO_xm   = list(
    BP = file.path(CFG$results_dir, "gene2GO_XM_BP.rds"),
    MF = file.path(CFG$results_dir, "gene2GO_XM_MF.rds"),
    CC = file.path(CFG$results_dir, "gene2GO_XM_CC.rds")
  ),
  # gene2GO por LOC (agregado por gen RefSeq) - key = LOC...
  gene2GO_loc  = list(
    BP = file.path(CFG$results_dir, "gene2GO_LOC_BP.rds"),
    MF = file.path(CFG$results_dir, "gene2GO_LOC_MF.rds"),
    CC = file.path(CFG$results_dir, "gene2GO_LOC_CC.rds")
  ),
  # Mapa transcript_id <-> gene_id extraido de merged.gtf
  gtf_map      = file.path(CFG$results_dir, "gtf_tx2gene_map.rds"),
  # DESeq2 por nivel
  dds          = list(
    gene       = file.path(CFG$results_dir, "dds_modelB_gene.rds"),
    transcript = file.path(CFG$results_dir, "dds_modelB_transcript.rds")
  ),
  vsd          = list(
    gene       = file.path(CFG$results_dir, "vsd_modelB_gene.rds"),
    transcript = file.path(CFG$results_dir, "vsd_modelB_transcript.rds")
  ),
  de_results   = list(
    gene       = file.path(CFG$results_dir, "DE_results_modelB_gene.rds"),
    transcript = file.path(CFG$results_dir, "DE_results_modelB_transcript.rds")
  ),
  topgo_results = list(
    gene       = file.path(CFG$results_dir, "topGO_results_modelB_gene.rds"),
    transcript = file.path(CFG$results_dir, "topGO_results_modelB_transcript.rds")
  ),
  rrvgo_results = list(
    gene       = file.path(CFG$results_dir, "rrvgo_semantic_modelB_gene.rds"),
    transcript = file.path(CFG$results_dir, "rrvgo_semantic_modelB_transcript.rds")
  )
)

# ---- Helpers ---------------------------------------------------------------
require_file <- function(path, hint = "") {
  if (!file.exists(path)) {
    stop(sprintf(
      "\n[config] Falta archivo requerido:\n  %s\n  %s\n  Edita R/config.R.",
      path, hint), call. = FALSE)
  }
  invisible(TRUE)
}

backup_write <- function(obj, path) {
  if (file.exists(path)) {
    bdir <- file.path(CFG$backup_dir, CFG$stamp)
    dir.create(bdir, recursive = TRUE, showWarnings = FALSE)
    file.copy(path, file.path(bdir, basename(path)), overwrite = TRUE)
  }
  readr::write_rds(obj, path)
  invisible(path)
}

message("[config] repo_dir = ", CFG$repo_dir)
