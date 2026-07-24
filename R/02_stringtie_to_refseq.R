# =============================================================================
# 02_stringtie_to_refseq.R
# Del merged.gtf (paso 04b) construye el mapa
#   transcript_id (XM_/XR_/NM_/NR_/MSTRG.*.*) <-> gene_id (LOC*|... / MSTRG.*)
# usando la misma logica que prepDE_counts.py (gene_id|gene_name).
#
# Con ese mapa, colapsa gene2GO XM_ (paso 01) a nivel LOC (gen RefSeq).
#
# Genera:
#   - results/gtf_tx2gene_map.rds
#   - results/gene2GO_LOC_{BP,MF,CC}.rds
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

here <- normalizePath(dirname(sys.frame(1)$ofile %||% "R/02_stringtie_to_refseq.R"),
                      mustWork = FALSE)
if (!file.exists(file.path(here, "config.R"))) here <- "R"
source(file.path(here, "config.R"))
source(file.path(here, "functions.R"))

log_msg("02_stringtie_to_refseq: iniciando")

# ---- 1) Parsea merged.gtf --------------------------------------------------
gtf_df <- parse_stringtie_gtf(CFG$stringtie_gtf)
tx2g   <- build_tx2gene_from_gtf(gtf_df)

# Diagnosticos
n_tx      <- nrow(tx2g)
n_refseq  <- sum(grepl("^[XN][MR]_", tx2g$transcript_id))
n_novel   <- sum(grepl("^MSTRG\\.",  tx2g$transcript_id))
n_genes   <- dplyr::n_distinct(tx2g$gene_id)

log_msg("  transcritos totales   : ", n_tx)
log_msg("  transcritos RefSeq XM_: ", n_refseq,
        sprintf(" (%.1f%%)", 100 * n_refseq / max(1, n_tx)))
log_msg("  transcritos novel (MSTRG): ", n_novel)
log_msg("  gene_id unicos        : ", n_genes)

backup_write(tx2g, CFG$out$gtf_map)
log_msg("  guardado mapa: ", CFG$out$gtf_map)

# ---- 2) Colapsa emapper GOs a gene-level (LOC) -----------------------------
require_file(CFG$out$emapper_tidy,
             "Corre primero 01_parse_emapper.R")
MAPPER_DB <- readr::read_rds(CFG$out$emapper_tidy)

for (onto in CFG$ontologies) {
  log_msg("  gene2GO LOC / ", onto)
  g2go_loc <- build_gene2GO_loc(MAPPER_DB, tx2g, ontology = onto)
  log_msg("    ", length(g2go_loc), " genes con >=1 GO ", onto)
  backup_write(g2go_loc, CFG$out$gene2GO_loc[[onto]])
}

log_msg("02_stringtie_to_refseq: OK")
