# =============================================================================
# 04_go_enrichment.R
# Enriquecimiento GO con topGO por (nivel x ontologia x contraste x direccion)
# + reduccion semantica de terminos con rrvgo.
#
# Requiere:
#   01_parse_emapper.R    -> gene2GO_XM_*.rds
#   02_stringtie_to_refseq.R -> gene2GO_LOC_*.rds + gtf_tx2gene_map.rds
#   03_deseq2_cohort.R    -> DE_results_modelB_{gene,transcript}.rds
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(topGO)
})

here <- normalizePath(dirname(sys.frame(1)$ofile %||% "R/04_go_enrichment.R"),
                      mustWork = FALSE)
if (!file.exists(file.path(here, "config.R"))) here <- "R"
source(file.path(here, "config.R"))
source(file.path(here, "functions.R"))

log_msg("04_go_enrichment: iniciando")

run_topgo_level <- function(level = c("gene","transcript")) {
  level <- match.arg(level)
  log_msg("  === nivel: ", level, " ===")

  require_file(CFG$out$de_results[[level]], "Corre 03_deseq2_cohort.R")
  DE <- readr::read_rds(CFG$out$de_results[[level]])

  # gene2GO segun nivel
  gene2GO_all <- list()
  for (onto in CFG$ontologies) {
    src <- switch(level,
      gene       = CFG$out$gene2GO_loc[[onto]],
      transcript = CFG$out$gene2GO_xm[[onto]])
    require_file(src, paste0("Falta gene2GO ", level, "/", onto))
    gene2GO_all[[onto]] <- readr::read_rds(src)
    log_msg("  gene2GO ", level, "/", onto, ": ",
            length(gene2GO_all[[onto]]), " IDs")
  }

  contrasts  <- unique(DE$contrast)
  directions <- c("up", "down")
  all_res <- list(); counter <- 0

  for (ctr in contrasts) {
    for (dir_ in directions) {
      degs <- select_deg(DE %>% filter(contrast == ctr),
                         direction = dir_,
                         alpha = CFG$alpha, lfc = CFG$lfc_thresh)
      log_msg("  ", ctr, " / ", dir_, ": ", nrow(degs), " DEGs")
      if (nrow(degs) < 10) { log_msg("    (saltando: <10 DEGs)"); next }

      query.p     <- degs$padj
      query.names <- degs$ids

      for (onto in CFG$ontologies) {
        counter <- counter + 1
        log_msg("    -> topGO ", onto)
        df <- tryCatch(
          GOenrichment(query.p, query.names,
                       gene2GO = gene2GO_all[[onto]],
                       cons    = CFG$topgo_conservative,
                       onto    = onto, Nodes = CFG$topgo_nodes),
          error = function(e) {
            warning("topGO fallo (", ctr, "/", dir_, "/", onto, "): ",
                    conditionMessage(e)); NULL })
        if (!is.null(df) && nrow(df) > 0) {
          all_res[[counter]] <- df %>%
            mutate(contrast = ctr, direction = dir_, level = level)
        }
      }
    }
  }

  TOPGO <- dplyr::bind_rows(all_res)
  log_msg("  TOPGO ", level, ": ", nrow(TOPGO), " filas")
  backup_write(TOPGO, CFG$out$topgo_results[[level]])

  # rrvgo por ontologia
  RRVGO_ALL <- list()
  for (onto in CFG$ontologies) {
    go_ids <- TOPGO %>% filter(ontology == onto) %>%
      distinct(GO.ID) %>% pull(GO.ID)
    if (length(go_ids) < 3) {
      log_msg("  rrvgo ", onto, ": <3 GO IDs, salto"); next }
    log_msg("  rrvgo ", onto, ": ", length(go_ids), " GO IDs")
    rvv <- tryCatch(
      SEMANTIC_SEARCH(go_ids, orgdb = CFG$rrvgo_orgdb,
                      ontology = onto, threshold = CFG$rrvgo_threshold),
      error = function(e) {
        warning("rrvgo ", onto, ": ", conditionMessage(e)); NULL })
    if (!is.null(rvv)) RRVGO_ALL[[onto]] <- rvv
  }
  RRVGO <- dplyr::bind_rows(RRVGO_ALL)
  if (nrow(RRVGO) > 0) {
    RRVGO <- RRVGO %>% dplyr::mutate(level = level)
    backup_write(RRVGO, CFG$out$rrvgo_results[[level]])
    log_msg("  rrvgo guardado: ", nrow(RRVGO), " filas")
  } else {
    # No escribir un rds vacio (rompe el left_join en 05).
    if (file.exists(CFG$out$rrvgo_results[[level]])) {
      file.remove(CFG$out$rrvgo_results[[level]])
    }
    log_msg("  rrvgo vacio para ", level, " -> NO se guarda rds.")
  }
}

for (lvl in CFG$levels) run_topgo_level(lvl)

log_msg("04_go_enrichment: OK")
