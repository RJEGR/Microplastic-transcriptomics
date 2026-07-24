# =============================================================================
# functions.R
# Utilidades para parseo eggNOG-mapper, mapeo StringTie<->RefSeq,
# DESeq2, topGO enrichment y semantic similarity (rrvgo).
#
# Reutiliza funciones adaptadas de:
#   https://github.com/RJEGR/Small-RNASeq-data-analysis/blob/main/FUNCTIONS.R
#   https://github.com/RJEGR/Cancer_sete_T_assembly/blob/main/functions.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# Blindaje: cuando topGO/AnnotationDbi (via rrvgo/org.Hs.eg.db) se cargan en
# la misma sesion R (run_all.R hace source secuencial), AnnotationDbi::select
# enmascara dplyr::select y rompe tbl_df. Forzamos el default correcto.
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

# -----------------------------------------------------------------------------
# 1) PARSEO EMAPPER
# -----------------------------------------------------------------------------

#' Lee un archivo `.emapper.annotations` y devuelve un tibble limpio.
#' `#query` viene como "REF|XM_..." (RefSeq mRNA). Derivamos `xm_id` sin prefijo.
read_emapper <- function(path) {
  require_file(path, "Salida de eggNOG-mapper (paso 10_emapper)")
  OUT <- readr::read_tsv(path, comment = "##",
                         show_col_types = FALSE, progress = FALSE) %>%
    dplyr::rename(query = `#query`)

  OUT %>%
    mutate(
      query   = query,
      xm_id   = gsub("^REF\\|", "", query)
    ) %>%
    select(xm_id, query, everything()) %>%
    arrange(xm_id)
}

paste_col <- function(x) {
  x <- x[!is.na(x)]
  x <- unique(sort(x))
  paste(x, collapse = ";")
}
paste_go <- paste_col

# -----------------------------------------------------------------------------
# 2) CONSTRUCCION gene2GO
# -----------------------------------------------------------------------------

#' Explota columna GOs (coma-separada) y devuelve tibble long id/GO.
explode_emapper_go <- function(emapper_df, id_col = "xm_id") {
  emapper_df %>%
    filter(!is.na(GOs), GOs != "-") %>%
    select(all_of(id_col), GOs) %>%
    mutate(GOs = strsplit(GOs, ",")) %>%
    tidyr::unnest(GOs) %>%
    dplyr::rename(id = all_of(id_col)) %>%
    filter(nzchar(GOs))
}

#' Filtra GO IDs por ontologia usando GO.db. Devuelve mismo tibble long.
filter_go_by_ontology <- function(long_df, ontology) {
  stopifnot(ontology %in% c("BP", "MF", "CC"))
  if (!requireNamespace("GO.db", quietly = TRUE))
    stop("GO.db requerido. BiocManager::install('GO.db')")
  go_onto <- AnnotationDbi::select(GO.db::GO.db,
    keys = unique(long_df$GOs), columns = "ONTOLOGY", keytype = "GOID")
  long_df %>%
    inner_join(go_onto, by = c("GOs" = "GOID")) %>%
    filter(ONTOLOGY == ontology) %>%
    select(-ONTOLOGY)
}

#' Lista named id -> vector GO ids (formato topGO).
long_to_gene2GO <- function(long_df) {
  df <- long_df %>%
    group_by(id) %>%
    summarise(GOs = list(unique(GOs)), .groups = "drop")
  setNames(df$GOs, df$id)
}

#' Construye gene2GO a nivel XM_ (transcript RefSeq) por ontologia.
build_gene2GO_xm <- function(emapper_df, ontology = NULL) {
  long <- explode_emapper_go(emapper_df, id_col = "xm_id")
  if (!is.null(ontology)) long <- filter_go_by_ontology(long, ontology)
  long_to_gene2GO(long)
}

#' Construye gene2GO a nivel LOC (gen RefSeq) agregando GOs por gen.
#' Requiere el mapa transcript_id -> gene_id extraido de merged.gtf.
build_gene2GO_loc <- function(emapper_df, tx2gene, ontology = NULL) {
  long <- explode_emapper_go(emapper_df, id_col = "xm_id")
  if (!is.null(ontology)) long <- filter_go_by_ontology(long, ontology)

  # tx2gene: tibble con transcript_id (XM_), gene_id (LOC*)
  long %>%
    dplyr::rename(transcript_id = id) %>%
    inner_join(tx2gene, by = "transcript_id",
               relationship = "many-to-many") %>%
    filter(!is.na(gene_id)) %>%
    dplyr::rename(id = gene_id) %>%
    distinct(id, GOs) %>%
    long_to_gene2GO()
}

# -----------------------------------------------------------------------------
# 3) MAPEO merged.gtf
# -----------------------------------------------------------------------------

#' Extrae del GTF merged de StringTie un tibble
#'   transcript_id / gene_id / ref_gene_id / gene_name
#'
#' Este GTF (04b) contiene:
#'   gene_id "MSTRG.X" | transcript_id "XM_..." | ref_gene_id "LOC..." | gene_name "..."
#' o para novel: gene_id "MSTRG.X" | transcript_id "MSTRG.X.Y"
parse_stringtie_gtf <- function(gtf_path) {
  require_file(gtf_path, "merged.gtf de StringTie (04b_stringtie_merge.slurm)")

  extract_attr <- function(attr, key) {
    m <- stringr::str_match(attr, sprintf('%s "([^"]+)"', key))[, 2]
    m
  }

  cn <- c("seqnames","source","type","start","end",
          "score","strand","frame","attributes")

  df <- readr::read_tsv(gtf_path, comment = "#", col_names = cn,
                        show_col_types = FALSE, progress = FALSE) %>%
    filter(type == "transcript")

  df %>%
    transmute(
      transcript_id = extract_attr(attributes, "transcript_id"),
      gene_id_stringtie = extract_attr(attributes, "gene_id"),
      ref_gene_id   = extract_attr(attributes, "ref_gene_id"),
      gene_name     = extract_attr(attributes, "gene_name")
    ) %>%
    distinct()
}

#' Construye mapa "canonico" transcript_id -> gene_id que iguala la logica
#' de prepDE_counts.py:
#'   - Si hay gene_id + gene_name -> "<gene_id>|<gene_name>"
#'   - Si solo hay gene_id -> gene_id
#'   - Si el gene_id es MSTRG y hay ref_gene_id -> preferir ref_gene_id
#' Esto asegura que los IDs coincidan con gene_count_matrix.csv.
build_tx2gene_from_gtf <- function(gtf_df) {
  gtf_df %>%
    mutate(
      gene_id_used = case_when(
        # Si StringTie asigno un ref_gene_id (RefSeq LOC*), usamos ese
        !is.na(ref_gene_id) & !is.na(gene_name) ~ paste0(ref_gene_id, "|", gene_name),
        !is.na(ref_gene_id)                     ~ ref_gene_id,
        # fallback: gene_id de StringTie (MSTRG.*)
        TRUE                                    ~ gene_id_stringtie
      )
    ) %>%
    transmute(transcript_id, gene_id = gene_id_used) %>%
    distinct()
}

# -----------------------------------------------------------------------------
# 4) DESeq2 helpers
# -----------------------------------------------------------------------------

get_res <- function(dds, contrast, alpha_cutoff = 0.05) {
  stopifnot(length(contrast) == 2)
  sA <- contrast[1]; sB <- contrast[2]
  contrast_var <- as.character(DESeq2::design(dds))[2]
  meta  <- as.data.frame(SummarizedExperiment::colData(dds))
  keepA <- meta[[contrast_var]] == sA
  keepB <- meta[[contrast_var]] == sB
  res <- DESeq2::results(dds,
    contrast = c(contrast_var, sA, sB), alpha = alpha_cutoff)
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  baseMeanA <- rowMeans(norm_counts[, keepA, drop = FALSE])
  baseMeanB <- rowMeans(norm_counts[, keepB, drop = FALSE])
  out <- as.data.frame(res) %>%
    cbind(baseMeanA = baseMeanA, baseMeanB = baseMeanB, .) %>%
    cbind(sampleA = sA, sampleB = sB, .) %>%
    as_tibble(rownames = "ids") %>%
    mutate(padj = ifelse(is.na(padj), 1, padj))

  # Redondeo explicito por columna (evita quirks de across()+matches en algunas versiones)
  round_cols <- intersect(
    c("baseMean","baseMeanA","baseMeanB","log2FoldChange","lfcSE","stat"),
    colnames(out))
  for (cc in round_cols) out[[cc]] <- round(out[[cc]], 3)
  out
}

# -----------------------------------------------------------------------------
# 5) topGO wrappers
# -----------------------------------------------------------------------------

runtopGO <- function(topGOdata, topNodes = 20, conservative = TRUE) {
  suppressPackageStartupMessages(require(topGO))
  RFisher <- runTest(topGOdata, algorithm = "classic", statistic = "fisher")
  if (conservative) {
    RKS      <- runTest(topGOdata, algorithm = "classic", statistic = "ks")
    RKS.elim <- runTest(topGOdata, algorithm = "elim",    statistic = "ks")
    GenTable(topGOdata,
      classicFisher = RFisher, classicKS = RKS, elimKS = RKS.elim,
      orderBy = "elimKS", ranksOf = "classicFisher", topNodes = topNodes)
  } else {
    RKS <- runTest(topGOdata, algorithm = "classic", statistic = "ks")
    test.stat <- new("weightCount",
      testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
    weights <- getSigGroups(topGOdata, test.stat)
    GenTable(topGOdata,
      classic = RFisher, KS = RKS, weight = weights,
      orderBy = "weight", ranksOf = "classic", topNodes = topNodes)
  }
}

GOenrichment <- function(query.p, query.names, gene2GO,
                         cons = TRUE, onto = "BP", Nodes = Inf) {
  suppressPackageStartupMessages(require(topGO))
  names(query.p) <- query.names

  keep <- names(gene2GO) %in% names(query.p)
  gene2GO <- gene2GO[keep]
  keep <- names(query.p) %in% names(gene2GO)
  query.p <- query.p[keep]

  if (length(query.p) == 0) {
    warning("No hay overlap query <-> gene2GO para ", onto)
    return(NULL)
  }

  topGOdata <- new("topGOdata",
    ontology = onto, description = "topGO (emapper GOs)",
    allGenes = query.p, geneSel = function(x) x,
    annot = annFUN.gene2GO, gene2GO = gene2GO)

  allGO <- usedGO(topGOdata)
  topNodes <- if (is.infinite(Nodes)) length(allGO) else min(Nodes, length(allGO))

  allRes <- runtopGO(topGOdata, topNodes = topNodes, conservative = cons)
  p.adj.ks <- p.adjust(allRes$classicKS, method = "BH")
  allRes <- cbind(allRes, p.adj.ks)
  allRes$Term <- gsub(" [a-z]*\\.\\.\\.$", "", allRes$Term)
  allRes$Term <- gsub("\\.\\.\\.$", "", allRes$Term)
  as_tibble(allRes) %>% mutate(ontology = onto)
}

# -----------------------------------------------------------------------------
# 6) rrvgo semantic similarity
# -----------------------------------------------------------------------------

SEMANTIC_SEARCH <- function(go_ids, orgdb = "org.Hs.eg.db",
                            ontology = "BP", threshold = 0.9, semdata = NULL) {
  if (!requireNamespace("rrvgo", quietly = TRUE))
    stop("Instala rrvgo (BiocManager::install('rrvgo'))")
  if (!requireNamespace(orgdb, quietly = TRUE))
    stop("Instala ", orgdb, " (BiocManager::install('", orgdb, "'))")

  go_ids <- sort(unique(go_ids))
  if (length(go_ids) < 2) return(NULL)

  # Prefiltrar GO IDs a los que existen en la ontologia solicitada (GO.db).
  # Esto evita que calculateSimMatrix elimine silenciosamente muchos IDs y
  # devuelva una matriz vacia o con warnings.
  if (requireNamespace("GO.db", quietly = TRUE)) {
    go_ann <- suppressMessages(AnnotationDbi::select(GO.db::GO.db,
      keys = go_ids, columns = "ONTOLOGY", keytype = "GOID"))
    go_ids <- go_ann$GOID[!is.na(go_ann$ONTOLOGY) & go_ann$ONTOLOGY == ontology]
  }
  if (length(go_ids) < 2) {
    message("[SEMANTIC_SEARCH] <2 GO IDs validos para ", ontology, " tras filtrar.")
    return(NULL)
  }

  SimMatrix <- tryCatch(
    suppressWarnings(rrvgo::calculateSimMatrix(go_ids,
      orgdb = orgdb, ont = ontology, semdata = semdata, method = "Wang")),
    error = function(e) { message("[SEMANTIC_SEARCH] calculateSimMatrix: ",
                                  conditionMessage(e)); NULL })
  if (is.null(SimMatrix) || !is.matrix(SimMatrix) || nrow(SimMatrix) < 2) {
    message("[SEMANTIC_SEARCH] SimMatrix vacia para ", ontology,
            " (GO IDs ", length(go_ids), " -> filas: ",
            if (is.null(SimMatrix)) 0 else nrow(SimMatrix), ")")
    return(NULL)
  }

  reduced <- tryCatch(
    suppressWarnings(rrvgo::reduceSimMatrix(SimMatrix,
      threshold = threshold, orgdb = orgdb)),
    error = function(e) { message("[SEMANTIC_SEARCH] reduceSimMatrix: ",
                                  conditionMessage(e)); NULL })
  if (is.null(reduced) || nrow(reduced) < 1) return(NULL)

  y <- tryCatch(
    cmdscale(as.matrix(as.dist(1 - SimMatrix)), eig = TRUE, k = 2),
    error = function(e) NULL)
  if (!is.null(y)) {
    pts <- as.data.frame(y$points); colnames(pts) <- c("V1", "V2")
    reduced <- cbind(pts, reduced[match(rownames(pts), reduced$go), ])
  }
  as_tibble(reduced) %>% dplyr::mutate(ontology = ontology)
}

# -----------------------------------------------------------------------------
# 7) Utilerias generales
# -----------------------------------------------------------------------------

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
}

select_deg <- function(res, direction = c("up", "down", "all"),
                       alpha = 0.05, lfc = 1) {
  direction <- match.arg(direction)
  x <- res %>% dplyr::filter(padj < alpha)
  switch(direction,
    up   = x %>% dplyr::filter(log2FoldChange >  lfc),
    down = x %>% dplyr::filter(log2FoldChange < -lfc),
    all  = x %>% dplyr::filter(abs(log2FoldChange) > lfc))
}
