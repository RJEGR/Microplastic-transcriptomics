# =============================================================================
# 05_visualize_go.R
# Plots por nivel (gene / transcript):
#   - go_bar_<level>_<ontology>_<direction>.png
#   - go_tornado_<level>_<ontology>.png (up vs down)
#   - go_parentTerm_scatter_<level>_<ontology>.png (MDS de rrvgo)
# Y tabla plana:
#   - results/topGO_summary_<level>.tsv
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

here <- normalizePath(dirname(sys.frame(1)$ofile %||% "R/05_visualize_go.R"),
                      mustWork = FALSE)
if (!file.exists(file.path(here, "config.R"))) here <- "R"
source(file.path(here, "config.R"))
source(file.path(here, "functions.R"))

log_msg("05_visualize_go: iniciando")

plot_level <- function(level = c("gene","transcript")) {
  level <- match.arg(level)
  log_msg("  === ", level, " ===")

  require_file(CFG$out$topgo_results[[level]], "Corre 04_go_enrichment.R")
  TOPGO <- readr::read_rds(CFG$out$topgo_results[[level]])

  rrvgo_ok <- FALSE
  if (file.exists(CFG$out$rrvgo_results[[level]])) {
    RRVGO_raw <- readr::read_rds(CFG$out$rrvgo_results[[level]])
    if (nrow(RRVGO_raw) > 0 &&
        all(c("go","ontology") %in% colnames(RRVGO_raw))) {
      RRVGO <- RRVGO_raw %>%
        dplyr::select(dplyr::any_of(c("go","parent","parentTerm","cluster","ontology")))
      TOPGO <- TOPGO %>%
        dplyr::left_join(RRVGO, by = c("GO.ID" = "go", "ontology" = "ontology"))
      rrvgo_ok <- TRUE
    } else {
      log_msg("  rrvgo rds sin columnas 'go'/'ontology' -> saltando join semantico")
    }
  }
  if (!rrvgo_ok) {
    TOPGO <- TOPGO %>%
      dplyr::mutate(parentTerm = NA_character_, cluster = NA_integer_)
  }

  TOPGO <- TOPGO %>%
    mutate(across(any_of(c("classicFisher","classicKS","elimKS","p.adj.ks")),
                  ~ suppressWarnings(as.numeric(gsub("<", "", .x)))))

  readr::write_tsv(TOPGO,
    file = file.path(CFG$results_dir, sprintf("topGO_summary_%s.tsv", level)))

  # ---- barra por contraste ------------------------------------------------
  for (onto in unique(TOPGO$ontology)) {
    for (dir_ in unique(TOPGO$direction)) {
      sub <- TOPGO %>%
        filter(ontology == onto, direction == dir_, !is.na(p.adj.ks)) %>%
        group_by(contrast) %>%
        slice_min(p.adj.ks, n = 15, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(Term = fct_reorder(Term, -log10(p.adj.ks + 1e-10)))
      if (nrow(sub) == 0) next
      p <- ggplot(sub, aes(y = Term, x = -log10(p.adj.ks + 1e-10),
                           fill = contrast)) +
        geom_col(width = 0.65) +
        facet_wrap(~ contrast, scales = "free_y", ncol = 1) +
        labs(x = expression(-log[10]~padj~(KS)),
             y = paste(onto, "-", dir_),
             title = sprintf("GO enrichment (%s-level) - %s / %s",
                             level, onto, dir_)) +
        theme_bw(base_size = 11) +
        theme(legend.position = "none",
              strip.background = element_rect(fill="grey92", color="white"),
              strip.text = element_text(hjust = 0),
              panel.grid.major.y = element_blank(),
              panel.grid.minor = element_blank())
      ggsave(p, filename = sprintf("go_bar_%s_%s_%s.png", level, onto, dir_),
             path = CFG$figures_dir, width = 8, height = 10, dpi = 300)
    }
  }

  # ---- tornado up vs down -------------------------------------------------
  for (onto in unique(TOPGO$ontology)) {
    sub <- TOPGO %>%
      filter(ontology == onto, !is.na(p.adj.ks)) %>%
      group_by(contrast, direction) %>%
      slice_min(p.adj.ks, n = 10, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(signed = ifelse(direction == "up",
                             -log10(p.adj.ks + 1e-10),
                              log10(p.adj.ks + 1e-10)))
    if (nrow(sub) == 0) next
    sub <- sub %>% group_by(contrast) %>%
      mutate(Term = fct_reorder(Term, abs(signed))) %>% ungroup()

    p <- ggplot(sub, aes(y = Term, x = signed, fill = direction)) +
      geom_col() + geom_vline(xintercept = 0, color = "grey50") +
      facet_wrap(~ contrast, scales = "free_y") +
      scale_fill_manual(values = c(up = "#d62728", down = "#1f77b4")) +
      labs(x = expression("<- down"~~-log[10]~padj~~"up ->"), y = onto,
           title = sprintf("GO tornado (%s-level) - ontology %s",
                           level, onto)) +
      theme_bw(base_size = 10) +
      theme(legend.position = "top",
            strip.background = element_rect(fill="grey92", color="white"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank())
    ggsave(p, filename = sprintf("go_tornado_%s_%s.png", level, onto),
           path = CFG$figures_dir, width = 12, height = 8, dpi = 300)
  }

  # ---- rrvgo MDS ----------------------------------------------------------
  if (rrvgo_ok) {
    RRVGO_FULL <- readr::read_rds(CFG$out$rrvgo_results[[level]])
    for (onto in unique(RRVGO_FULL$ontology)) {
      sub <- RRVGO_FULL %>% filter(ontology == onto)
      if (nrow(sub) < 3 || !all(c("V1","V2") %in% colnames(sub))) next
      top_labels <- sub %>% group_by(parentTerm) %>%
        slice_max(size, n = 1, with_ties = FALSE) %>% ungroup()
      p <- ggplot(sub, aes(V1, V2, color = parentTerm, size = size)) +
        geom_point(alpha = 0.7) +
        ggrepel::geom_text_repel(data = top_labels,
          aes(label = parentTerm), size = 3, show.legend = FALSE,
          max.overlaps = 15) +
        labs(title = paste("rrvgo -", level, "-", onto),
             x = "MDS 1", y = "MDS 2") +
        theme_bw(base_size = 10) + theme(legend.position = "none")
      ggsave(p,
        filename = sprintf("go_parentTerm_scatter_%s_%s.png", level, onto),
        path = CFG$figures_dir, width = 8, height = 6, dpi = 300)
    }
  }
}

for (lvl in CFG$levels) plot_level(lvl)

log_msg("05_visualize_go: OK")
