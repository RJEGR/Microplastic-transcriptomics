# =============================================================================
# 01_parse_emapper.R
# Parseo de eggnog_xbMagGiga1.1.emapper.annotations
#
# Genera:
#   - results/emapper_annotations_tidy.rds
#   - results/gene2GO_XM_{BP,MF,CC}.rds   (key = XM_... para transcript-level DE)
#   - figures/cog_nog_transcriptome.png
#
# gene2GO a nivel LOC se genera en 02_stringtie_to_refseq.R (necesita el mapa).
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

here <- normalizePath(dirname(sys.frame(1)$ofile %||% "R/01_parse_emapper.R"),
                      mustWork = FALSE)
if (!file.exists(file.path(here, "config.R"))) here <- "R"
source(file.path(here, "config.R"))
source(file.path(here, "functions.R"))

log_msg("01_parse_emapper: iniciando")

MAPPER_DB <- read_emapper(CFG$emapper_file)
log_msg("  ", nrow(MAPPER_DB), " filas / ",
        dplyr::n_distinct(MAPPER_DB$xm_id), " xm_id unicos")

eggNOG_cols <- c("xm_id","query","seed_ortholog","evalue","score",
                 "eggNOG_OGs","max_annot_lvl","COG_category",
                 "Description","Preferred_name","GOs","EC",
                 "PFAMs","BRITE","CAZy")
keep_cols <- c(eggNOG_cols[eggNOG_cols %in% colnames(MAPPER_DB)],
               grep("^KEGG", colnames(MAPPER_DB), value = TRUE))

backup_write(MAPPER_DB %>% select(all_of(keep_cols)), CFG$out$emapper_tidy)
log_msg("  guardado: ", CFG$out$emapper_tidy)

# ---- gene2GO XM por ontologia ----------------------------------------------
for (onto in CFG$ontologies) {
  log_msg("  gene2GO XM / ", onto)
  g2go <- build_gene2GO_xm(MAPPER_DB, ontology = onto)
  log_msg("    ", length(g2go), " transcritos con >=1 GO ", onto)
  backup_write(g2go, CFG$out$gene2GO_xm[[onto]])
}

# ---- COG/NOG plot (referencia transcriptomica) -----------------------------
NOG.col <- c(
  'J@Translation, ribosomal structure and biogenesis',
  'A@RNA processing and modification','K@Transcription',
  'L@Replication, recombination and repair',
  'B@Chromatin structure and dynamics',
  'D@Cell cycle control, cell division, chromosome partitioning',
  'Y@Nuclear structure','V@Defense mechanisms',
  'T@Signal transduction mechanisms',
  'M@Cell wall/membrane/envelope biogenesis',
  'N@Cell motility','Z@Cytoskeleton','W@Extracellular structures',
  'U@Intracellular trafficking, secretion, and vesicular transport',
  'O@Posttranslational modification, protein turnover, chaperones',
  'X@Mobilome: prophages, transposons',
  'C@Energy production and conversion',
  'G@Carbohydrate transport and metabolism',
  'E@Amino acid transport and metabolism',
  'F@Nucleotide transport and metabolism',
  'H@Coenzyme transport and metabolism',
  'I@Lipid transport and metabolism',
  'P@Inorganic ion transport and metabolism',
  'Q@Secondary metabolites biosynthesis, transport and catabolism',
  'R@General function prediction only','S@Function unknown'
) %>% tibble(NOG.col = .) %>%
  separate(NOG.col, sep = "@", into = c("EggNM.COG_category", "COG_name"))

plotdf <- MAPPER_DB %>%
  filter(!is.na(COG_category), COG_category != "-") %>%
  select(xm_id, COG_category) %>%
  mutate(COG_category = strsplit(COG_category, "")) %>%
  tidyr::unnest(COG_category) %>%
  right_join(NOG.col, by = c("COG_category" = "EggNM.COG_category")) %>%
  filter(!is.na(xm_id)) %>%
  mutate(COG_category = paste0(COG_category, ", ", COG_name)) %>%
  dplyr::count(COG_category, sort = TRUE) %>%
  mutate(frac = n / sum(n)) %>%
  arrange(desc(frac)) %>%
  mutate(COG_category = factor(COG_category, levels = unique(COG_category)),
         facet = "A) xbMagGiga1.1 expressed proteome")

p_cog <- plotdf %>%
  ggplot(aes(y = COG_category, x = frac)) +
  labs(y = "Nested Orthologous Gene Group (NOGs)",
       x = "Enrichment ratio (N transcripts / Total)") +
  facet_grid(~ facet, scales = "free_x", space = "free_x") +
  geom_col(fill = "black") +
  scale_x_continuous(labels = scales::percent_format()) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = 'grey95', color = 'white')) +
  geom_text(aes(label = paste0("(", n, ")")), size = 3.5, hjust = -0.1)

ggsave(p_cog, filename = "cog_nog_transcriptome.png",
       path = CFG$figures_dir, width = 8, height = 7, dpi = 300)

log_msg("01_parse_emapper: OK")
