#
# RUN GO ENRICHMENT ANALYSIS 
# USING TOPGO AND REVIGO 
# RICARDO GOMEZ-REYES
# 2025

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

dir <- "~/Documents/MANO_DELEON/DATOS_Paulina_dir/"

subdir <- ""

URL <- "https://raw.githubusercontent.com/RJEGR/Small-RNASeq-data-analysis/master/FUNCTIONS.R"

source(URL)

# RES.P <- read_rds(paste0(path, "/DESEQ2RES.rds")) %>% filter( padj < 0.05 & abs(log2FoldChange) > 2) 

RES.P <- read_rds(paste0(file.path(dir, subdir), "/p05_DESEQ_multiple_contrast_condition.rds")) %>%
  drop_na(padj) %>%
  filter(abs(log2FoldChange) > 2 & padj < 0.05)


QUERIES <- RES.P %>% 
  mutate(sampleX = ifelse( sign(log2FoldChange) == 1, sampleA, sampleB)) %>%
  # filter(log2FoldChange < 0 ) %>% # ONLY UP-EXPRESSED IN EXPERIMENTAL CNTRST (i.e. DOWN-EXP. IN Ref)
  group_by(ids) %>%
  summarise(across(sampleX, .fns = list), n = n()) 

QUERIES %>% dplyr::count(n)

# Unique genes

QUERIES <- QUERIES %>% filter(n == 1) %>% unnest(sampleX)

QUERIES %>% dplyr::count(sampleX)

# ANNOT ====

orgdb <- "org.Hs.eg.db"

# semdata <- GOSemSim::godata(orgdb, ont="BP")

# write_rds(semdata, file = paste0(dir, orgdb, ".PB.rds"))

semdata <- read_rds(paste0(dir, orgdb, ".PB.rds"))

go_file <- paste0(dir, 'Trinotate_transcripts.txt')

MAP <- topGO::readMappings(go_file)

# Continue w/....

NAME2GO <- data.frame(ids = rep(names(MAP),
  sapply(MAP, length)),
  GO.ID = unlist(MAP), row.names = NULL) %>% as_tibble() %>%
  right_join(QUERIES) %>%
  group_by(sampleX, ids) %>%
  summarise(across(GO.ID, .fns = paste_go), n = n())

gene2GO <- split(strsplit(NAME2GO$GO.ID, ";") , NAME2GO$ids)

gene2GO <- lapply(gene2GO, unlist)


print(which_sam <- QUERIES %>% distinct(sampleX) %>% pull())

# TOPGO ====

allRes <- list()

for (i in which_sam) {
  
  cat("\nRunning GO enrichment for ", i)
  
  query.names <- NAME2GO %>% filter(sampleX %in% i) %>% distinct(ids) %>% pull()
  
  cat("\nUsing ", length(query.names), "transcrits")
  
  query.p <- RES.P %>% filter(ids %in% query.names) %>%
    group_by(ids) %>% sample_n(1) %>%
    pull(padj, name = ids)
  
  query.p <- query.p[match(query.names, names(query.p))]
  
  # identical(names(query.p), query.names)
  
  df <- GOenrichment(query.p, query.names, gene2GO, Nodes = 25, onto = "BP")
  
  allRes[[i]] <- data.frame(df, sampleX = i)
}

data <- do.call(rbind, allRes) %>% as_tibble()

# write_rds(data, file = paste0(path, "DESEQ2TOPGO.rds"))

# data <- read_rds(paste0(path, "DESEQ2TOPGO.rds"))

# REVIGO ====
# Additionally, run semantic simmilarity analysis in order to split parent terms themes

GO.IDS <- data %>% distinct(GO.ID) %>% pull() %>% sort()

SEMANTIC_SEARCH <- function(x, orgdb = "org.Ce.eg.db", semdata = semdata) {
  
  
  require(rrvgo)
  
  # semdata <- read_rds(paste0(wd, orgdb, ".rds"))
  
  x <- sort(x)
  
  SimMatrix <- calculateSimMatrix(x, 
    orgdb = orgdb,
    ont="BP", 
    semdata = semdata,
    method = 'Wang')
  
  data <- reduceSimMatrix(SimMatrix, threshold = 0.9, orgdb = orgdb) 
  
  y <- cmdscale(as.matrix(as.dist(1 - SimMatrix)), eig = TRUE, k = 2)
  
  data <- cbind(as.data.frame(y$points), data[match(rownames(y$points), data$go),])
  
  return(data)
}

REVIGO <- SEMANTIC_SEARCH(GO.IDS, orgdb, semdata)

data <- left_join(data, REVIGO, by = c("GO.ID" = "go"))

# write_rds(data, file = paste0(dir, "p05_DESEQ_multiple_contrast_condition_deseq2topgo25_multiple_contrast_condition.rds"))

data <- read_rds(paste0(dir, "p05_DESEQ_multiple_contrast_condition_deseq2topgo25_multiple_contrast_condition.rds"))

data %>% 
  distinct(sampleX, parentTerm)

# PLOT

recode_to <- c(`Bas` = "A) Bas (560 up-genes)",
  `Cao` = "B) Cao (161 up-genes)",
  `Cte` = "C) Cte (242 up-genes)")


data %>%
  dplyr::mutate(sampleX = dplyr::recode_factor(sampleX, !!!recode_to)) %>%
  group_by(sampleX) %>%
  mutate(size = size / max(size)) %>%
  filter(size > 0) %>%
  # arrange(desc(size), .by_group = T) %>%
  # mutate(Term = factor(Term)) %>%
  mutate(Term = fct_reorder(Term, size, .desc = F)) %>%
  # mutate(term = fct_reorder2(Term, sampleB, size, .desc = F)) %>%
  ggplot(aes(y = Term, x = size, color = p.adj.ks)) + 
  ggforce::facet_col(sampleX ~ ., scales = "free_y") +
  geom_segment(aes(xend = 0, yend = Term), linewidth = 2) +
  labs(y = "Biological process", x = "Enrichment frac.") +
  scale_color_viridis_c("padj", option = "viridis", direction = 1) +
  theme_bw(base_family = "GillSans", base_size = 12)  +
  theme(
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    strip.text = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) -> p

# p
ggsave(p, filename = 'p05_DESEQ_multiple_contrast_condition_deseq2topgo25_multiple_contrast_condition.png', path = dir, width = 6, height = 10, device = png, dpi = 300)

# By design =====


RES.P <- read_rds(paste0(file.path(dir, subdir), "/p05_DESEQ_multiple_contrast.rds")) %>%
  filter(sampleA != "BLA_Reg" & !sampleB  %in% c("BLA_Reg", "LOL_Reg")) %>%
  drop_na(padj) %>%
  filter(abs(log2FoldChange) > 2 & padj < 0.05)


QUERIES <- RES.P %>% 
  mutate(sampleX = ifelse( sign(log2FoldChange) == 1, sampleA, sampleB)) %>%
  distinct(ids, sampleX) %>%
  group_by(ids) %>%
  summarise(across(sampleX, .fns = list), n = n()) 

QUERIES %>% dplyr::count(n)

QUERIES <- QUERIES %>% filter(n == 1) %>% unnest(sampleX)

QUERIES %>% dplyr::count(sampleX, sort = T)

# ANNOT ====

orgdb <- "org.Hs.eg.db"

semdata <- read_rds(paste0(dir, orgdb, ".PB.rds"))

go_file <- paste0(dir, 'Trinotate_transcripts.txt')

MAP <- topGO::readMappings(go_file)

NAME2GO <- data.frame(ids = rep(names(MAP),
  sapply(MAP, length)),
  GO.ID = unlist(MAP), row.names = NULL) %>% as_tibble() %>%
  right_join(QUERIES) %>%
  group_by(sampleX, ids) %>%
  summarise(across(GO.ID, .fns = paste_go), n = n())

gene2GO <- split(strsplit(NAME2GO$GO.ID, ";") , NAME2GO$ids)

gene2GO <- lapply(gene2GO, unlist)

print(which_sam <- QUERIES %>% distinct(sampleX) %>% pull())

allRes <- list()

for (i in which_sam) {
  
  cat("\nRunning GO enrichment for ", i)
  
  query.names <- NAME2GO %>% filter(sampleX %in% i) %>% distinct(ids) %>% pull()
  
  cat("\nUsing ", length(query.names), "transcrits")
  
  query.p <- RES.P %>% filter(ids %in% query.names) %>%
    group_by(ids) %>% sample_n(1) %>%
    pull(padj, name = ids)
  
  query.p <- query.p[match(query.names, names(query.p))]
  
  # identical(names(query.p), query.names)
  
  df <- GOenrichment(query.p, query.names, gene2GO, Nodes = 25, onto = "BP")
  
  allRes[[i]] <- data.frame(df, sampleX = i)
}

data <- do.call(rbind, allRes) %>% as_tibble()

GO.IDS <- data %>% distinct(GO.ID) %>% pull() %>% sort()

REVIGO <- SEMANTIC_SEARCH(GO.IDS, orgdb, semdata)

data <- left_join(data, REVIGO, by = c("GO.ID" = "go"))

write_rds(data, file = paste0(dir, "p05_DESEQ_multiple_contrast_deseq2topgo25.rds"))

data <- read_rds(paste0(dir, "p05_DESEQ_multiple_contrast_deseq2topgo25.rds"))

data %>% 
  distinct(sampleX)

recode_to <- c(`BLA_Bas` = "A) BLA_Bas (228 up-genes)",
  `LOL_Bas` = "B) LOL_Bas (150 up-genes)",
  `LOL_Cte` = "C) LOL_Cte (150 up-genes)",
  `BLA_Cte` = "D) BLA_Cte (126 up-genes)",
  `BLA_Cao` = "E) BLA_Cao (80 up-genes)",
  `LOL_Cao` = "F) LOL_Cao (66 up-genes)")

data %>%
  separate(sampleX, into = c("Site", "Condition"), sep = "_", remove = F) %>%
  dplyr::mutate(sampleX = dplyr::recode_factor(sampleX, !!!recode_to)) %>%
  group_by(sampleX) %>%
  mutate(size = size / max(size)) %>%
  filter(size > 0.05) %>% # & p.adj.ks < 0.05
  mutate(Term = fct_reorder(Term, size, .desc = F)) %>%
  ggplot(aes(y = Term, x = size, color = Condition)) + 
  ggforce::facet_col(Site ~., scales = "free_y", space = "free") +
  geom_segment(aes(xend = 0, yend = Term, alpha = -log10(p.adj.ks)), linewidth = 2) +
  labs(y = "Biological process", x = "Enrichment frac.") +
  # scale_color_viridis_c("padj", option = "viridis", direction = 1) +
  theme_bw(base_family = "GillSans", base_size = 12)  +
  theme(
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    strip.text = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) -> p

ggsave(p, filename = 'p05_DESEQ_multiple_contrast_deseq2topgo25.png', path = dir, width = 7, height = 10, device = png, dpi = 300)

recode_to <- c(`Bas` = "A) Bas 228/150 up-genes)",
  `Cao` = "B) Cao 80/66 up-genes)",
  `Cte` = "C) Cte 150/126 up-genes)")


DataViz <- data %>%
  separate(sampleX, into = c("Site", "Condition"), sep = "_", remove = F) %>%
  group_by(sampleX, Site, Condition) %>%
  # summarise(size = sum(size)) %>%
  mutate(frac = size / max(size)) %>%
  arrange(desc(frac), .by_group = T) %>% 
  mutate(Label = term, row_number = row_number(Label)) %>% 
  mutate(Label = factor(paste(Label, row_number, sep = "__"), 
    levels = rev(paste(Label, row_number, sep = "__")))) %>%
  ungroup() %>% filter(frac > 0) %>%
  dplyr::mutate(Condition = dplyr::recode_factor(Condition, !!!recode_to)) %>%
  mutate(Site = factor(Site, levels = c("BLA", "LOL")))

DataViz %>%
  mutate(frac = ifelse(Site == "BLA", frac*-1, frac*1)) %>% 
  ggplot(aes(y = Label, x = frac, color = Site)) + # 
  facet_wrap(~ Condition, scales = "free_y") +
  geom_segment(aes(xend = 0, yend = Label), linewidth = 2) +
  scale_y_discrete(labels = function(x) gsub("__.+$", "", x)) +
  labs(y = "Biological process", x = "Enrichment ratio") +
  theme_bw(base_family = "GillSans", base_size = 7) +
  scale_color_manual("",values = c("gray20","gray")) +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    # panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
    # axis.text.x = element_text(angle = 90, size = 7)
  ) -> p2

p2

data_text <- DataViz %>% 
  mutate(Text = gsub("__.+$", "", Label)) #%>%
# mutate(Text = stringr::str_to_sentence(Text))

p2 <- p2 + 
  # scale_y_discrete(labels = scale_y_disc) +
  geom_text(data = filter(data_text, Site %in% "BLA"), 
    aes(label=Text), x = 0.02, hjust=0, vjust = 0.1, size = 2.5, family = "GillSans", color = "gray35") + 
  geom_text(data = filter(data_text, Site %in% "LOL"), 
    aes(label=Text), x = -0.02, hjust=1, vjust = 0.1, size = 2.5, family = "GillSans", color = "gray35") +
  
  guides(color = guide_legend(title = "",
    keywidth = 1,keyheight = 1,
    label_size = 0.25, ncol = 2))

p2 <- p2 + scale_x_continuous(limits = c(-1.5, 1.5), breaks = seq(-1,1, by = 0.5))

ggsave(p2, filename = 'p05_DESEQ_multiple_contrast_deseq2topgo25_tornado.png', path = dir, width = 14, height = 5, device = png, dpi = 300)

# Rerunn Go semantic for those intersected genes ====

