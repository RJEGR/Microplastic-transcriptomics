# =============================================================================
# run_all.R
# Orquestador secuencial del pipeline GO/DE post-eggnog-mapper.
#
# Uso:
#   Rscript run_all.R                # corre todo desde 01
#   Rscript run_all.R 03 04 05       # corre solo esos pasos
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

steps_all <- c(
  "01" = "R/01_parse_emapper.R",
  "02" = "R/02_stringtie_to_refseq.R",
  "03" = "R/03_deseq2_cohort.R",
  "04" = "R/04_go_enrichment.R",
  "05" = "R/05_visualize_go.R"
)

steps <- if (length(args) == 0) steps_all else steps_all[args]

if (any(is.na(names(steps))) || any(is.na(steps))) {
  stop("Pasos invalidos. Usa dos digitos: 01 02 03 04 05")
}

dir.create("results/logs", recursive = TRUE, showWarnings = FALSE)

for (k in seq_along(steps)) {
  step_id <- names(steps)[k]
  script  <- steps[k]

  log_file <- file.path("results", "logs",
    sprintf("%s_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"), step_id))

  cat(sprintf("\n=== Paso %s: %s ===\n", step_id, script))
  t0 <- Sys.time()

  con <- file(log_file, open = "wt")
  sink(con, split = TRUE); sink(con, type = "message")
  err <- tryCatch(source(script, echo = FALSE, chdir = FALSE),
                  error = function(e) e)
  sink(type = "message"); sink()
  close(con)

  if (inherits(err, "error")) {
    cat(sprintf("[FAIL] paso %s: %s\n", step_id, conditionMessage(err)))
    cat("  Log completo: ", log_file, "\n")
    quit(status = 1)
  }

  cat(sprintf("[OK]   paso %s (%.1fs) log: %s\n",
              step_id, as.numeric(Sys.time() - t0, units = "secs"),
              log_file))
}

cat("\n=== Pipeline terminado ===\n")
