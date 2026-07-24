---
name: project-microplastics-chats
description: Índice cronológico de sesiones/chats relacionados con Microplastic-transcriptomics — qué se resolvió en cada uno y qué quedó pendiente.
metadata:
  type: project
---

# Sesiones del proyecto Microplastic-transcriptomics

Orden inverso cronológico (más reciente arriba). IDs de sesión utilizables con `mcp__session_info__read_transcript`.

## 1. `local_10270efa` — DESeq2 differential expression
**Resultado:** Análisis DE completo con PyDESeq2. Matriz 40,281 × 12 generada por `prepDE_counts.py` (workaround al StringTie sin `-e`).
**Entregables:** `deseq2_pipeline.{R,py}`, `run_de_only.py`, `prepDE_counts.py`, `ANALYSIS_NOTES.md`, 4 volcanos, heatmap top DEGs, PCA VST, 3 `ranked_*.rnk` para GSEA.
**Pendiente:** re-correr StringTie con `-e -B -G ref.gtf` para conteos canónicos.

## 2. `local_2587a71d` — Transcript biotype extraction and cohort analysis
**Resultado:** `workdir/06_biotype_annotation/results_summary.md` con overlap intra-range vía `pyranges`, asignación por max-overlap sobre catálogo de 60,992 transcritos, tablas por muestra y por cohorte (counts + %). Slides-ready para pptx.

## 3. `local_5ebe483a` — Microplastic transcriptomics count matrix
**Resultado:** `results_summary_slides.md` (12 slides) + `PCA_PC1_PC2.png`. Portada, resumen ejecutivo, diseño experimental, pipeline, construcción de matriz, varianza explicada, hallazgos, acciones recomendadas. Estructurado para la skill `pptx`.

## 4. `local_3165c792` — Emapper scripts and CDS extraction
**Resultado:** `10_emapper.slurm` con `mkdir -p "${TMPDIR_EMAPPER}"` antes del submit (emapper valida `--temp_dir` como `existing_dir`). Compatible con emapper 2.1.12 / DB 5.0.2.
**Nota runtime:** verificar `eggnog.db`, `eggnog_proteins.dmnd`, `eggnog.taxa.db` en `${EGGNOG_DB}` antes de lanzar `-m diamond`.

## 5. `local_0d3852dc` — Microplastic transcriptomics pipeline
**Resultado:** Arquitectura inicial del pipeline (11 scripts SLURM), parseo de Manifest en `samples.tsv` + `groups.tsv`, mapa CPU/RAM por etapa, propuesta científica priorizada. Base de todo el proyecto.

## Sesiones relacionadas (contexto adyacente, no del mismo repo)

- `local_83f31324` "SLURM script configuration" — RESCRIPt/QIIME2 (proyecto **LANCO**, no microplásticos). Rotó `NCBI_API_KEY` comprometida.
- `local_977d6a15` "RUN29 fastq workflow" — proyecto LANCO, deck de 11 slides sobre run multiplexado.

Cuando el usuario mencione "el chat de DE / de matriz / del emapper", mapear al ID correspondiente arriba y usar `mcp__session_info__read_transcript` antes de asumir el estado.

Ver también [[project-microplastics-pipeline]] para el estado técnico consolidado.
