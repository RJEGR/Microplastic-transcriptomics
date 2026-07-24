---
name: project-microplastics-pipeline
description: Proyecto Microplastic-transcriptomics — RNA-Seq de Magallana gigas expuesto a microplásticos. Estado del pipeline, cohortes, decisiones de diseño y hallazgos hasta 2026-07-16.
metadata:
  type: project
---

# Microplastic-transcriptomics — Estado del proyecto

**Objetivo:** caracterizar la respuesta transcriptómica del ostión *Magallana gigas* a exposición a microplásticos, extendiendo el reporte 2023 (Tripp Valdez, CIBNOR/CICESE) con nuevo genoma anclado a cromosoma, ensamblaje de novo y modelo de DE re-especificado.

## Datos

- **12 librerías paired-end (150 pb)** en `Manifest.tsv`, agrupadas en 4 cohortes:
  - `Ch14` — Chapala, alta carga MP (invierno) — lote 2022
  - `EP08` — Punta Banda, alta carga MP — lote 2022
  - `Se17` — Sesma, baja carga MP — lote 2022
  - `SES` (SES1-3) — lote 2024, sin clasificación cuantitativa previa
- **Genoma de referencia:** NCBI `GCF_963853765.1` `xbMagGiga1.1` (10 cromosomas anclados, RefSeq). Reemplaza al `GCA902806645v1` fragmentado del 2023.

## Arquitectura del pipeline (11 pasos SLURM en `pipeline/slurm/`)

Mapa etapa→partición (ver [[hpc-omica-slurm]] para justificación):

| Paso | Script | Partición | CPU | RAM | Array |
|---|---|---|---:|---:|---|
| 01 | `01_fastp.slurm` | cicese | 8 | 24 G | 1-12%6 |
| 02 | `02_hisat2_build.slurm` | cicese | 8 | 32 G | — (síncrono con `sbatch --wait`) |
| 03 | `03_hisat2_align.slurm` | d30 | 20 | 90 G | 1-12%12 |
| 04 | `04_stringtie.slurm` | cicese | 8 | 16 G | 1-12%6 |
| 04b | `04b_stringtie_merge.slurm` | cicese | 8 | 24 G | — |
| 04c | `04c_stringtie_quant.slurm` | cicese | 8 | 16 G | 1-12%6 |
| 05 | `05_trinity.slurm` | d30 | 20 | 40 G | 1-4 |
| 06 | `06_rnaspades.slurm` | d30 | 20 | 40 G | 1-4 |
| 07 | `07_multiqc.slurm` | cicese | 4 | 24 G | — |
| 08 | `08_extract_cds.slurm` | cicese | 4 | 16 G | — |
| 09 | `09_transdecoder.slurm` | cicese | 8 | 24 G | — |
| 10 | `10_emapper.slurm` | cicese | 20 | 100 G | — |

**Cadena de dependencias:** `parse_manifest → fastp → hisat2_align → stringtie → merge → quant → multiqc`; `trinity` y `rnaspades` corren pooled por grupo en paralelo. `hisat2_build` se ejecuta SÍNCRONO (`sbatch --wait`) desde `master_launcher.sh` para evitar `DependencyNeverSatisfied`.

**Why:** cluster hpc_omica tiene sólo `cicese` (6d) y `d30` (30d) como particiones prácticas; el resto (d15/d45/d60/mayor) es reserva táctica.
**How to apply:** al crear/editar cualquier `.slurm` de este pipeline, respetar el mapa. Sólo migrar a `d30` si tiempo estimado > 5 h o memoria > 40 G.

## Decisiones de diseño registradas

1. **DE model (DESeq2):** `SES` (lote 2024) está perfectamente confundido con `batch` → `~batch + cohort` no es full-rank. Se corren DOS modelos:
   - (B) **Limpio** — `~cohort` sobre las 9 muestras 2022, 3 contrastes pareados.
   - (C) **Exploratorio** — `SES vs Se17`, etiquetado `BATCH_CONFOUNDED` en filename.
2. **StringTie sin `-e`:** el corrido inicial fue en modo assembly, `prepDE.py3` oficial rechaza esos GTF. Se escribió `prepDE_counts.py` con la fórmula `ceil(cov × len / 150)`. **Acción pendiente:** re-correr StringTie con `-e -B -G ref.gtf` para conteos canónicos antes del análisis funcional final.
3. **Diseño 2023 subóptimo:** modelo "una localidad vs. promedio de las otras dos" mezcla localidad × microplástico. Reestructurar a `~ localidad + carga_MP + localidad:carga_MP` cuando haya mediciones cuantitativas de MP; mientras tanto usar `~ batch + grupo`.
4. **fastp reemplaza Trimmomatic** con `--low_complexity_filter --complexity_threshold 30` para mitigar pico GC≈70% sin la pérdida masiva de SortMeRNA. MINLEN=36 retrocompatible con 2023.
5. **Ensamblaje de novo pooled por grupo** (Trinity + rnaSPAdes) — triplica cobertura por transcrito y reduce 12 jobs a 4.

## Hallazgos DE (padj<0.05, |LFC|>1)

| Contraste | DEGs | Up | Down |
|---|---:|---:|---:|
| EP08 vs Ch14 | 1,075 | 350 | 725 |
| Se17 vs Ch14 | 875 | 415 | 460 |
| Se17 vs EP08 | 1,093 | 722 | 371 |
| SES vs Se17 ⚠ | 3,509 | 3,008 | 501 |

Inflado de SES vs Se17 (~3× contrastes limpios) = firma típica de batch effect, no necesariamente biología. **EP08 y Ch14 (ambos altos en MP) comparten más estado transcriptómico entre sí** → respuesta convergente a microplásticos.

## Próximos pasos priorizados

a. IsoformSwitchAnalyzeR sobre GTF merged (Cytochrome c, XIAP, BIRC3, TLR6).
b. EnTAP/Trinotate focalizado en CYP450, GST (≥20 familias en *Magallana*), sulfotransferasas, ABC transporters.
c. WGCNA estructurado por carga de MP.
d. Validación con PCR digital de 8-10 genes centinela (GADD45α, XIAP, Trehalase, BMP-BER, TLR6).
e. Búsqueda de miRNAs en reads `unmapped`.

## Rutas clave

- Repo local: `/Users/rjegr/Documents/GitHub/Microplastic-transcriptomics/`
- Datos HPC: `/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics/`
- Software HPC: `/LUSTRE/apps/bioinformatica/` y `/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/`
- eggNOG DB: `/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Trinotate/TRINOTATE_DB/EGGNOG_DATA_DIR/`

Ver también [[project-microplastics-chats]] para el historial de sesiones y [[feedback-slurm-reporting-style]] para preferencias de formato.
