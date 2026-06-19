## 1. ARQUITECTURA DEL PIPELINE Y PARSING DEL MANIFEST

**Datos detectados en `Manifest.tsv`:** 12 librerías paired-end agrupadas en 4 cohortes biológicas — `Ch14` (Chapala, alta en invierno), `EP08` (Punta Banda, alta), `Se17` (Sesma, baja) del lote 2022, más `SES1-3` (lote 2024 nuevo, sin clasificar en el reporte previo).

**`00_parse_manifest.sh`** transforma el manifest crudo en dos artefactos:
- `samples.tsv` (12 filas) → consumido por los array jobs `SLURM_ARRAY_TASK_ID`-indexed (fastp, HISAT2 align, StringTie).
- `groups.tsv` (4 filas) → R1/R2 concatenados con coma por grupo, requerido por Trinity (`--left A,B,C`) y rnaSPAdes (`--pe<N>-1/2`). Esto **pooled-assembly** triplica la cobertura por transcrito y reduce 12 jobs masivos a 4.

**Distribución de carga SLURM** — granularidad por etapa. El clúster local
asignado a este pipeline tiene 48 nodos dedicados de 24 cores / 100 GB RAM
cada uno; `hisat2 align` se dimensionó para ocupar ~1 nodo casi completo por
tarea:

| Etapa | Tipo | Concurrencia | CPU/tarea | Memoria/tarea | Razón |
|---|---|---|---|---|---|
| fastp | array 1-12%3 | 3 | 8 | 24 GB | I/O ligado |
| hisat2-build | único, **síncrono** (`sbatch --wait`) | — | 8 | 32 GB | una sola vez; bloquea el launcher para evitar `DependencyNeverSatisfied` aguas abajo si falla |
| hisat2 align | array 1-12%12 | 12 | 20 | 90 GB | un nodo casi completo por tarea; hilos repartidos entre `hisat2` y `samtools sort` (corren en paralelo en el mismo pipe) |
| stringtie | array 1-12%3 | 3 | — | 16 GB | ligero |
| Trinity | array 1-4 | 4 | — | **256 GB** | pico Inchworm/Chrysalis |
| rnaSPAdes | array 1-4 | 4 | — | 128 GB | gráfica de Bruijn |

## 2. SCRIPTS DE LANZAMIENTO

Se generaron 11 archivos en `pipeline/` (visibles arriba). Puntos críticos:

- **fastp** (paso 3 que solicitaste): sustituye Trimmomatic del reporte 2023. Incluye `--low_complexity_filter --complexity_threshold 30` para mitigar el pico GC≈70% (presunto rRNA) **sin la pérdida masiva que sufrió SortMeRNA**, además de `--correction` (overlap PE) y `--detect_adapter_for_pe`. MINLEN=36 retro-compatible.
- **HISAT2** usa `--dta` (downstream transcript assembly) y `--new-summary` (compatible MultiQC). El índice se construye con `--ss/--exon` extraídos del GTF GCF_963853765.1 → ganancia de sensibilidad ≈5% sobre splice junctions canónicos.
- **StringTie** sigue el protocolo Pertea 2016 en tres pasos: per-sample → `--merge` → `-eB` (re-quantification con GTF consolidado). Esto sí permite descubrir **isoformas novel** que no estaban en el `gffcmp` del 2023.
- **Trinity**: `--SS_lib_type RF`, `--normalize_max_read_cov 50` (corta picos rRNA residuales), `--full_cleanup` libera disco intermedio. Pooled por grupo.
- **rnaSPAdes**: k-mers `25,45,65` para reads 150 pb (cromatografía multi-k), `--ss rf`. **Verifica el stranding real de tu librería antes de lanzar**: si es unstranded omite `--ss rf` y `--SS_lib_type RF`.

## 3. MONITOREO Y FLUJO DE TRABAJO

```bash
# Lanzar todo el pipeline (1 comando). master_launcher.sh ahora BLOQUEA
# durante el build síncrono del índice HISAT2 (sbatch --wait) — usar
# nohup/tmux/screen para que sobreviva si se cierra la sesión.
nohup bash pipeline/scripts/master_launcher.sh /LUSTRE/.../microplastics/Manifest.tsv \
    > pipeline/logs/launcher.log 2>&1 &

# Monitoreo
squeue -u $USER -o '%.10i %.20j %.8T %.10M %.6D %R'
sacct  -j <jobid>  --format=JobID,JobName,State,Elapsed,MaxRSS,ReqMem
tail -f pipeline/logs/fastp_*.err pipeline/logs/launcher.log

# Reanudar tras fallo de UNA muestra (ej. tarea 7 del array fastp)
sbatch --array=7 pipeline/slurm/01_fastp.slurm

# Cancelar cadena completa preservando builds
scancel --name=hisat2_aln,stringtie,trinity,rnaspades
```

**Cadena de dependencias** (`afterok` garantiza fallback limpio):

```
parse_manifest ─► fastp ─┬─► hisat2_align ─► stringtie ─► merge ─► quant ─► multiqc
                          │
                          ├─► trinity   (pooled por grupo)
                          └─► rnaspades (pooled por grupo)

hisat2_build ya NO se encadena por --dependency=afterok (si fallaba dejaba
a hisat2_align bloqueado para siempre en DependencyNeverSatisfied). Ahora
corre SÍNCRONO vía `sbatch --wait` antes de encolar hisat2_align; el
launcher aborta si el build falla.
```

## 4. ANÁLISIS CIENTÍFICO Y PROPUESTA EXPLORATORIA

**Lectura crítica del reporte 2023.** El estudio identificó 1356 DEGs con un diseño tipo "una localidad vs. promedio de las otras dos", revelando enriquecimiento consistente en **MAPK cascade, phospholipid metabolism, ECM-receptor interaction, CS/DS degradation, apoptosis vía apoptosoma y cytokine signaling**. Las firmas son biológicamente coherentes con estrés por microplásticos (disrupción de membrana → MAPK; daño tisular → ECM remodeling; señal de muerte celular vía Cytochrome c / XIAP / BIRC3 que el reporte ya destaca). Sin embargo, hay **tres limitaciones críticas** que este nuevo pipeline está diseñado para resolver:

**(1) Subutilización del genoma.** El 2023 mapeó contra `GCA902806645v1` (236 scaffolds, L50=5) y obtuvo 74-88% mapping. El nuevo `xbMagGiga1.1` (10 cromosomas anclados, anotación NCBI RefSeq) debería incrementar el mapeo a ≥90% y, más importante, recuperar isoformas en cromosomas previamente fragmentados. El paso `04b_stringtie_merge` + `gffcompare` cuantificará explícitamente cuántos transcritos son **novel** respecto a la anotación oficial — esto es la primera prioridad exploratoria.

**(2) Pipeline mono-enfoque.** Solo se usó alineamiento guiado. La sospecha de "contaminación con pico GC≈70%" del reporte podría no ser rRNA sino **transcritos microbianos del microbioma intestinal del ostión** o **transcritos quiméricos asociados a estrés**. El ensamblaje de novo (Trinity + rnaSPAdes) permitirá: (a) BLAST de transcritos no-mapeables contra `nr` → discriminar contaminación real vs. transcritos novel del hospedero; (b) detectar **lncRNAs** ausentes del GTF actual, particularmente relevantes en respuestas inmunes de moluscos.

**(3) Diseño contrastivo subóptimo.** El modelo "una vs. promedio de las otras dos" mezcla efecto-localidad con efecto-microplástico. Con la incorporación de `SES1-3` (lote 2024) recomiendo **reestructurar el modelo en DESeq2** como `~ localidad + carga_MP + localidad:carga_MP` cuando estén disponibles las mediciones cuantitativas de MP por sitio, o como mínimo `~ batch + grupo` para controlar el efecto de lote secuenciación.

**Próximos pasos prioritarios (en orden de retorno científico):**

a. **Análisis de isoformas diferenciales con IsoformSwitchAnalyzeR** sobre el GTF merged → particularmente en `Cytochrome c`, `XIAP`, `BIRC3` y receptores Toll-like (TLR6), donde el reporte 2023 sólo reportó gene-level. Cambios de isoforma sin cambio en expresión total son una firma clásica de estrés xenobiótico subletal.

b. **Anotación funcional específica para detoxificación** — pipeline EnTAP/Trinotate sobre los ensamblajes de novo, focalizando familias **CYP450, GST, sulfotransferasas y ABC transporters**. El reporte 2023 sólo reportó `GST omega-1`; existen ≥20 familias GST en *Magallana* que podrían responder diferencialmente.

c. **Co-expresión WGCNA** estructurada por carga de MP — esperaría módulos hub centrados en MAPK14/JNK y NF-κB que el ORA no detectó porque opera gen-a-gen.

d. **Validación dirigida con PCR digital** de 8-10 genes "centinela" (GADD45α, XIAP, Trehalase, BMP-binding endothelial regulator, Toll-like 6) para confirmar las direcciones de cambio antes de cualquier publicación.

e. **Búsqueda activa de small RNAs / miRNAs** si quedan lecturas en los outputs `unmapped` — la regulación post-transcripcional vía miRNAs es un mecanismo bien descrito de respuesta a xenobióticos en bivalvos.

Sources: `Reporte_Mgigas_2023.pdf` (Tripp Valdez, 2023, CIBNOR/CICESE); `Manifest.tsv` (lote 2022 + lote 2024); NCBI Assembly GCF_963853765.1 (xbMagGiga1.1).