# Pipeline R — Parseo emapper + DE + Enriquecimiento GO

**Proyecto:** Microplastic-transcriptomics (*Magallana gigas*, RNA-Seq).
**Entrada base:** `workdir/10_emapper/eggnog_xbMagGiga1.1.emapper.annotations`.
**Ejecucion:** R interactivo o `Rscript` local (no requiere SLURM).

## Resumen

Cinco scripts modulares, un orquestador y un `config.R` unico donde viven **todas** las rutas y umbrales. Se reutilizan las funciones que ya tenias (`GOenrichment`, `runtopGO`, `SEMANTIC_SEARCH`, `paste_go`) adaptadas a la salida de eggNOG-mapper v2.1 (GOs separados por coma, sin sufijo `.pN`).

## Estructura

```
Microplastic-transcriptomics/
|-- R/
|   |-- config.R                    # rutas, umbrales, contrastes
|   |-- functions.R                 # utilidades reusables
|   |-- 01_parse_emapper.R          # emapper -> gene2GO_{BP,MF,CC}.rds
|   |-- 02_stringtie_to_refseq.R    # MSTRG.X -> XM_ desde GTF merged
|   |-- 03_deseq2_cohort.R          # modelo B ~cohort (9 libs 2022)
|   |-- 04_go_enrichment.R          # topGO x rrvgo
|   `-- 05_visualize_go.R           # plots
|-- run_all.R                       # orquestador
|-- results/                        # .rds tidy + backups por fecha
`-- figures/                        # png de todas las etapas
```

## Decisiones de diseno codificadas

- **Emapper** IDs: `#query = REF|XM_...`. Se deriva `xm_id = gsub("REF|","",query)`.
- **StringTie merged.gtf** aporta el mapa `transcript_id (XM_/XR_/NM_/NR_/MSTRG.*.*) <-> gene_id (LOC*|<name> / MSTRG.*)`. Este GTF **no** usa `reference_id`; el ID RefSeq va directamente en `transcript_id` y el gen en `ref_gene_id`.
- **Conteos regenerados desde `ballgown/*.gtf`** (paso 04c, `-e -B -G merged.gtf`). El backup del run anterior (que salio del assembly per-sample con IDs `STRG.N`) esta en `workdir/05_matrix_pca/backup_YYYYMMDD/`.
- **Metadata**: no hay `Manifest.tsv`. Se deriva la cohorte del `sample_id` (`Ch14_1..3 -> Ch14`, `EP08_1..3 -> EP08`, `Se17_1..3 -> Se17`, `SES{1..3} -> SES`).
- **Dos corridas paralelas:** gene-level (LOC*) y transcript-level (XM_*). Todos los outputs llevan sufijo `_gene` o `_transcript`.
- **Modelo B**: `~cohort` sobre 9 libs 2022 (Ch14, EP08, Se17). SES excluido.
- **Contrastes**: `EP08_vs_Ch14`, `Se17_vs_Ch14`, `Se17_vs_EP08`.
- **Umbrales DE**: `padj < 0.05 & |log2FC| > 1`.
- **Ontologias**: BP + MF + CC.
- **rrvgo** con `org.Hs.eg.db` solo como fuente de la estructura semantica GO; no implica homologia con humano.
- **Backups automaticos** a `results/backup/YYYY-MM-DD/` antes de sobreescribir cualquier `.rds`.

## Antes de correr

Verifica que existan (todas ya presentes en tu repo tras la regeneracion de counts):

- `workdir/10_emapper/eggnog_xbMagGiga1.1.emapper.annotations`
- `workdir/04_stringtie/merged/merged.gtf`
- `workdir/05_matrix_pca/gene_count_matrix.csv`  (LOC*)
- `workdir/05_matrix_pca/transcript_count_matrix.csv`  (XM_*)

## Dependencias R

CRAN: `tidyverse`, `ggrepel`, `ggforce`, `scales`.
Bioconductor: `DESeq2`, `topGO`, `GO.db`, `rrvgo`, `org.Hs.eg.db`, `rtracklayer` (opcional pero recomendado).

Instalacion:
```r
install.packages(c("tidyverse", "ggrepel", "ggforce", "scales"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","topGO","GO.db","rrvgo","org.Hs.eg.db","rtracklayer"))
```

## Uso

Todo el pipeline (en la raiz del repo):

```bash
Rscript run_all.R
```

Solo unos pasos:

```bash
Rscript run_all.R 01 02        # solo emapper + mapa MSTRG
Rscript run_all.R 03 04 05     # DE + enriquecimiento + plots
```

Logs por paso quedan en `results/logs/`.

## Salidas clave

| Archivo | Descripcion |
|---|---|
| `results/emapper_annotations_tidy.rds` | Tabla emapper limpia (por xm_id) |
| `results/gene2GO_XM_{BP,MF,CC}.rds`    | gene2GO named por transcrito RefSeq (XM_) |
| `results/gene2GO_LOC_{BP,MF,CC}.rds`   | gene2GO named por gen RefSeq (LOC), agregando isoformas |
| `results/gtf_tx2gene_map.rds`          | Mapa transcript_id <-> gene_id desde merged.gtf |
| `results/dds_modelB_{gene,transcript}.rds` | DESeqDataSet por nivel |
| `results/DE_results_modelB_{gene,transcript}.rds` | Tibble long DE por nivel |
| `results/topGO_results_modelB_{gene,transcript}.rds` | Enriquecimiento long por nivel |
| `results/topGO_summary_{gene,transcript}.tsv` | Version plana para Excel |
| `results/rrvgo_semantic_modelB_{gene,transcript}.rds` | Clusters semanticos rrvgo |
| `figures/*.png` | COG/NOG, PCA, MA, barras, tornado, scatter MDS |

## Notas y limitaciones

1. **Los conteos originales estaban mal alineados**: `prepDE_counts.py` habia corrido sobre `04_stringtie/per_sample/*.gtf` (IDs `STRG.N` pre-merge). Se regeneraron desde `04_stringtie/ballgown/*.gtf` (paso 04c, con `-e -B -G merged.gtf`) para obtener IDs LOC/XM_ que coinciden con la anotacion emapper. Backup del run anterior en `workdir/05_matrix_pca/backup_YYYYMMDD/`.
2. **Transcritos novel MSTRG.*.*** no tienen anotacion GO (no estan en el proteoma expresado que anoto emapper). Quedan en las matrices de conteo pero se pierden en el enriquecimiento — no es un bug, es lo esperado.
3. **rrvgo con org.Hs.eg.db** solo agrupa terminos GO por semantica. Ninguna anotacion genica de humano se transfiere a *M. gigas*.
4. **Enriquecimiento up y down por separado** en cada contraste.

## Regeneracion de counts (por si necesitas repetirla)

```bash
cd Microplastic-transcriptomics
# construir lista sample_id -> ballgown gtf
> workdir/05_matrix_pca/ballgown_sample_list.tsv
for d in workdir/04_stringtie/ballgown/*/; do
  sid=$(basename "$d")
  gtf=$(find "$d" -name "*.gtf" | head -1)
  echo -e "${sid}\t${gtf}" >> workdir/05_matrix_pca/ballgown_sample_list.tsv
done
python3 workdir/05_matrix_pca/prepDE_counts.py \
   workdir/05_matrix_pca/ballgown_sample_list.tsv \
   workdir/05_matrix_pca/
```
