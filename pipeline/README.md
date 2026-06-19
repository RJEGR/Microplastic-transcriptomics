# Pipeline RNA-Seq Microplastics — *Magallana gigas*

Pipeline modular SLURM para reproducir y extender el análisis transcriptómico
del reporte 2023 (`Reporte_Mgigas_2023.pdf`) sobre la respuesta de *Magallana
gigas* a microplásticos en tres localidades de Baja California (Chapala / Punta
Banda / Sesma), incorporando el lote nuevo `SES1-3` (2024).

## Estructura

```
pipeline/
├── scripts/
│   ├── 00_parse_manifest.sh      # genera samples.tsv y groups.tsv
│   └── master_launcher.sh        # encadena todos los jobs con --dependency
├── slurm/
│   ├── 01_fastp.slurm            # ARRAY: limpieza paired-end
│   ├── 02_hisat2_build.slurm     # único: índice del genoma
│   ├── 03_hisat2_align.slurm     # ARRAY: alineamiento RNA-Seq
│   ├── 04_stringtie.slurm        # ARRAY: ensamblaje por muestra (guiado)
│   ├── 04b_stringtie_merge.slurm # único: GTF consolidado no-redundante
│   ├── 04c_stringtie_quant.slurm # ARRAY: re-cuantificación (-eB)
│   ├── 05_trinity.slurm          # ARRAY por grupo: de novo
│   ├── 06_rnaspades.slurm        # ARRAY por grupo: de novo
│   └── 07_multiqc.slurm          # reporte unificado
└── logs/                         # stdout/stderr de cada job
```

## Lanzamiento

```bash
PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"
#PROJECT="/Users/rjegr/Documents/GitHub/Microplastic-transcriptomics"
EXPORT="${PROJECT}/pipeline/slurm"
chmod +x "${EXPORT}/"*.slurm
export PATH=$PATH:$EXPORT
EXPORT="${PROJECT}/pipeline/scripts"
chmod +x "${EXPORT}/"*.sh
export PATH=$PATH:$EXPORT
```

```bash
# Una sola línea — el master_launcher hace parsing + sbatch + dependencias
ls -1 *gz > filenames.txt
cat filenames.txt | bash 00_parse_filenames.sh > Manifest.tsv

# IMPORTANTE: master_launcher.sh ahora bloquea durante el build del indice
# HISAT2 (sbatch --wait, ver "Cambios recientes" abajo). Lanzar bajo
# nohup/tmux/screen si la sesion puede cerrarse antes de que termine.
nohup bash master_launcher.sh Manifest.tsv > pipeline/logs/launcher.log 2>&1 &
```

## Cambios recientes (2026-06-19)

- **`02_hisat2_build.slurm`**: el filtro `zcat "${GTF_GZ}" || grep -v 'transcript_id ""'`
  usaba `||` en vez de `|`, así que `grep` nunca recibía el GTF descomprimido y
  generaba un archivo vacío/corrupto. Corregido a un pipe real; el GTF
  filtrado sigue guardándose en `${IDX_DIR}/${PREFIX}.gtf`.
- **`master_launcher.sh`**: el build del índice ya no se encadena con
  `--dependency=afterok`. Si fallaba, el job de alineamiento quedaba
  bloqueado para siempre en estado `DependencyNeverSatisfied` (la condición
  nunca podía cumplirse). Ahora se ejecuta de forma **síncrona** con
  `sbatch --wait`, y el launcher aborta con un mensaje claro si el build
  falla, antes de gastar cupo en fastp/align con un índice roto.
- **`03_hisat2_align.slurm`**: dimensionado para nodos dedicados de 24
  cores / 100 GB RAM — `--cpus-per-task=20 --mem=90G`, cupo de array
  `%12` (un nodo casi completo por tarea). Se corrigió además que
  `samtools sort -m` es memoria **por hilo**: con 4 hilos a 20G cada uno
  el job pedía hasta 80 GB solo para el sort, sobre una reserva total de
  32 GB, causando OOM-kill. Los hilos ahora se reparten entre `hisat2`
  (alineador) y `samtools sort` (que corren en paralelo dentro del mismo
  pipe) en vez de asignarse el total de `cpus-per-task` a ambos a la vez.
- Los tres pasos de `03_hisat2_align.slurm` (alineamiento, índice BAM,
  flagstat) ahora son idempotentes — cada uno se omite si su salida ya
  existe, igual que en `01_fastp.slurm`.

## Genoma de referencia

`GCF_963853765.1_xbMagGiga1.1` (assembly cromosomal, 10 cromosomas) — generación
posterior a la `GCA902806645v1` usada en 2023. Mejora la sensibilidad de
HISAT2 (más contigs anclados) y la calidad de la anotación GTF.

Ruta esperada en LUSTRE:
`/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics/Reference_dir/`
