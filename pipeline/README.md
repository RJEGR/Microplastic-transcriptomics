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
# Una sola línea — el master_launcher hace parsing + sbatch + dependencias
bash pipeline/scripts/master_launcher.sh /ruta/a/Manifest.tsv
```

## Genoma de referencia

`GCF_963853765.1_xbMagGiga1.1` (assembly cromosomal, 10 cromosomas) — generación
posterior a la `GCA902806645v1` usada en 2023. Mejora la sensibilidad de
HISAT2 (más contigs anclados) y la calidad de la anotación GTF.

Ruta esperada en LUSTRE:
`/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics/Reference_dir/`
