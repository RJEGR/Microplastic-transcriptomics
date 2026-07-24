---
name: hpc-omica-slurm
description: Configuración de particiones SLURM del clúster hpc_omica (CICESE) — límites de tiempo, nodos disponibles y política de asignación por etapa del pipeline.
metadata:
  type: reference
---

# Clúster `hpc_omica` (CICESE) — Particiones SLURM

Salida real de `sinfo` (2026-07-16). Todas las particiones comparten el mismo pool de 30 nodos físicos `nodo[1-30]`, salvo `large` y `docker` que están restringidas.

| PARTITION | TIMELIMIT   | NODES(A/I/O/T) | NODELIST     | Notas                                          |
|-----------|-------------|----------------|--------------|------------------------------------------------|
| cicese*   | 6-00:00:00  | 5/19/6/30      | nodo[1-30]   | Default. Elegir para trabajos < 6 días.        |
| d15       | 15-00:00:00 | 5/19/6/30      | nodo[1-30]   | Trabajos ≤ 15 días.                            |
| d30       | 30-00:00:00 | 5/19/6/30      | nodo[1-30]   | Ensamblajes de novo, DE largos.                |
| d45       | 45-00:00:00 | 5/19/6/30      | nodo[1-30]   | Sólo si realmente lo requiere el análisis.     |
| d60       | 60-00:00:00 | 5/19/6/30      | nodo[1-30]   | Idem. Justificar en el header del script.      |
| mayor     | 30-00:00:00 | 5/19/6/30      | nodo[1-30]   | Alias de facto para 30 días.                   |
| large     | 30-00:00:00 | 1/0/0/1        | nodo30       | 1 nodo dedicado — trabajos monolíticos RAM alta|
| docker    | 30-00:00:00 | 0/3/0/3        | nodo[26-28]  | Requiere contenedores. Reservar sólo si se usa docker. |

## Regla de asignación por etapa (aplicar en `pipeline/slurm/`)

- Etapas ligeras (< 6 h) y arrays con concurrencia moderada → **cicese**: `fastp`, `hisat2-build`, `stringtie`, `stringtie_merge`, `stringtie_quant`, `multiqc`, `extract_cds`, `transdecoder`, `emapper`.
- Etapas largas o de alto pico de RAM (> 12 h, > 40 GB) → **d30**: `hisat2_align`, `Trinity`, `rnaSPAdes`.
- Nunca usar `d45`/`d60` sin justificar tiempo en el header del `.slurm`.
- `large` sólo si un solo proceso monolítico necesita > 200 GB RAM o hilos concentrados en un único nodo — Trinity fase Butterfly califica.
- `docker` sólo si el binario corre por contenedor. Actualmente ningún paso del pipeline lo requiere.

## Estado del clúster al 2026-07-16

- 30 nodos totales (`A/I/O/T = 5/19/6/30`). 19 idle → ventana buena para lanzar arrays grandes.
- 6 nodos en drain/offline: revisar antes de asumir capacidad plena.
- **Discrepancia detectada:** el `README.md` del pipeline afirma "48 nodos de 24 cores / 100 GB RAM". El clúster real tiene 30 nodos. Verificar cores/RAM por nodo con `scontrol show node nodo1` y actualizar el README antes de re-dimensionar los `--cpus-per-task` / `--mem`.

## Header SLURM canónico para este proyecto

```bash
#SBATCH --job-name=<etapa>_mp
#SBATCH --partition=<cicese|d30>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=<n>
#SBATCH --mem=<X>G
#SBATCH --time=<HH:MM:SS>          # nunca dejar `time` sin declarar
#SBATCH --array=1-12%<K>            # sólo en etapas por muestra
#SBATCH --output=logs/<etapa>_%A_%a.out
#SBATCH --error=logs/<etapa>_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rgomez41@uabc.edu.mx
```

## Cheatsheet de comandos

```bash
sinfo -o "%.10P %.11l %.15F %.20N"                       # ver estado por partición
squeue -u $USER -o '%.10i %.20j %.8T %.10M %.6D %R'      # cola personal
sacct -j <jobid> --format=JobID,JobName,State,Elapsed,MaxRSS,ReqMem
scontrol show node nodo1                                  # cores/RAM reales por nodo
scancel --name=<job-name>                                 # cancelar por nombre
sbatch --array=<i> pipeline/slurm/0X_step.slurm           # relanzar tarea puntual
```

Ver también [[project-microplastics-pipeline]] para el mapa etapa→partición usado en producción.
