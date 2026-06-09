#!/bin/bash
###############################################################################
# 04_stringtie.slurm  -- Reconstruccion de transcritos por muestra (ARRAY)
# ----------------------------------------------------------------------------
# Esta es la fase 1 del protocolo Pertea 2016:
#   1) StringTie por muestra (este script)
#   2) StringTie --merge (script 04b)
#   3) StringTie -eB por muestra para conteos (script 04c) -- opcional segun
#      si se usa DESeq2/tximport o Ballgown downstream.
###############################################################################
#SBATCH --job-name=stringtie
#SBATCH --partition=cicese
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --array=1-12%6
#SBATCH --output=logs/stringtie_%A_%a.out
#SBATCH --error=logs/stringtie_%A_%a.err

set -euo pipefail
module load stringtie/2.2.1 || true

PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"
SAMPLES="${PROJECT}/pipeline/scripts/samples.tsv"
BAMDIR="${PROJECT}/03_hisat2"
GTF_REF="${PROJECT}/Reference_dir/hisat2_index/xbMagGiga1.1.gtf"
OUTDIR="${PROJECT}/04_stringtie/per_sample"
mkdir -p "${OUTDIR}"

IFS=$'\t' read -r GRUPO SAMPLE R1 R2 < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES}")
BAM="${BAMDIR}/${SAMPLE}.sorted.bam"
GTF_OUT="${OUTDIR}/${SAMPLE}.gtf"

[[ -s "${GTF_OUT}" ]] && { echo "[skip] ${GTF_OUT}"; exit 0; }

stringtie "${BAM}" \
    -G "${GTF_REF}" \
    -o "${GTF_OUT}" \
    -p "${SLURM_CPUS_PER_TASK}" \
    -A "${OUTDIR}/${SAMPLE}.gene_abund.tab"

echo "[OK] ${SAMPLE} -> ${GTF_OUT}"
