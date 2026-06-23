#!/bin/bash
###############################################################################
# 04c_stringtie_quant.slurm -- Re-cuantifica con el GTF merged usando -eB
# Genera tablas de conteo listas para prepDE.py (DESeq2/edgeR via tximport).
# Dependencia: afterok:<jobid_stringtie_merge>
###############################################################################
#SBATCH --job-name=stie_quant
#SBATCH --partition=cicese
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --array=1-12%6
#SBATCH --output=logs/stie_quant_%A_%a.out
#SBATCH --error=logs/stie_quant_%A_%a.err

set -euo pipefail

EXPORT=/LUSTRE/apps/bioinformatica/stringtie/
export PATH=$PATH:$EXPORT

PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"

WORDIR="${PROJECT}/workdir"

SAMPLES="${WORDIR}/samples.tsv"

BAMDIR="${WORDIR}/03_hisat2"
MERGED_GTF="${WORDIR}/04_stringtie/merged/merged.gtf"
OUTDIR="${WORDIR}/04_stringtie/ballgown"
mkdir -p "${OUTDIR}"

IFS=$'\t' read -r GRUPO SAMPLE R1 R2 < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES}")
BAM="${BAMDIR}/${SAMPLE}.sorted.bam"
SAMP_DIR="${OUTDIR}/${SAMPLE}"
GTF_OUT="${SAMP_DIR}/${SAMPLE}.gtf"
mkdir -p "${SAMP_DIR}"

if [[ -s "${GTF_OUT}" ]]; then
    echo "[skip] ${GTF_OUT} ya existe"
else
    stringtie -e -B \
        -G "${MERGED_GTF}" \
        -o "${GTF_OUT}" \
        -A "${SAMP_DIR}/${SAMPLE}.gene_abund.tab" \
        -p "${SLURM_CPUS_PER_TASK}" \
        "${BAM}"

    echo "[OK] ${SAMPLE} cuantificado (Ballgown / prepDE compatible)"
fi
