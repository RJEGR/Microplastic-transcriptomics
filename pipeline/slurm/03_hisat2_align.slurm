#!/bin/bash
###############################################################################
# 03_hisat2_align.slurm  -- Alineamiento guiado (ARRAY JOB)
# ----------------------------------------------------------------------------
# Entrada: FASTQ limpios de fastp (paso 01).
# Salida : BAM ordenado + indice + summary de alineamiento.
# Dependencia: afterok:<jobid_hisat2_build> AND afterok:<jobid_fastp_array>
###############################################################################
#SBATCH --job-name=hisat2_aln
#SBATCH --partition=cicese
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=1-12%4
#SBATCH --output=logs/hisat2_aln_%A_%a.out
#SBATCH --error=logs/hisat2_aln_%A_%a.err

set -euo pipefail
#module load hisat2/2.2.1 samtools/1.18 || true


EXPORT=/LUSTRE/apps/bioinformatica/hisat2/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3
export PATH=$PATH:$EXPORT

PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"

WORDIR="${PROJECT}/workdir"
SAMPLES="${WORDIR}/samples.tsv"

CLEANDIR="${WORDIR}/01_fastp"
IDX="${PROJECT}/Reference_dir/hisat2_index/xbMagGiga1.1"
OUTDIR="${WORDIR}/03_hisat2"
mkdir -p "${OUTDIR}"

IFS=$'\t' read -r GRUPO SAMPLE R1 R2 < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES}")
echo "[$(date +%F_%T)] Aln HISAT2 | ${GRUPO} | ${SAMPLE}"

CLEAN_R1="${CLEANDIR}/${SAMPLE}_R1.clean.fq.gz"
CLEAN_R2="${CLEANDIR}/${SAMPLE}_R2.clean.fq.gz"
BAM_OUT="${OUTDIR}/${SAMPLE}.sorted.bam"
SUMMARY="${OUTDIR}/${SAMPLE}.align.summary"

[[ -s "${BAM_OUT}" ]] && { echo "  [skip] BAM existente ${BAM_OUT}"; exit 0; }

# --- Alineamiento RNA-Seq stranded? (revisar libreria; default: unstranded) -
# Si la libreria es dUTP (ScriptSeq, NEBNext Ultra II Directional): --rna-strandness RF
# Si Illumina TruSeq Stranded mRNA: --rna-strandness RF tambien.
# Por defecto, retro-compatible con el reporte 2023: SIN strand.

hisat2 \
    -x "${IDX}" \
    -1 "${CLEAN_R1}" -2 "${CLEAN_R2}" \
    -p "${SLURM_CPUS_PER_TASK}" \
    --dta \
    --no-unal \
    --new-summary --summary-file "${SUMMARY}" \
    | samtools sort -@ 4 -m 4G -o "${BAM_OUT}" -

samtools index -@ 4 "${BAM_OUT}"
samtools flagstat "${BAM_OUT}" > "${OUTDIR}/${SAMPLE}.flagstat"

echo "[$(date +%F_%T)] OK ${SAMPLE}"
