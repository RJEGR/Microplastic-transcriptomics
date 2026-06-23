#!/bin/bash
###############################################################################
# 03_hisat2_align.slurm  -- Alineamiento guiado (ARRAY JOB)
# ----------------------------------------------------------------------------
# Entrada: FASTQ limpios de fastp (paso 01).
# Salida : BAM ordenado + indice + summary de alineamiento.
# Dependencia: afterok:<jobid_fastp_array>  (el indice HISAT2 ya se construye
#              de forma sincrona en master_launcher.sh antes de llegar aqui)
# Nodos del cluster local: 24 cores / 100GB RAM cada uno. Cada tarea del
# array ocupa un nodo casi completo (20 cores, 90GB) para evitar OOM-kill
# por el sort de samtools, que reserva memoria POR HILO.
###############################################################################
#SBATCH --job-name=hisat2_aln
#SBATCH --partition=d30
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=90G
#SBATCH --time=12:00:00
#SBATCH --array=1-12%12
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
BAM_BAI="${OUTDIR}/${SAMPLE}.sorted.bam.bai"
SUMMARY="${OUTDIR}/${SAMPLE}.align.summary"
FLAGSTAT="${OUTDIR}/${SAMPLE}.flagstat"

# --- Alineamiento RNA-Seq stranded? (revisar libreria; default: unstranded) -
# Si la libreria es dUTP (ScriptSeq, NEBNext Ultra II Directional): --rna-strandness RF
# Si Illumina TruSeq Stranded mRNA: --rna-strandness RF tambien.
# Por defecto, retro-compatible con el reporte 2023: SIN strand.

if [[ -s "${BAM_OUT}" && -s "${SUMMARY}" ]]; then
    echo "  [SKIP] BAM y summary ya existen para ${SAMPLE}, omitiendo hisat2"
else
    # hisat2 y samtools sort corren en paralelo dentro del mismo pipe, asi
    # que los hilos deben repartirse entre ambos (no asignar CPUS_PER_TASK
    # completo a cada uno). samtools sort -m es memoria POR HILO: con
    # SORT_THREADS=5 y -m 3G el sort usa <=15G, dejando margen dentro de
    # los 90G/20 cores reservados arriba para hisat2 + overhead del SO.
    SORT_THREADS=$(( SLURM_CPUS_PER_TASK / 4 ))
    ALN_THREADS=$(( SLURM_CPUS_PER_TASK - SORT_THREADS ))
    hisat2 \
        -x "${IDX}" \
        -1 "${CLEAN_R1}" -2 "${CLEAN_R2}" \
        -p "${ALN_THREADS}" \
        --dta \
        --no-unal \
        --new-summary --summary-file "${SUMMARY}" \
        | samtools sort -@ "${SORT_THREADS}" -m 3G -T "${OUTDIR}/${SAMPLE}.tmp" -o "${BAM_OUT}" -
fi

if [[ -s "${BAM_BAI}" ]]; then
    echo "  [SKIP] Indice BAM ya existe para ${SAMPLE}, omitiendo samtools index"
else
    samtools index "${BAM_OUT}"
fi

if [[ -s "${FLAGSTAT}" ]]; then
    echo "  [SKIP] flagstat ya existe para ${SAMPLE}, omitiendo samtools flagstat"
else
    samtools flagstat "${BAM_OUT}" > "${FLAGSTAT}"
fi

echo "[$(date +%F_%T)] OK ${SAMPLE}"
