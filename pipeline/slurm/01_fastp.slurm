#!/bin/bash
###############################################################################
# 01_fastp.slurm  -- Limpieza paired-end con fastp (ARRAY JOB)
# ----------------------------------------------------------------------------
# - Lee samples.tsv (col1=grupo, col2=sample_id, col3=R1, col4=R2)
# - Cada tarea del array procesa UNA libreria.
# - Corrige el problema reportado en 2023 (contaminacion 70% GC -> rRNA):
#     * deteccion automatica de adaptadores + correccion de overlap
#     * filtro low-complexity (mitiga rRNA repetitivo sin perdida masiva)
#     * MINLEN 36 (igual a referencia 2023 para retro-comparabilidad)
###############################################################################
#SBATCH --job-name=fastp_mp
#SBATCH --partition=cicese
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --array=1-12%6                      # 12 muestras, max 6 concurrentes
#SBATCH --output=logs/fastp_%A_%a.out
#SBATCH --error=logs/fastp_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rgomez41@uabc.edu.mx

set -euo pipefail

FASTQC=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software
export PATH=$PATH:$FASTQC

# --- Configuracion ---------------------------------------------------------

PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"

WORDIR="${PROJECT}/workdir"
SAMPLES="${WORDIR}/samples.tsv"
OUTDIR="${WORDIR}/01_fastp"
QCDIR="${OUTDIR}/qc_reports"
mkdir -p "${OUTDIR}" "${QCDIR}"

# --- Extraer fila correspondiente al SLURM_ARRAY_TASK_ID -------------------
IFS=$'\t' read -r GRUPO SAMPLE R1 R2 < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES}")
echo "[$(date +%F_%T)] TASK ${SLURM_ARRAY_TASK_ID} | ${GRUPO} | ${SAMPLE}"
echo "  R1: ${R1}"
echo "  R2: ${R2}"

OUT_R1="${OUTDIR}/${SAMPLE}_R1.clean.fq.gz"
OUT_R2="${OUTDIR}/${SAMPLE}_R2.clean.fq.gz"

# --- Idempotencia: skip si ya esta hecho -----------------------------------
if [[ -s "${OUT_R1}" && -s "${OUT_R2}" ]]; then
    echo "  [skip] outputs ya existen para ${SAMPLE}"
    exit 0
fi

# --- Ejecucion fastp -------------------------------------------------------
fastp \
    --in1  "${R1}"  --in2  "${R2}" \
    --out1 "${OUT_R1}" --out2 "${OUT_R2}" \
    --detect_adapter_for_pe \
    --correction \
    --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 15 \
    --qualified_quality_phred 20 \
    --length_required 36 \
    --low_complexity_filter --complexity_threshold 30 \
    --json "${QCDIR}/${SAMPLE}.fastp.json" \
    --html "${QCDIR}/${SAMPLE}.fastp.html" \
    --thread "${SLURM_CPUS_PER_TASK}" \
    2>&1 | tee "${QCDIR}/${SAMPLE}.fastp.log"

echo "[$(date +%F_%T)] FIN ${SAMPLE} (task ${SLURM_ARRAY_TASK_ID})"
