#!/bin/bash
###############################################################################
# 05_trinity.slurm  -- Ensamblaje de Novo con Trinity (ARRAY por GRUPO)
# ----------------------------------------------------------------------------
# Trinity es MUY intensivo en memoria y CPU. Estrategia:
#   - Pooled por grupo biologico (Ch14, EP08, Se17, SES) en lugar de por
#     muestra -> incrementa cobertura por transcrito y reduce 4 corridas vs 12.
#   - groups.tsv tiene (col1=grupo, col2=R1_csv, col3=R2_csv).
#   - Trinity acepta lecturas concatenadas con coma.
###############################################################################
#SBATCH --job-name=trinity
#SBATCH --partition=cicese
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --array=1-4                         # n grupos biologicos
#SBATCH --output=logs/trinity_%A_%a.out
#SBATCH --error=logs/trinity_%A_%a.err

set -euo pipefail

module load trinityrnaseq-v2.15.1 || true


PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"

WORDIR="${PROJECT}/workdir"
GROUPS="${WORDIR}/groups.tsv"

CLEANDIR="${WORDIR}/01_fastp"
OUTBASE="${WORDIR}/05_trinity"
mkdir -p "${OUTBASE}"

IFS=$'\t' read -r GRUPO R1_CSV R2_CSV < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "${GROUPS}")

# --- Sustituir rutas crudas por las limpias del paso 01 --------------------
# samples.tsv contiene los nombres limpios; re-derivamos los paths de fastp.
# Lista de sample_ids para este grupo:
SAMP_IDS=$(awk -F'\t' -v g="${GRUPO}" '$1==g{print $2}' \
            "${WORDIR}/samples.tsv")

LEFT=""; RIGHT=""
for S in ${SAMP_IDS}; do
    LEFT="${LEFT}${LEFT:+,}${CLEANDIR}/${S}_R1.clean.fq.gz"
    RIGHT="${RIGHT}${RIGHT:+,}${CLEANDIR}/${S}_R2.clean.fq.gz"
done

OUTDIR="${OUTBASE}/trinity_${GRUPO}"
mkdir -p "${OUTDIR}"

echo "[$(date +%F_%T)] Trinity ${GRUPO}"
echo "  LEFT : ${LEFT}"
echo "  RIGHT: ${RIGHT}"

Trinity \
    --seqType fq \
    --max_memory $(( SLURM_MEM_PER_NODE / 1024 - 16 ))G \
    --left  "${LEFT}" \
    --right "${RIGHT}" \
    --CPU "${SLURM_CPUS_PER_TASK}" \
    --SS_lib_type RF \
    --min_kmer_cov 2 \
    --normalize_max_read_cov 50 \
    --output "${OUTDIR}" \
    --full_cleanup

# --- Estadisticas N50 -------------------------------------------------------
TrinityStats.pl "${OUTDIR}.Trinity.fasta" > "${OUTDIR}.stats.txt"
echo "[OK] Trinity ${GRUPO} finalizado"
