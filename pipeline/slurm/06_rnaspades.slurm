#!/bin/bash
###############################################################################
# 06_rnaspades.slurm  -- Ensamblaje de Novo con rnaSPAdes (ARRAY por GRUPO)
# ----------------------------------------------------------------------------
# rnaSPAdes es complementario a Trinity:
#   - Mas rapido y menos memoria
#   - Mejor recuperacion de isoformas raras / con bajo soporte
#   - Permite k-mer multiple (default: -k 25,55 para reads 150 pb)
###############################################################################
#SBATCH --job-name=rnaspades
#SBATCH --partition=d30
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --time=30:00:00
#SBATCH --array=1-4
#SBATCH --output=logs/rnaspades_%A_%a.out
#SBATCH --error=logs/rnaspades_%A_%a.err

set -euo pipefail
#module load spades/4.0.0 || true     # rnaSPAdes viene en SPAdes >=3.13

EXPORT=/LUSTRE/apps/bioinformatica/SPAdes-3.15.5-Linux/bin/
export PATH=$PATH:$EXPORT

PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"

WORDIR="${PROJECT}/workdir"
GROUPS_TSV="${WORDIR}/groups.tsv"
CLEANDIR="${WORDIR}/01_fastp"
OUTBASE="${WORDIR}/06_rnaspades"
mkdir -p "${OUTBASE}"

IFS=$'\t' read -r GRUPO R1_CSV R2_CSV < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "${GROUPS_TSV}")

# --- Construir argumentos --pe1-1 / --pe1-2 para cada libreria del grupo ---
SAMP_IDS=$(awk -F'\t' -v g="${GRUPO}" '$1==g{print $2}' \
            "${WORDIR}/samples.tsv")

ARGS=()
i=1
for S in ${SAMP_IDS}; do
    ARGS+=( "--pe${i}-1" "${CLEANDIR}/${S}_R1.clean.fq.gz" )
    ARGS+=( "--pe${i}-2" "${CLEANDIR}/${S}_R2.clean.fq.gz" )
    i=$((i+1))
done

OUTDIR="${OUTBASE}/rnaspades_${GRUPO}"
mkdir -p "${OUTDIR}"

echo "[$(date +%F_%T)] rnaSPAdes ${GRUPO}"
rnaspades.py \
    "${ARGS[@]}" \
    -t "${SLURM_CPUS_PER_TASK}" \
    -m $(( SLURM_MEM_PER_NODE / 1024 - 8 )) \
    -k 25,45,65 \
    --ss rf \
    -o "${OUTDIR}"

echo "[OK] rnaSPAdes ${GRUPO} -> ${OUTDIR}/transcripts.fasta"
