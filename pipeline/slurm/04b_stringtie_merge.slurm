#!/bin/bash
###############################################################################
# 04b_stringtie_merge.slurm  -- Consolida un GTF no-redundante a partir de
# todos los GTFs por muestra. UNICO (no array).
# Dependencia: afterok:<jobid_stringtie_array>
###############################################################################
#SBATCH --job-name=stie_merge
#SBATCH --partition=cicese
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=logs/stringtie_merge_%j.out
#SBATCH --error=logs/stringtie_merge_%j.err

set -euo pipefail
module load stringtie/2.2.1 gffcompare/0.12.6 || true

PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"
PER_SAMP="${PROJECT}/04_stringtie/per_sample"
MERGEDIR="${PROJECT}/04_stringtie/merged"
GTF_REF="${PROJECT}/Reference_dir/hisat2_index/xbMagGiga1.1.gtf"
mkdir -p "${MERGEDIR}"

# Lista de GTFs por muestra
ls "${PER_SAMP}"/*.gtf > "${MERGEDIR}/gtf_list.txt"

# --- Merge -----------------------------------------------------------------
stringtie --merge \
    -G "${GTF_REF}" \
    -o "${MERGEDIR}/merged.gtf" \
    "${MERGEDIR}/gtf_list.txt"

# --- Compara assembly merged vs anotacion oficial (sensibilidad isoformas) -
gffcompare -r "${GTF_REF}" -o "${MERGEDIR}/gffcmp" "${MERGEDIR}/merged.gtf"

echo "[OK] merged.gtf con $(grep -c $'\ttranscript\t' "${MERGEDIR}/merged.gtf") transcritos"
