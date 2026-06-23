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
#SBATCH --time=06:00:00
#SBATCH --output=logs/stringtie_merge_%j.out
#SBATCH --error=logs/stringtie_merge_%j.err

set -euo pipefail
#module load stringtie/2.2.1 gffcompare/0.12.6 || true

EXPORT=/LUSTRE/apps/bioinformatica/stringtie/
export PATH=$PATH:$EXPORT

#EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/apps/gffread
#export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/apps/gffcompare
export PATH=$PATH:$EXPORT


PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"
WORDIR="${PROJECT}/workdir"

PER_SAMP="${WORDIR}/04_stringtie/per_sample"
MERGEDIR="${WORDIR}/04_stringtie/merged"
GTF_REF="${PROJECT}/Reference_dir/hisat2_index/xbMagGiga1.1.gtf"
mkdir -p "${MERGEDIR}"

MERGED_GTF="${MERGEDIR}/merged.gtf"

if [[ -s "${MERGED_GTF}" ]]; then
    echo "[skip] ${MERGED_GTF} ya existe"
else
    # Lista de GTFs por muestra
    ls "${PER_SAMP}"/*.gtf > "${MERGEDIR}/gtf_list.txt"

    # --- Merge ---------------------------------------------------------------
    stringtie --merge \
        -G "${GTF_REF}" \
        -o "${MERGED_GTF}" \
        "${MERGEDIR}/gtf_list.txt"

    # --- Compara assembly merged vs anotacion oficial (sensibilidad isoformas)
    gffcompare -r "${GTF_REF}" -o "${MERGEDIR}/gffcmp" "${MERGED_GTF}"

    echo "[OK] merged.gtf con $(grep -c $'\ttranscript\t' "${MERGED_GTF}") transcritos"
fi
