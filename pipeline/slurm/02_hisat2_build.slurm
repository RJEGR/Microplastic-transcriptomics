#!/bin/bash
###############################################################################
# 02_hisat2_build.slurm  -- Construye indice HISAT2 (UNICO, no array)
# ----------------------------------------------------------------------------
# Genoma de referencia: GCF_963853765.1 xbMagGiga1.1 (assembly cromosomal).
# IMPORTANTE: este indice es 10x mas eficiente si se hace ANTES de los array
# jobs de alineamiento. Usar --dependency=afterok:<jobid> al lanzar 03.
###############################################################################
#SBATCH --job-name=hisat2_build
#SBATCH --partition=cicese
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/hisat2_build_%j.out
#SBATCH --error=logs/hisat2_build_%j.err

set -euo pipefail

#module load hisat2/2.2.1 || true

EXPORT=/LUSTRE/apps/bioinformatica/hisat2/
export PATH=$PATH:$EXPORT


PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"
REF_GZ="${PROJECT}/Reference_dir/GCF_963853765.1_xbMagGiga1.1_genomic.fna.gz"
GTF_GZ="${PROJECT}/Reference_dir/GCF_963853765.1_xbMagGiga1.1_genomic.gtf.gz"
IDX_DIR="${PROJECT}/Reference_dir/hisat2_index"
PREFIX="xbMagGiga1.1"

# issue w/ stringtie Error: no valid ID found for GFF record
# cp xbMagGiga1.1.gtf xbMagGiga1.1.gtf.bkp
# gffread xbMagGiga1.1.gtf.bkp -T -o xbMagGiga1.1.gtf OR
# grep -v 'transcript_id ""' xbMagGiga1.1.gtf.bkp > xbMagGiga1.1.gtf

mkdir -p "${IDX_DIR}"
cd "${IDX_DIR}"

# --- 1) Decompresion del FASTA (HISAT2-build requiere texto plano) ---------
REF_FA="${IDX_DIR}/${PREFIX}.fa"
[[ -s "${REF_FA}" ]] || zcat "${REF_GZ}" > "${REF_FA}"

# --- 2) Sitios de splicing y exones para mejorar la sensibilidad -----------
GTF_TXT="${IDX_DIR}/${PREFIX}.gtf"
[[ -s "${GTF_TXT}" ]] || zcat "${GTF_GZ}" | grep -v 'transcript_id ""' > "${GTF_TXT}"

SS_FILE="${IDX_DIR}/${PREFIX}.ss"
EXON_FILE="${IDX_DIR}/${PREFIX}.exon"

if [[ -s "${SS_FILE}" && -s "${EXON_FILE}" ]]; then
    echo "[SKIP] Archivos .ss y .exon ya existen, omitiendo extraccion"
else
    hisat2_extract_splice_sites.py "${GTF_TXT}" > "${SS_FILE}"
    hisat2_extract_exons.py        "${GTF_TXT}" > "${EXON_FILE}"
fi

# --- 3) Construccion del indice -------------------------------------------
if ls "${IDX_DIR}/${PREFIX}".*.ht2 &>/dev/null; then
    echo "[SKIP] Indice HISAT2 ya existe en ${IDX_DIR}/${PREFIX}.*.ht2, omitiendo build"
else
    hisat2-build \
        -p "${SLURM_CPUS_PER_TASK}" \
        --ss   "${SS_FILE}" \
        --exon "${EXON_FILE}" \
        "${REF_FA}" \
        "${IDX_DIR}/${PREFIX}"
    echo "[OK] Indice HISAT2 generado en ${IDX_DIR}/${PREFIX}.*.ht2"
fi
