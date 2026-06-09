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
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/hisat2_build_%j.out
#SBATCH --error=logs/hisat2_build_%j.err

set -euo pipefail
module load hisat2/2.2.1 || true

PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"
REF_GZ="${PROJECT}/Reference_dir/GCF_963853765.1_xbMagGiga1.1_genomic.fna.gz"
GTF_GZ="${PROJECT}/Reference_dir/GCF_963853765.1_xbMagGiga1.1_genomic.gtf.gz"
IDX_DIR="${PROJECT}/Reference_dir/hisat2_index"
PREFIX="xbMagGiga1.1"

mkdir -p "${IDX_DIR}"
cd "${IDX_DIR}"

# --- 1) Decompresion del FASTA (HISAT2-build requiere texto plano) ---------
REF_FA="${IDX_DIR}/${PREFIX}.fa"
[[ -s "${REF_FA}" ]] || zcat "${REF_GZ}" > "${REF_FA}"

# --- 2) Sitios de splicing y exones para mejorar la sensibilidad -----------
GTF_TXT="${IDX_DIR}/${PREFIX}.gtf"
[[ -s "${GTF_TXT}" ]] || zcat "${GTF_GZ}" > "${GTF_TXT}"

hisat2_extract_splice_sites.py "${GTF_TXT}" > "${IDX_DIR}/${PREFIX}.ss"
hisat2_extract_exons.py        "${GTF_TXT}" > "${IDX_DIR}/${PREFIX}.exon"

# --- 3) Construccion del indice -------------------------------------------
hisat2-build \
    -p "${SLURM_CPUS_PER_TASK}" \
    --ss   "${IDX_DIR}/${PREFIX}.ss" \
    --exon "${IDX_DIR}/${PREFIX}.exon" \
    "${REF_FA}" \
    "${IDX_DIR}/${PREFIX}"

echo "[OK] Indice HISAT2 generado en ${IDX_DIR}/${PREFIX}.*.ht2"
