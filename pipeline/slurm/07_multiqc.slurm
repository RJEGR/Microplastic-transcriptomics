#!/bin/bash
###############################################################################
# 07_multiqc.slurm -- Reporte unificado QC + alineamiento
# Dependencia: afterok:<fastp_array>,<hisat2_align_array>
###############################################################################
#SBATCH --job-name=multiqc
#SBATCH --partition=cicese
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=logs/multiqc_%j.out
#SBATCH --error=logs/multiqc_%j.err

set -euo pipefail

#module load multiqc/1.25 || true

WDM=/LUSTRE/apps/Anaconda/2023/miniconda3/bin/

#EXPORT=/LUSTRE/apps/Anaconda/2023/miniconda3/bin/
#export PATH=$PATH:$EXPORT

PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"
WORDIR="${PROJECT}/workdir"

REPORT_DIR="${WORDIR}/00_multiqc"
mkdir -p "${REPORT_DIR}"

$WDM/multiqc \
    "${WORDIR}/01_fastp/qc_reports" \
    "${WORDIR}/03_hisat2" \
    "${WORDIR}/04_stringtie/merged" \
    -o "${REPORT_DIR}" \
    -n "microplastics_multiqc" \
    --force
