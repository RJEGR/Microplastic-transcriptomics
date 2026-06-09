#!/bin/bash
###############################################################################
# 07_multiqc.slurm -- Reporte unificado QC + alineamiento
# Dependencia: afterok:<fastp_array>,<hisat2_align_array>
###############################################################################
#SBATCH --job-name=multiqc
#SBATCH --partition=cicese
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/multiqc_%j.out
#SBATCH --error=logs/multiqc_%j.err

set -euo pipefail
module load multiqc/1.25 || true

PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"
REPORT_DIR="${PROJECT}/00_multiqc"
mkdir -p "${REPORT_DIR}"

multiqc \
    "${PROJECT}/01_fastp/qc_reports" \
    "${PROJECT}/03_hisat2" \
    "${PROJECT}/04_stringtie/merged" \
    -o "${REPORT_DIR}" \
    -n "microplastics_multiqc" \
    --force
