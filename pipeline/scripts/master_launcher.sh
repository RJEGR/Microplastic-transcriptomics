#!/usr/bin/env bash
###############################################################################
# master_launcher.sh  -- Orquestador del pipeline con dependencias SLURM
# ----------------------------------------------------------------------------
# Encadena: parse_manifest -> fastp -> [hisat2_build || trinity || rnaspades]
#                                  \-> hisat2_align -> stringtie -> merge
#                                                                -> quant
#                                                              \-> multiqc
#
# Uso:
#   bash master_launcher.sh /ruta/Manifest.tsv
###############################################################################
set -euo pipefail

PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"
WORDIR="${PROJECT}/workdir"

SLURM_DIR="${PROJECT}/pipeline/slurm"
SCRIPTS="${PROJECT}/pipeline/scripts"
LOGS="${PROJECT}/pipeline/logs"
MANIFEST="${1:-${PROJECT}/Manifest.tsv}"

mkdir -p "${LOGS}"
cd "${PROJECT}/pipeline"

echo "===================================================================="
echo "MICROPLASTICS RNA-Seq PIPELINE  | $(date)"
echo "Manifest: ${MANIFEST}"
echo "===================================================================="

# -- 0) Parsing local (no SLURM, segundos) ----------------------------------

bash "${SCRIPTS}/00_parse_manifest.sh" ../"${MANIFEST}" "${WORDIR}"

NSAMP=$(wc -l < "${WORDIR}/samples.tsv")
NGRP=$(wc -l < "${WORDIR}/groups.tsv")

echo "  -> ${NSAMP} muestras, ${NGRP} grupos"

# Ajustar dinamicamente el rango del array si difiere de 12
sed -i.bak "s|^#SBATCH --array=1-[0-9]*%[0-9]*|#SBATCH --array=1-${NSAMP}%6|" \
    "${SLURM_DIR}/01_fastp.slurm" \
    "${SLURM_DIR}/03_hisat2_align.slurm" \
    "${SLURM_DIR}/04_stringtie.slurm" \
    "${SLURM_DIR}/04c_stringtie_quant.slurm"
sed -i.bak "s|^#SBATCH --array=1-[0-9]*|#SBATCH --array=1-${NGRP}|" \
    "${SLURM_DIR}/05_trinity.slurm" \
    "${SLURM_DIR}/06_rnaspades.slurm"

# -- 1) FASTP --------------------------------------------------------------
JID_FASTP=$(sbatch --parsable "${SLURM_DIR}/01_fastp.slurm")
echo "  [submit] fastp           jobid=${JID_FASTP}"

# -- 2) HISAT2 build (independiente de fastp) ------------------------------
JID_IDX=$(sbatch --parsable "${SLURM_DIR}/02_hisat2_build.slurm")
echo "  [submit] hisat2_build    jobid=${JID_IDX}"

# -- 3) HISAT2 align (depende de fastp Y de indice) ------------------------
JID_ALN=$(sbatch --parsable \
    --dependency=afterok:${JID_FASTP}:${JID_IDX} \
    "${SLURM_DIR}/03_hisat2_align.slurm")
echo "  [submit] hisat2_align    jobid=${JID_ALN}    [dep ${JID_FASTP},${JID_IDX}]"

# -- 4) StringTie per-sample -----------------------------------------------
JID_STIE=$(sbatch --parsable --dependency=afterok:${JID_ALN} \
    "${SLURM_DIR}/04_stringtie.slurm")
echo "  [submit] stringtie       jobid=${JID_STIE}   [dep ${JID_ALN}]"

# -- 4b) StringTie merge ---------------------------------------------------
JID_MRG=$(sbatch --parsable --dependency=afterok:${JID_STIE} \
    "${SLURM_DIR}/04b_stringtie_merge.slurm")
echo "  [submit] stringtie_mrg   jobid=${JID_MRG}    [dep ${JID_STIE}]"

# -- 4c) Re-quantification para DESeq2 -------------------------------------
JID_QNT=$(sbatch --parsable --dependency=afterok:${JID_MRG} \
    "${SLURM_DIR}/04c_stringtie_quant.slurm")
echo "  [submit] stringtie_qnt   jobid=${JID_QNT}    [dep ${JID_MRG}]"

# -- 5) Trinity de novo (depende de fastp) ---------------------------------
JID_TRI=$(sbatch --parsable --dependency=afterok:${JID_FASTP} \
    "${SLURM_DIR}/05_trinity.slurm")
echo "  [submit] trinity         jobid=${JID_TRI}    [dep ${JID_FASTP}]"

# -- 6) rnaSPAdes de novo (depende de fastp) -------------------------------
JID_SPD=$(sbatch --parsable --dependency=afterok:${JID_FASTP} \
    "${SLURM_DIR}/06_rnaspades.slurm")
echo "  [submit] rnaspades       jobid=${JID_SPD}    [dep ${JID_FASTP}]"

# -- 7) MultiQC final ------------------------------------------------------
JID_MQC=$(sbatch --parsable \
    --dependency=afterok:${JID_QNT}:${JID_ALN} \
    "${SLURM_DIR}/07_multiqc.slurm")
echo "  [submit] multiqc         jobid=${JID_MQC}    [dep ${JID_QNT},${JID_ALN}]"

echo "===================================================================="
echo "Cola enviada. Monitorear con:"
echo "   squeue -u \$USER -o '%.10i %.20j %.8T %.10M %R'"
echo "   sacct  -j ${JID_FASTP},${JID_ALN},${JID_STIE},${JID_MRG},${JID_TRI},${JID_SPD}"
echo "===================================================================="
