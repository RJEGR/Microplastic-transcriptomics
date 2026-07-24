#!/bin/bash
###############################################################################
# 09_transdecoder.slurm  -- Predice ORFs/proteinas en transcritos NOVEL
# ----------------------------------------------------------------------------
# Solo aplica al subset 'novel' del merged.gtf (MSTRG.* sin ref_gene_id).
# Los CDS de RefSeq (paso 08) ya vienen traducidos por gffread y no necesitan
# prediccion de ORFs.
#
# Protocolo TransDecoder v5.7.0 estandar:
#   1) TransDecoder.LongOrfs    (ORFs >=100 aa por defecto)
#   2) TransDecoder.Predict     (--single_best_only para deduplicar)
#
# Dependencia: afterok:<jobid_extract_cds>  (08_extract_cds.slurm)
###############################################################################
#SBATCH --job-name=transdecoder_mp
#SBATCH --partition=d45
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rgomez41@uabc.edu.mx
#SBATCH --output=logs/transdecoder_%j.out
#SBATCH --error=logs/transdecoder_%j.err
# Nota: la cuenta 'clarissa' solo tiene asociacion SLURM con d15 y d45.

set -euo pipefail

EXPORT=/LUSTRE/apps/bioinformatica/TransDecoder-v5.7.0/
export PATH=$PATH:$EXPORT

PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"
WORDIR="${PROJECT}/workdir"

CDS_DIR="${WORDIR}/08_extract_cds"
NOVEL_FA="${CDS_DIR}/xbMagGiga1.1.novel_transcripts.fa"

OUTDIR="${WORDIR}/09_transdecoder"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# Early exit si no hay novel
if [[ ! -s "${NOVEL_FA}" ]]; then
    echo "[INFO] No hay novel transcripts (${NOVEL_FA} ausente o vacio). Nada que predecir."
    exit 0
fi

LOCAL_FA="${OUTDIR}/$(basename "${NOVEL_FA}")"
[[ -e "${LOCAL_FA}" ]] || ln -s "${NOVEL_FA}" "${LOCAL_FA}"

# --- 1) LongOrfs ----------------------------------------------------------- #
if [[ -d "${LOCAL_FA}.transdecoder_dir" ]]; then
    echo "[SKIP] LongOrfs ya ejecutado."
else
    TransDecoder.LongOrfs -t "${LOCAL_FA}" --output_dir "${LOCAL_FA}.transdecoder_dir"
fi

# --- 2) Predict ------------------------------------------------------------ #
# --single_best_only: una ORF por transcrito (la mejor); evita ORFs anidadas
# que inflarian el input de eggnog-mapper y duplicarian hits.
FINAL_PEP="${LOCAL_FA}.transdecoder.pep"

if [[ -s "${FINAL_PEP}" ]]; then
    echo "[SKIP] ${FINAL_PEP} ya existe."
else
    TransDecoder.Predict \
        -t "${LOCAL_FA}" \
        --single_best_only \
        --output_dir "${LOCAL_FA}.transdecoder_dir"
fi

# --- 3) Resumen ------------------------------------------------------------ #
N_PEP=$(grep -c '^>' "${FINAL_PEP}" || echo 0)
echo "[OK] Proteinas novel predichas: ${N_PEP}"
echo "[OUT] ${FINAL_PEP}"
echo "[OUT] ${LOCAL_FA}.transdecoder.cds"
echo "[OUT] ${LOCAL_FA}.transdecoder.gff3"
echo "[OUT] ${LOCAL_FA}.transdecoder.bed"
