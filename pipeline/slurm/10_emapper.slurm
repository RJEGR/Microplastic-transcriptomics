#!/bin/bash
###############################################################################
# 10_emapper.slurm  -- Anotacion funcional con eggNOG-mapper v2 (DIAMOND mode)
# ----------------------------------------------------------------------------
# Input: proteoma "consenso" del experimento =
#        (A) proteinas RefSeq filtradas por gene_ids expresados (paso 08)
#      + (B) proteinas novel predichas por TransDecoder         (paso 09)
#
# eggNOG-mapper genera:
#   - <prefix>.emapper.annotations           (tabla principal: GO, KEGG, COG, ...)
#   - <prefix>.emapper.seed_orthologs        (mejor hit ortologo)
#   - <prefix>.emapper.hits                  (DIAMOND raw)
#
# Modo DIAMOND (default --m diamond): rapido y suficiente para >95% de hits;
# si requieres mayor sensibilidad para huerfanos cambia a --m mmseqs o --m hmmer
# (este ultimo requiere mas memoria/tiempo).
#
# Dependencia: afterok:<jobid_transdecoder>  (09_transdecoder.slurm)
###############################################################################
#SBATCH --job-name=emapper
#SBATCH --partition=cicese
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rgomez@uabc.edu.mx
#SBATCH --output=logs/emapper_%j.out
#SBATCH --error=logs/emapper_%j.err

set -euo pipefail

module load conda-2024 || true

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/eggnog-mapper-master
export PATH=$PATH:$EXPORT

# ---------------------------- Rutas ---------------------------------------- #
PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"
WORDIR="${PROJECT}/workdir"

CDS_DIR="${WORDIR}/08_extract_cds"
TD_DIR="${WORDIR}/09_transdecoder"

REF_PEP="${CDS_DIR}/xbMagGiga1.1.expressed_ref.pep.fa"
NOVEL_PEP="${TD_DIR}/xbMagGiga1.1.novel_transcripts.fa.transdecoder.pep"

EGGNOG_DB="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Trinotate/TRINOTATE_DB/EGGNOG_DATA_DIR"

OUTDIR="${WORDIR}/10_emapper"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# ---------------------------- (1) Combinar proteomas ----------------------- #
COMBINED_PEP="${OUTDIR}/xbMagGiga1.1.expressed_proteome.pep.fa"

if [[ ! -s "${COMBINED_PEP}" ]]; then
    : > "${COMBINED_PEP}"
    if [[ -s "${REF_PEP}" ]]; then
        # Prefijo 'REF|' para distinguir trazabilidad en la tabla emapper.
        sed -E 's/^>/>REF|/' "${REF_PEP}" >> "${COMBINED_PEP}"
    else
        echo "[WARN] No existe ${REF_PEP} (no se ejecuto 08?)."
    fi

    if [[ -s "${NOVEL_PEP}" ]]; then
        # TransDecoder cabecera tipo: MSTRG.X.Y.p1 type:complete len:N gc:univ ...
        # Limpiar y prefijo 'NOV|'
        sed -E 's/^>([^ ]+).*/>NOV|\1/' "${NOVEL_PEP}" >> "${COMBINED_PEP}"
    else
        echo "[INFO] No hay proteinas novel (sin TransDecoder output)."
    fi

    # Sanea stops internos '*' que rompen DIAMOND
    sed -i 's/\*//g' "${COMBINED_PEP}"

    N_TOTAL=$(grep -c '^>' "${COMBINED_PEP}")
    echo "[OK] Proteoma combinado: ${N_TOTAL} secuencias -> ${COMBINED_PEP}"
else
    echo "[SKIP] ${COMBINED_PEP} ya existe."
fi

# ---------------------------- (2) eggNOG-mapper ---------------------------- #
PREFIX="${OUTDIR}/eggnog_xbMagGiga1.1"

if [[ -s "${PREFIX}.emapper.annotations" ]]; then
    echo "[SKIP] Anotacion eggnog ya existe (${PREFIX}.emapper.annotations)."
else
    emapper.py \
        -i "${COMBINED_PEP}" \
        --itype proteins \
        -m diamond \
        --cpu "${SLURM_CPUS_PER_TASK}" \
        --data_dir "${EGGNOG_DB}" \
        --output "$(basename "${PREFIX}")" \
        --output_dir "${OUTDIR}" \
        --temp_dir "${OUTDIR}/tmp" \
        --go_evidence non-electronic \
        --tax_scope auto \
        --target_orthologs all \
        --override
fi

# ---------------------------- (3) Quick stats ------------------------------ #
ANN="${PREFIX}.emapper.annotations"
if [[ -s "${ANN}" ]]; then
    TOTAL=$(grep -cv '^#' "${ANN}" || true)
    REF_HITS=$(grep -cE '^REF\|' "${ANN}"  || true)
    NOV_HITS=$(grep -cE '^NOV\|' "${ANN}"  || true)
    echo "============================================================"
    echo "[STATS] Anotaciones totales: ${TOTAL}"
    echo "[STATS]   - desde REF (RefSeq): ${REF_HITS}"
    echo "[STATS]   - desde NOV (novel) : ${NOV_HITS}"
    echo "============================================================"
fi

echo "[DONE] eggNOG-mapper terminado. Outputs en ${OUTDIR}/"
