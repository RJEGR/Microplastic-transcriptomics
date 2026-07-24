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
#SBATCH --job-name=emapper_mp
#SBATCH --partition=cicese
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/emapper_%j.out
#SBATCH --error=logs/emapper_%j.err
#
# ----------------------------------------------------------------------------
# Dimensionamiento y particion (2026-07-16, cluster hpc_omica CICESE):
#
#   Nodos observados: 24 CPU / 128 GB c/u. ~4 idle (nodo[8,10-11,13...]),
#   5+ en DOWN/DRAIN, 4 en ALLOC. Ver sinfo del dia.
#
#   Se cambio --partition=cicese -> d45 porque el usuario 'rgomez' bajo la
#   cuenta 'clarissa' SOLO tiene asociacion SLURM con d15 y d45 (verificado
#   con: sacctmgr show assoc user=$USER). En cualquier otra particion SLURM
#   rechaza con "Nodes required for job are DOWN, DRAINED or reserved for
#   jobs in higher priority partitions" -- mensaje enganoso: los nodos
#   existen pero la cuenta no tiene privilegio para usarlos ahi.
#     d15: MaxJobs=2 concurrentes, 15 dias timelimit.
#     d45: sin limite de jobs, 45 dias timelimit.  <- eleccion segura.
#
#   eggNOG-mapper -m diamond: ~4-6 GB/thread + .dmnd v5 ~10 GB.
#   Con 16 threads: RSS pico esperado 40-55 GB -> 64 GB deja ~15% de margen.
#
#   Diagnostico si sigue en PENDING:
#     squeue -u \$USER -o '%.10i %.10P %.20j %.8T %.10M %.20R'
#     sinfo -o '%.15N %.10P %.6t %.10c %.10m' | grep -E 'd30|idle'
#     sprio -j \$JOBID          # score de prioridad
#     scontrol show job \$JOBID # razon detallada del PENDING
# ----------------------------------------------------------------------------

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
TMPDIR_EMAPPER="${OUTDIR}/tmp"
mkdir -p "${OUTDIR}" "${TMPDIR_EMAPPER}"   # --temp_dir exige un dir EXISTENTE
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

    # Sanea caracteres invalidos SOLO en lineas de secuencia (no en headers):
    #   '*' = stop codon (TransDecoder los inserta en ORFs '5prime_partial')
    #   '.' = placeholder por codon ambiguo/incompleto (gffread + TransDecoder)
    # NOTA: usar '/^>/!' es critico -- los headers 'NOV|MSTRG.1.2.p1' llevan puntos
    # legitimos que NO deben tocarse; borrarlos rompe la trazabilidad al gene_id.
    sed -i '/^>/! s/[*.]//g' "${COMBINED_PEP}"

    N_TOTAL=$(grep -c '^>' "${COMBINED_PEP}")
    echo "[OK] Proteoma combinado: ${N_TOTAL} secuencias -> ${COMBINED_PEP}"
else
    # Reaplica la sanitizacion aunque el archivo ya exista: cubre reruns tras
    # cambios en la regex (evita quedarnos con un pep con '.' sin limpiar).
    sed -i '/^>/! s/[*.]//g' "${COMBINED_PEP}"
    echo "[SKIP] ${COMBINED_PEP} ya existe (sanitizacion re-aplicada por seguridad)."
fi

# ---------------------------- (2) eggNOG-mapper ---------------------------- #
PREFIX="${OUTDIR}/eggnog_xbMagGiga1.1"

SLURM_CPUS_PER_TASK=64

if [[ -s "${PREFIX}.emapper.annotations" ]]; then
    echo "[SKIP] Anotacion eggnog ya existe (${PREFIX}.emapper.annotations)."
else
    emapper.py \
        -i "${COMBINED_PEP}" \
        --itype proteins \
        -m diamond \
        --cpu "${SLURM_CPUS_PER_TASK:-16}" \
        --data_dir "${EGGNOG_DB}" \
        --output "$(basename "${PREFIX}")" \
        --output_dir "${OUTDIR}" \
        --temp_dir "${TMPDIR_EMAPPER}" \
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
