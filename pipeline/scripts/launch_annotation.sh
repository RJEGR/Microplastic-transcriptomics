#!/usr/bin/env bash
###############################################################################
# launch_annotation.sh -- Encadena 08_extract_cds -> 09_transdecoder -> 10_emapper
# ----------------------------------------------------------------------------
# Uso:
#   bash pipeline/scripts/launch_annotation.sh              # sin dependencia
#   bash pipeline/scripts/launch_annotation.sh <JID_MERGE>  # con dependencia
#
# Resuelve el error clasico "sbatch: Batch job submission failed: Job dependency
# problem" que ocurre cuando:
#   (a) la variable $JID_MERGE viene vacia / no numerica,
#   (b) el job de referencia ya termino en estado FAILED/CANCELLED (afterok
#       nunca se satisface -> SLURM aborta el submit),
#   (c) el job ya termino COMPLETED hace tiempo y salio del accounting.
#
# Estrategia: chequea el estado real con `sacct` y OMITE --dependency cuando
# ya no aplica; asi el submit no falla y el pipeline arranca inmediatamente.
###############################################################################
set -euo pipefail

SLURM_DIR="$(cd "$(dirname "$0")/../slurm" && pwd)"
cd "$(dirname "$SLURM_DIR")"   # pipeline/
mkdir -p logs

# ---------------------------- helpers -------------------------------------- #
_is_num() { [[ "${1:-}" =~ ^[0-9]+$ ]]; }

# Devuelve el estado principal del job (COMPLETED / RUNNING / PENDING / FAILED / ...)
_job_state() {
    local jid="$1"
    # sacct puede tardar unos segundos en registrar jobs recien enviados
    sacct -j "${jid}" -X -n -o State 2>/dev/null \
        | awk 'NR==1{gsub(/ /,""); print; exit}'
}

# Devuelve la flag --dependency=afterok:<jid> SOLO si tiene sentido usarla.
# Imprime cadena vacia (y avisa por stderr) cuando no aplica.
_dep_flag() {
    local jid="${1:-}"
    if ! _is_num "${jid}"; then
        echo "[WARN] JID de dependencia vacio o no numerico ('${jid}'); envio sin --dependency." >&2
        return 0
    fi

    local state; state="$(_job_state "${jid}")"
    case "${state}" in
        PENDING|RUNNING|CONFIGURING|SUSPENDED|REQUEUED|RESIZING)
            echo "--dependency=afterok:${jid}"
            ;;
        COMPLETED)
            echo "[INFO] Job ${jid} ya COMPLETED; envio sin --dependency." >&2
            ;;
        FAILED|CANCELLED*|TIMEOUT|NODE_FAIL|OUT_OF_MEMORY|BOOT_FAIL|DEADLINE|PREEMPTED)
            echo "[ERROR] Job ${jid} termino en estado '${state}'. afterok NUNCA se satisfara." >&2
            echo "        Revisa logs/ y relanza el paso previo antes de continuar." >&2
            exit 2
            ;;
        "" )
            echo "[WARN] No pude leer estado de job ${jid} via sacct; envio sin --dependency." >&2
            ;;
        *)
            echo "[INFO] Job ${jid} en estado '${state}'; encadeno afterok por si acaso." >&2
            echo "--dependency=afterok:${jid}"
            ;;
    esac
}

# ---------------------------- lanzamiento ---------------------------------- #
JID_MERGE="${1:-}"

DEP08="$(_dep_flag "${JID_MERGE}")"
J08=$(sbatch --parsable ${DEP08} "${SLURM_DIR}/08_extract_cds.slurm")
echo "[SUBMIT] 08_extract_cds  -> JobID ${J08}  (dep: ${DEP08:-none})"

DEP09="--dependency=afterok:${J08}"    # J08 acaba de ser encolado -> siempre valido
J09=$(sbatch --parsable ${DEP09} "${SLURM_DIR}/09_transdecoder.slurm")
echo "[SUBMIT] 09_transdecoder -> JobID ${J09}  (dep: ${DEP09})"

DEP10="--dependency=afterok:${J09}"
J10=$(sbatch --parsable ${DEP10} "${SLURM_DIR}/10_emapper.slurm")
echo "[SUBMIT] 10_emapper      -> JobID ${J10}  (dep: ${DEP10})"

echo ""
echo "============================================================"
echo "Monitorea con:  squeue -u \$USER -o '%.10i %.20j %.8T %.10M %R'"
echo "                sacct  -j ${J08},${J09},${J10} --format=JobID,JobName,State,Elapsed"
echo "============================================================"
