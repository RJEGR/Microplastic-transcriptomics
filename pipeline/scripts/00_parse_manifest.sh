#!/usr/bin/env bash
###############################################################################
# 00_parse_manifest.sh
# ----------------------------------------------------------------------------
# Lee el archivo Manifest.tsv (5 columnas: idx, grupo, replica, R1, R2) y
# genera dos artefactos:
#   1) samples.tsv  -> tabla limpia (grupo, sample_id, R1, R2)  para los
#                      array jobs de SLURM (lectura por SLURM_ARRAY_TASK_ID).
#   2) groups.tsv   -> tabla colapsada por grupo biologico (Ch14, EP08, Se17,
#                      SES) usada por Trinity y rnaSPAdes en modo pooled.
#
# Uso:
#   bash 00_parse_manifest.sh <ruta/Manifest.tsv> <ruta/output_dir>
#
# Salida esperada (samples.tsv):
#   grupo<TAB>sample_id<TAB>R1<TAB>R2
###############################################################################
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${1:-${SCRIPT_DIR}/../../Manifest.tsv}"
OUTDIR="${2:-./}"
mkdir -p "${OUTDIR}"

SAMPLES="${OUTDIR}/samples.tsv"
GROUPS="${OUTDIR}/groups.tsv"

# --- Validacion ------------------------------------------------------------
[[ -s "${MANIFEST}" ]] || { echo "ERROR: Manifest vacio o inexistente: ${MANIFEST}"; exit 1; }

# --- samples.tsv : una linea por libreria paired-end ------------------------
awk -F'\t' 'BEGIN{OFS="\t"}
    # Col1=idx, Col2=grupo, Col3=sample_id, Col4=R1, Col5=R2
    NF>=5 && $4 ~ /\.f(ast)?q\.gz$/ && $5 ~ /\.f(ast)?q\.gz$/ {
        print $2, $3, $4, $5
    }' "${MANIFEST}" > "${SAMPLES}"

N=$(wc -l < "${SAMPLES}")
echo "[parse] ${N} librerias paired-end registradas en ${SAMPLES}"

# --- groups.tsv : concatenacion comma-separated R1/R2 por grupo -------------
# Necesario para Trinity --left A,B,C --right A',B',C' y rnaSPAdes pooled
awk -F'\t' 'BEGIN{OFS="\t"}
    { r1[$1]=(r1[$1]?r1[$1]",":"")$3; r2[$1]=(r2[$1]?r2[$1]",":"")$4 }
    END { for (g in r1) print g, r1[g], r2[g] }' "${SAMPLES}" \
  | sort > "${GROUPS}"

NG=$(wc -l < "${GROUPS}")
echo "[parse] ${NG} grupos biologicos consolidados en ${GROUPS}"

# --- Reporte rapido --------------------------------------------------------
printf "\n--- samples.tsv (preview) ---\n"
column -t -s $'\t' "${SAMPLES}" | head
printf "\n--- groups.tsv  (preview) ---\n"
column -t -s $'\t' "${GROUPS}" | head
