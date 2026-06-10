#!/usr/bin/env bash
###############################################################################
# 00_parse_manifest.sh
# ----------------------------------------------------------------------------
# Reads Manifest.tsv produced by 00_parse_filenames.sh
# (columns: idx, grupo, replica, R1, R2) and writes two artifacts:
#
#   1) samples.tsv  -> grupo, replica, R1, R2  (one row per library;
#                      indexed by SLURM_ARRAY_TASK_ID in array jobs)
#   2) groups.tsv   -> grupo, R1_csv, R2_csv   (paths joined per biological
#                      group for Trinity --left/--right and rnaSPAdes pooled)
#
# Usage:
#   bash 00_parse_manifest.sh <Manifest.tsv> [workdir]
#   bash 00_parse_filenames.sh < filenames.txt | bash 00_parse_manifest.sh -
###############################################################################
set -euo pipefail

MANIFEST="${1:--}"          # accepts a file path or '-' for stdin
WORKDIR="${2:-${PWD}/workdir}"

mkdir -p "${WORKDIR}"
SAMPLES="${WORKDIR}/samples.tsv"
GROUPS_TSV="${WORKDIR}/groups.tsv"

if [[ "${MANIFEST}" != "-" ]]; then
    [[ -s "${MANIFEST}" ]] || { echo "ERROR: Manifest vacio o inexistente: ${MANIFEST}"; exit 1; }
fi

# --- samples.tsv : one row per paired-end library ---------------------------
# Input cols: idx(1) grupo(2) replica(3) R1(4) R2(5); skip header row
awk -F'\t' 'BEGIN{OFS="\t"}
    NR==1 { next }
    NF>=5 && $4 ~ /\.f(ast)?q\.gz$/ && $5 ~ /\.f(ast)?q\.gz$/ {
        print $2, $3, $4, $5
    }' "${MANIFEST}" > "${SAMPLES}"

N=$(wc -l < "${SAMPLES}")
echo "[parse] ${N} paired-end libraries -> ${SAMPLES}"

# --- groups.tsv : comma-joined R1/R2 paths per biological group -------------
# Required for Trinity --left A,B,C --right X,Y,Z and rnaSPAdes pooled mode
awk -F'\t' 'BEGIN{OFS="\t"}
    { r1[$1]=(r1[$1]?r1[$1]",":"")$3; r2[$1]=(r2[$1]?r2[$1]",":"")$4 }
    END { for (g in r1) print g, r1[g], r2[g] }' "${SAMPLES}" \
  | sort > "${GROUPS_TSV}"

NG=$(wc -l < "${GROUPS_TSV}")
echo "[parse] ${NG} biological groups -> ${GROUPS}"

# --- Quick preview ----------------------------------------------------------
printf "\n--- samples.tsv ---\n"
column -t -s $'\t' "${SAMPLES}" | head
printf "\n--- groups.tsv ---\n"
column -t -s $'\t' "${GROUPS_TSV}" | head
