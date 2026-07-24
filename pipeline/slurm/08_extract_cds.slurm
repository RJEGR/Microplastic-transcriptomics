#!/bin/bash
###############################################################################
# 08_extract_cds.slurm  -- Extrae CDS/proteinas del genoma xbMagGiga1.1
# filtrando por los gene_ids expresados en el ensamblaje guiado de StringTie.
# ----------------------------------------------------------------------------
# Estrategia de dos vias (UNICO, no array):
#
#   (A) Conocidos: gene_ids del merged.gtf que conservan 'ref_gene_id' (anotacion
#       RefSeq NCBI). Se filtra la GTF de referencia y se usa gffread para
#       extraer CDS coordinadas + traduccion oficial.
#
#   (B) Novel: transcritos MSTRG.* sin 'ref_gene_id' (descubiertos por StringTie).
#       Se extrae la secuencia spliced de mRNA via gffread; las ORFs se
#       predicen en 09_transdecoder.slurm.
#
# El subset (A+B novel translated) se concatena en 10_emapper.slurm para anotar
# con eggNOG-mapper solo lo que realmente se expresa en las muestras 2022/2024.
#
# Dependencia: afterok:<jobid_stringtie_merge>  (04b_stringtie_merge.slurm)
###############################################################################
#SBATCH --job-name=extract_cds_mp
#SBATCH --partition=d45
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rgomez41@uabc.edu.mx
#SBATCH --output=logs/extract_cds_%j.out
#SBATCH --error=logs/extract_cds_%j.err
# Nota: la cuenta 'clarissa' solo tiene asociacion SLURM con d15 (MaxJobs=2)
# y d45 (sin limite). Ver sacctmgr show assoc user=$USER.

set -euo pipefail

# ---------------------------- Herramientas --------------------------------- #
EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/apps/gffread
export PATH=$PATH:$EXPORT

# seqkit (filtrado por ID) ya disponible en el cluster; sino:
# EXPORT=/LUSTRE/apps/bioinformatica/seqkit/
# export PATH=$PATH:$EXPORT

# ---------------------------- Rutas ---------------------------------------- #
PROJECT="/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/microplastics"
WORDIR="${PROJECT}/workdir"

REF_DIR="${PROJECT}/Reference_dir/hisat2_index"
GENOME_FA="${REF_DIR}/xbMagGiga1.1.fa"
GTF_REF="${REF_DIR}/xbMagGiga1.1.gtf"

MERGED_GTF="${WORDIR}/04_stringtie/merged/merged.gtf"

OUTDIR="${WORDIR}/08_extract_cds"
mkdir -p "${OUTDIR}"

# ---------------------------- Sanity checks -------------------------------- #
for f in "${GENOME_FA}" "${GTF_REF}" "${MERGED_GTF}"; do
    [[ -s "${f}" ]] || { echo "[ERROR] No existe o esta vacio: ${f}"; exit 1; }
done

# ---------------------------- (1) Listas de IDs ---------------------------- #
# Solo lineas tipo 'transcript' del merged.gtf.
EXPRESSED_REF_GENES="${OUTDIR}/expressed_ref_geneids.txt"
NOVEL_TX_IDS="${OUTDIR}/novel_transcriptids.txt"
EXPRESSED_REF_TX="${OUTDIR}/expressed_ref_transcriptids.txt"

awk -F'\t' '$3=="transcript"' "${MERGED_GTF}" \
    | grep -oP 'ref_gene_id "\K[^"]+' \
    | sort -u > "${EXPRESSED_REF_GENES}"

awk -F'\t' '$3=="transcript"' "${MERGED_GTF}" \
    | awk '/ref_gene_id/ {match($0,/transcript_id "([^"]+)"/,a); print a[1]}' \
    | sort -u > "${EXPRESSED_REF_TX}"

awk -F'\t' '$3=="transcript"' "${MERGED_GTF}" \
    | awk '!/ref_gene_id/ {match($0,/transcript_id "([^"]+)"/,a); print a[1]}' \
    | sort -u > "${NOVEL_TX_IDS}"

N_REF_GENES=$(wc -l < "${EXPRESSED_REF_GENES}")
N_REF_TX=$(wc -l   < "${EXPRESSED_REF_TX}")
N_NOVEL=$(wc -l    < "${NOVEL_TX_IDS}")
echo "[INFO] Genes RefSeq expresados:      ${N_REF_GENES}"
echo "[INFO] Transcritos RefSeq expresados: ${N_REF_TX}"
echo "[INFO] Transcritos novel (MSTRG.*):   ${N_NOVEL}"

# ---------------------------- (2) CDS/Prot REFERENCIA (filtrados) ---------- #
# gffread acepta --ids para filtrar por transcript_id.
REF_CDS_FA="${OUTDIR}/xbMagGiga1.1.expressed_ref.cds.fa"
REF_PEP_FA="${OUTDIR}/xbMagGiga1.1.expressed_ref.pep.fa"
REF_GTF_SUB="${OUTDIR}/xbMagGiga1.1.expressed_ref.gtf"

if [[ ! -s "${REF_PEP_FA}" ]]; then
    gffread "${GTF_REF}" \
        --ids "${EXPRESSED_REF_TX}" \
        -g "${GENOME_FA}" \
        -T  -o "${REF_GTF_SUB}" \
        -x "${REF_CDS_FA}" \
        -y "${REF_PEP_FA}"

    # Sanea cabeceras: gffread agrega coordenadas/atributos; deja solo
    # transcript_id para trazabilidad rio abajo (emapper, DE).
    sed -i -E 's/^(>[^ ]+).*/\1/' "${REF_CDS_FA}" "${REF_PEP_FA}"

    echo "[OK] Reference CDS:     $(grep -c '^>' "${REF_CDS_FA}")"
    echo "[OK] Reference proteins: $(grep -c '^>' "${REF_PEP_FA}")"
else
    echo "[SKIP] ${REF_PEP_FA} ya existe."
fi

# ---------------------------- (3) mRNA NOVEL (spliced) --------------------- #
# Para los novel (sin CDS anotado) extraemos la secuencia spliced; las ORFs
# se predicen en 09_transdecoder.slurm.
NOVEL_FA="${OUTDIR}/xbMagGiga1.1.novel_transcripts.fa"
NOVEL_GTF="${OUTDIR}/xbMagGiga1.1.novel_transcripts.gtf"

if [[ ! -s "${NOVEL_FA}" && "${N_NOVEL}" -gt 0 ]]; then
    gffread "${MERGED_GTF}" \
        --ids "${NOVEL_TX_IDS}" \
        -g "${GENOME_FA}" \
        -T  -o "${NOVEL_GTF}" \
        -w "${NOVEL_FA}"

    sed -i -E 's/^(>[^ ]+).*/\1/' "${NOVEL_FA}"
    echo "[OK] Novel transcripts: $(grep -c '^>' "${NOVEL_FA}")"
else
    [[ "${N_NOVEL}" -eq 0 ]] && echo "[INFO] No hay transcritos novel para extraer."
    [[ -s "${NOVEL_FA}" ]]    && echo "[SKIP] ${NOVEL_FA} ya existe."
fi

# ---------------------------- (4) Manifest --------------------------------- #
cat > "${OUTDIR}/README.txt" <<EOF
Inputs:
  - Genoma:        ${GENOME_FA}
  - GTF ref RefSeq: ${GTF_REF}
  - GTF merged:    ${MERGED_GTF}

Outputs:
  - ${EXPRESSED_REF_GENES}  (${N_REF_GENES} gene_ids RefSeq expresados)
  - ${EXPRESSED_REF_TX}     (${N_REF_TX} transcript_ids RefSeq expresados)
  - ${NOVEL_TX_IDS}         (${N_NOVEL} transcritos novel MSTRG.*)
  - ${REF_CDS_FA}           CDS RefSeq filtrados (nucleotido)
  - ${REF_PEP_FA}           Proteinas RefSeq filtradas (amino, traduccion oficial)
  - ${NOVEL_FA}             mRNA spliced de novel (para TransDecoder en paso 09)

Siguiente paso:
  sbatch --dependency=afterok:\$SLURM_JOB_ID 09_transdecoder.slurm
EOF

echo "[DONE] Extraccion completa. Resumen en ${OUTDIR}/README.txt"
