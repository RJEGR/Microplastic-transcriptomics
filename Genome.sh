#!/usr/bin/env bash

BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/853/765/GCF_963853765.1_xbMagGiga1.1"

echo "==> Descargando en paralelo usando xargs..."

# Listamos los archivos y los procesamos en paralelo (Máximo 3 a la vez)
printf "%s\n" \
    "GCF_963853765.1_xbMagGiga1.1_genomic.gtf.gz" \
    "GCF_963853765.1_xbMagGiga1.1_cds_from_genomic.fna.gz" \
    "GCF_963853765.1_xbMagGiga1.1_genomic.fna.gz" | \
xargs -I {} -P 3 curl -C - -O "${BASE_URL}/{}"

echo "==> Completado!"