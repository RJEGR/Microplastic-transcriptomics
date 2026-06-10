#!/usr/bin/env bash
# Genera Manifest.tsv con 5 columnas: idx, grupo, replica, R1, R2
# Patron de nombre: {grupo}_{rep}_{CKDL...}-..._L{n}_1.fq.gz

printf 'idx\tgrupo\treplica\tR1\tR2\n' > Manifest.tsv

for file in *_1.fq.gz; do
    [[ -e "$file" ]] || continue

    base="${file%_1.fq.gz}"

    IFS='_' read -r -a parts <<< "$base"

    grupo="${parts[0]}"
    replica="${parts[0]}_${parts[1]}"
    idx="${parts[0]}_${parts[1]}_${parts[2]%%-*}"

    printf '%s\t%s\t%s\t%s/%s_1.fq.gz\t%s/%s_2.fq.gz\n' \
        "$idx" "$grupo" "$replica" "$PWD" "$base" "$PWD" "$base"
done >> Manifest.tsv