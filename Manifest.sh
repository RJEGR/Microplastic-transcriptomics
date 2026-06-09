#!/usr/bin/env bash

for file in *_1.fq.gz; do
    [[ -e "$file" ]] || continue

    base="${file%_1.fq.gz}"

    # Split the base string into an array using '_' as the delimiter
    IFS='_' read -r -a parts <<< "$base"
    
    # parts[0] is "EP08", parts[1] is "2"
    bs1="${parts[0]}"
    bs2="${parts[0]}_${parts[1]}"

    # Output bs1 and bs2 as the first two columns
    printf '%s\t%s\t%s/%s_1.fq.gz\t%s/%s_2.fq.gz\n' "$bs1" "$bs2" "$PWD" "$base" "$PWD" "$base"
done > Manifest.tsv