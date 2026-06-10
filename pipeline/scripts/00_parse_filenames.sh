#!/usr/bin/env bash

# Declare associative arrays to group our data
declare -A grupo_map replica_map r1_map r2_map

# Read filenames line-by-line from STDIN
while IFS= read -r file || [[ -n "$file" ]]; do
    [[ -z "$file" ]] && continue

    # Resolve absolute path using readlink (fallback to PWD if run blindly on strings)
    abs_path=$(readlink -f "$file" 2>/dev/null || echo "$PWD/$file")

    # Isolate the basename
    basename="${file##*/}"

    # 1. Extract 'idx' using extended regex in sed
    idx=$(echo "$basename" | sed -E 's/(-[0-9]+)?_[A-Z0-9]+_L[0-9]+_[12]\.fq\.gz$//')

    # 2. Extract 'replica' (strip from _CKDL to the end)
    replica="${idx%%_CKDL*}"

    # 3. Extract 'grupo': strip _[digits] separator (Ch14_1 → Ch14);
    #    if no underscore separator exists, strip bare trailing digits (SES1 → SES)
    grupo=$(echo "$replica" | sed -E 's/_[0-9]+$//')
    [[ "${grupo}" == "${replica}" ]] && grupo=$(echo "$replica" | sed -E 's/[0-9]+$//')

    # 4. Map the reads based on filename endings
    if [[ "$basename" == *_1.fq.gz ]]; then
        r1_map["$idx"]="$abs_path"
        grupo_map["$idx"]="$grupo"
        replica_map["$idx"]="$replica"
    elif [[ "$basename" == *_2.fq.gz ]]; then
        r2_map["$idx"]="$abs_path"
    fi
done

# Print the TSV header
printf "%s\t%s\t%s\t%s\t%s\n" "idx" "grupo" "replica" "R1" "R2"

# Extract map keys, sort them, and print the tabular rows
while IFS= read -r idx; do
    [[ -z "$idx" ]] && continue
    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$idx" \
        "${grupo_map[$idx]}" \
        "${replica_map[$idx]}" \
        "${r1_map[$idx]}" \
        "${r2_map[$idx]}"
done < <(printf "%s\n" "${!grupo_map[@]}" | sort)