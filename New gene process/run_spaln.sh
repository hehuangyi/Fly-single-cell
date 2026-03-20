#!/bin/bash
set -e
set -o pipefail

# ===============================
# Usage
# ===============================
if [ "$#" -lt 1 ]; then
    echo "Usage: bash run_spaln.sh <species_dir1> [species_dir2] ..."
    exit 1
fi

# ===============================
# Fixed input: Dpse protein fasta directory
# ===============================
F1="Dpse_denovo_faa"

# ===============================
# Read gene pairs
# ===============================
genes=($(cat denovolist.txt))

# Ensure the list contains an even number of entries (pairs)
if (( ${#genes[@]} % 2 != 0 )); then
    echo "Error: list size is not even. Expect gene pairs in denovolist.txt"
    exit 1
fi

# ===============================
# Loop over species provided by user
# ===============================
for F2 in "$@"; do

    # Remove .fa suffix to obtain species name
    SP=$(basename "$F2" .fa)

    echo ""
    echo "=============================="
    echo "Running species: ${SP}"
    echo "=============================="

    OUT="spaln_results_${SP}"
    mkdir -p "$OUT"

    for ((i=0; i<${#genes[@]}; i+=2)); do
        g1=${genes[$i]}
        g2=${genes[$i+1]}

        f1="${F1}/${g1}.fa"
        f2="${F2}/${g2}.fa"

        if [[ ! -f "$f1" ]]; then
            echo "Missing $f1"
            continue
        fi
        if [[ ! -f "$f2" ]]; then
            echo "Missing $f2"
            continue
        fi

        out_file="${OUT}/${g1}_vs_${g2}.spaln.out"

        echo "Running ${SP}: $g1 vs $g2 ..."
        spaln -O0 -H10 "$f2" "$f1" > "$out_file"
    done

    echo "Finished ${SP}"
done

echo ""
echo "All done."
