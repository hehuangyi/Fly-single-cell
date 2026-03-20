# ===============================
# Step1Cactus alignment
# ===============================
cactus jobStore seqFile.txt output.hal --maxCores 8

# ===============================
# Step2  preprocess
# ===============================
sh getbed.sh
awk 'BEGIN{OFS="\t"} {$5=($5=="." ? 0 : $5); print}' denovo_genes.bed > genes.fixed.bed
sh liftover.sh
python mergebed.py
python batch_extract_lifted.py

# ===============================
# Step3  Run genewise for each species
# ===============================
species=(Dmel Dpse Dsub Dvir Dwil Dsim Dere Dsuz Dana Dmoj Dgri)
for sp in "${species[@]}"; do
    echo "Running genewise for $sp ..."

    target="./${sp}_lifted_sequences.fa"

    if [[ ! -f "$target" ]]; then
        echo "  Missing $target , skip."
        continue
    fi

    genewise "$REF" "$target" -sum -gff > "${sp}_genewise.gff"
    genewise "$REF" "$target" -sum -pep > "${sp}_genewise.pep"

done

echo "All genewise jobs finished."

cut -f4 *.lifted.bed | sort | uniq > merged_uniq_gene.txt

# ===============================
# Step4  spaln
# ===============================

awk '/^>/ {f=substr($0,2); gsub(/[^A-Za-z0-9_.-]/,"",f); print > f".fa"; next} {print >> f".fa"}' ../Dpse_denovo.faa
species=(Dmel Dsub Dvir Dwil Dsim Dere Dsuz Dana Dmoj Dgri)


for SP in "${species[@]}"; do
   bash run_spaln.sh ${SP}
    cd spaln_results_${SP}
    cat *.out > all.gff
    sed 's/Target=[^; ]*//g' all.gff > all.clean.gff
    gffread all.clean.gff \
        -g ../../02.denovo_v2/${SP}_lifted_sequences.fa \
        -y all.pep.fa

    cd ..

     python merge_and_mafft.py ${SP}
     
    python batch_denovo.py merged_fa_${SP} denovo_results_${SP}

    echo "Finished ${SP}"
done
