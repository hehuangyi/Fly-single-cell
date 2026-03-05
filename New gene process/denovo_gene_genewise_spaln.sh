{\rtf1\ansi\ansicpg936\cocoartf2707
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh15760\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # ===============================\
# Step1Cactus alignment\
# ===============================\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 cactus jobStore seqFile.txt output.hal --maxCores 8\
\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 # ===============================\
# Step2  preprocess\
# ===============================\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 sh getbed.sh\
awk 'BEGIN\{OFS="\\t"\} \{$5=($5=="." ? 0 : $5); print\}' denovo_genes.bed > genes.fixed.bed\
sh liftover.sh\
python mergebed.py\
python batch_extract_lifted.py\
\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 # ===============================\
# Step3  Run genewise for each species\
# ===============================\
species=(Dmel Dpse Dsub Dvir Dwil Dsim Dere Dsuz Dana Dmoj Dgri)\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 for sp in "$\{species[@]\}"; do\
    echo "Running genewise for $sp ..."\
\
    target="./$\{sp\}_lifted_sequences.fa"\
\
    if [[ ! -f "$target" ]]; then\
        echo "  Missing $target , skip."\
        continue\
    fi\
\
    genewise "$REF" "$target" -sum -gff > "$\{sp\}_genewise.gff"\
    genewise "$REF" "$target" -sum -pep > "$\{sp\}_genewise.pep"\
\
done\
\
echo "All genewise jobs finished."\
\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 cut -f4 *.lifted.bed | sort | uniq > merged_uniq_gene.txt\
\
# ===============================\
# Step4  spaln\
# ===============================\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 \
awk '/^>/ \{f=substr($0,2); gsub(/[^A-Za-z0-9_.-]/,"",f); print > f".fa"; next\} \{print >> f".fa"\}' ../Dpse_denovo.faa\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 species=(Dmel Dsub Dvir Dwil Dsim Dere Dsuz Dana Dmoj Dgri)\
\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 \
for SP in "$\{species[@]\}"; do\
   bash run_spaln.sh $\{SP\}\
    cd spaln_results_$\{SP\}\
    cat *.out > all.gff\
    sed 's/Target=[^; ]*//g' all.gff > all.clean.gff\
    gffread all.clean.gff \\\
        -g ../../02.denovo_v2/$\{SP\}_lifted_sequences.fa \\\
        -y all.pep.fa\
\
    cd ..\
\
     python merge_and_mafft.py $\{SP\}\
     \
    python batch_denovo.py merged_fa_$\{SP\} denovo_results_$\{SP\}\
\
    echo "Finished $\{SP\}"\
done\
}