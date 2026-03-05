{\rtf1\ansi\ansicpg936\cocoartf2707
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/bin/bash\
set -euo pipefail\
\
READ_DIR=~/reads/\
\
HISAT2= PATH_TO_HISAT2\
INDEX=Dpse.fa\
\
GTF= Dpse.Convert.gtf\
FEATURECOUNTS=PATH_TO_FEATURECOUNTS\
\
THREADS=6\
MAPQ=30\
\
# =============== OUTPUT =================\
OUTDIR=RNAseq_results\
mkdir -p $\{OUTDIR\}/\{bam,sam,count\}\
\
############################################################\
# 1. Alignment\
############################################################\
\
for r1 in $\{READ_DIR\}/*_1.fq.gz; do\
    sample=$(basename $\{r1\} _1.fq.gz)\
    r2=$\{READ_DIR\}/$\{sample\}_2.fq.gz\
\
    echo "=== Processing $\{sample\} ==="\
\
    $\{HISAT2\} \\\
        -x $\{INDEX\} \\\
        -1 $\{r1\} \\\
        -2 $\{r2\} \\\
        --phred33 \\\
        -p $\{THREADS\} \\\
    | samtools view -q $\{MAPQ\} -b -@ $\{THREADS\} - \\\
    > $\{OUTDIR\}/bam/$\{sample\}.bam\
done\
\
############################################################\
# 2. Sort BAM\
############################################################\
\
for bam in $\{OUTDIR\}/bam/*.bam; do\
    sample=$(basename $\{bam\} .bam)\
    samtools sort $\{bam\} -@ $\{THREADS\} -O BAM \\\
        -o $\{OUTDIR\}/bam/$\{sample\}.sort.bam\
done\
\
############################################################\
# 3. Keep uniquely mapped reads\
############################################################\
\
for bam in $\{OUTDIR\}/bam/*.sort.bam; do\
    sample=$(basename $\{bam\} .sort.bam)\
    samtools view $\{bam\} \\\
        | grep 'NH:i:1' \\\
        | grep -v 'ZS:i' \\\
        > $\{OUTDIR\}/sam/$\{sample\}.unique.sam\
done\
\
############################################################\
# 4. featureCounts\
############################################################\
\
for sam in $\{OUTDIR\}/sam/*.sam; do\
    sample=$(basename $\{sam\} .unique.sam)\
\
    $\{FEATURECOUNTS\} \\\
        -T $\{THREADS\} \\\
        -p \\\
        -t exon \\\
        -g gene_id \\\
        -a $\{GTF\} \\\
        -o $\{OUTDIR\}/count/$\{sample\}.txt \\\
        $\{sam\}\
done\
\
############################################################\
# 5. Build count matrix\
############################################################\
\
cd $\{OUTDIR\}/count\
\
for f in *.txt; do\
    sample=$(basename $\{f\} .txt)\
    cut -f 7 $\{f\} \\\
        | grep -v '^#' \\\
        | grep -v '^Geneid' \\\
        | sed "1i$\{sample\}" \\\
        > $\{sample\}.cut\
done\
\
paste -d "\\t" *.cut > all_samples.count\
\
echo "Done."\
}