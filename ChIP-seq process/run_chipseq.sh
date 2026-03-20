#!/bin/bash
set -e
set -o pipefail

THREADS=8
INDEX=Dpse_index
GTF=Dpse_sorted.gtf
FASTQ_DIR=../00.fastq

BAM_DIR=bam
BW_DIR=bw
SIG_DIR=signal
LOG_DIR=logs

mkdir -p ${BAM_DIR} ${BW_DIR} ${SIG_DIR} ${LOG_DIR}

# =====================================
# sample table
# format:
# SAMPLE   IP_R1 IP_R2 INPUT_R1 INPUT_R2
# =====================================

SAMPLES=(
"H3K9me3 SRR24006472_1.fastq SRR24006472_2.fastq SRR24006438_1.fastq SRR24006438_2.fastq"
)

# =====================================
# function
# =====================================
process_sample () {

    SAMPLE=$1
    IP_R1=$2
    IP_R2=$3
    IN_R1=$4
    IN_R2=$5

    echo ""
    echo "=============================="
    echo "Processing ${SAMPLE}"
    echo "=============================="

    IP=${SAMPLE}_IP
    INPUT=${SAMPLE}_Input

    # -----------------
    # alignment
    # -----------------
    bowtie2 -p ${THREADS} -x ${INDEX} \
        -1 ${FASTQ_DIR}/${IP_R1} \
        -2 ${FASTQ_DIR}/${IP_R2} \
        -S ${BAM_DIR}/${IP}.sam \
        > ${LOG_DIR}/${IP}.bowtie2.log 2>&1

    bowtie2 -p ${THREADS} -x ${INDEX} \
        -1 ${FASTQ_DIR}/${IN_R1} \
        -2 ${FASTQ_DIR}/${IN_R2} \
        -S ${BAM_DIR}/${INPUT}.sam \
        > ${LOG_DIR}/${INPUT}.bowtie2.log 2>&1

    # -----------------
    # sam -> bam
    # -----------------
    samtools view -bS ${BAM_DIR}/${IP}.sam > ${BAM_DIR}/${IP}.bam
    samtools view -bS ${BAM_DIR}/${INPUT}.sam > ${BAM_DIR}/${INPUT}.bam

    # -----------------
    # fixmate
    # -----------------
    samtools sort -n -o ${BAM_DIR}/${IP}.namesort.bam ${BAM_DIR}/${IP}.bam
    samtools sort -n -o ${BAM_DIR}/${INPUT}.namesort.bam ${BAM_DIR}/${INPUT}.bam

    samtools fixmate -m ${BAM_DIR}/${IP}.namesort.bam ${BAM_DIR}/${IP}.fix.bam
    samtools fixmate -m ${BAM_DIR}/${INPUT}.namesort.bam ${BAM_DIR}/${INPUT}.fix.bam

    # -----------------
    # coordinate sort
    # -----------------
    samtools sort -o ${BAM_DIR}/${IP}.sort.bam ${BAM_DIR}/${IP}.fix.bam
    samtools sort -o ${BAM_DIR}/${INPUT}.sort.bam ${BAM_DIR}/${INPUT}.fix.bam

    # -----------------
    # remove duplicates
    # -----------------
    samtools markdup -r ${BAM_DIR}/${IP}.sort.bam ${BAM_DIR}/${IP}.rmdup.bam
    samtools markdup -r ${BAM_DIR}/${INPUT}.sort.bam ${BAM_DIR}/${INPUT}.rmdup.bam

    samtools index ${BAM_DIR}/${IP}.rmdup.bam
    samtools index ${BAM_DIR}/${INPUT}.rmdup.bam

    # -----------------
    # bigwig
    # -----------------
    bamCompare \
        -b1 ${BAM_DIR}/${IP}.rmdup.bam \
        -b2 ${BAM_DIR}/${INPUT}.rmdup.bam \
        --operation ratio \
        --scaleFactorsMethod SES \
        -bs 1 \
        -p ${THREADS} \
        -o ${BW_DIR}/${SAMPLE}.ratio.bw

    # -----------------
    # gene signal
    # -----------------
    bigWigToWig ${BW_DIR}/${SAMPLE}.ratio.bw ${SIG_DIR}/${SAMPLE}.ratio.wig

    cat ${SIG_DIR}/${SAMPLE}.ratio.wig | \
    bedtools map -c 4 -o mean -b - -a ${GTF} \
    > ${SIG_DIR}/${SAMPLE}.gene_mean.txt

    echo "Finished ${SAMPLE}"
}

# =====================================
# run all samples
# =====================================
for entry in "${SAMPLES[@]}"; do
    process_sample $entry
done

echo ""
echo "All samples completed."
