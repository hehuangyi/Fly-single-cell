#!/bin/bash
set -e
set -o pipefail

# ===============================
# configuration
# ===============================
QUERY="query/Dpse.faa"
DBROOT="db"
OUTDIR="results"
EVALUE="1e-6"

mkdir -p ${OUTDIR}

# species list
SPECIES=(
Dmel
Dpse
Dsub
Dvir
Dwil
Dsim
Dere
Dsuz
Dana
Dmoj
Dgri
Uniprot
)

# ===============================
# function: make db if missing
# ===============================
make_db () {
    fasta=$1
    dbname=$2
    dbtype=$3

    if [ ! -f "${dbname}.pin" ] && [ ! -f "${dbname}.phr" ]; then
        echo "Building database: ${dbname}"
        makeblastdb -in ${fasta} -dbtype ${dbtype} -parse_seqids -out ${dbname}
    fi
}

# ===============================
# start
# ===============================
for sp in "${SPECIES[@]}"; do

    echo ""
    echo "=============================="
    echo "Processing ${sp}"
    echo "=============================="

    PEP=${DBROOT}/${sp}/${sp}.pep
    CDS=${DBROOT}/${sp}/${sp}.CDS.fa

    # output names
    blastp_out=${OUTDIR}/Dpse_to_${sp}.blastp.out
    blastn_out=${OUTDIR}/Dpse_to_${sp}.tblastn.out

    # -----------------
    # BLASTP
    # -----------------
    if [ -f ${PEP} ]; then
        make_db ${PEP} ${DBROOT}/${sp}/${sp}.pep prot

        blastp -query ${QUERY} \
               -db ${DBROOT}/${sp}/${sp}.pep \
               -outfmt 7 \
               -evalue ${EVALUE} \
               > ${blastp_out}

        if [ "${sp}" = "Dpse" ]; then
            # keep the 2nd hit (skip self)
            awk 'a[$1]++ == 1 {print}' ${blastp_out} > ${blastp_out/.out/.uniq}
        else
            # keep the best hit
            awk '!a[$1]++{print}' ${blastp_out} > ${blastp_out/.out/.uniq}
        fi

    else
        echo "Protein fasta missing for ${sp}"
    fi

    # -----------------
    # TBLASTN
    # -----------------
    if [ -f ${CDS} ]; then
        make_db ${CDS} ${DBROOT}/${sp}/${sp}.CDS nucl

        tblastn -query ${QUERY} \
                -db ${DBROOT}/${sp}/${sp}.CDS \
                -outfmt 7 \
                -evalue ${EVALUE} \
                > ${blastn_out}

        if [ "${sp}" = "Dpse" ]; then
            awk 'a[$1]++ == 1 {print}' ${blastn_out} > ${blastn_out/.out/.uniq}
        else
            awk '!a[$1]++{print}' ${blastn_out} > ${blastn_out/.out/.uniq}
        fi

    else
        echo "CDS fasta missing for ${sp}"
    fi

    echo "Finished ${sp}"
done

echo ""
echo "All species finished."
