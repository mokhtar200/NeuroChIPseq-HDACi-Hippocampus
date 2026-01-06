#!/bin/bash

# ---------------------------
# 0️⃣ Set directories
# ---------------------------
BASE_DIR=~/Desktop/histone
GENOME=$BASE_DIR/mm10.fa
THREADS=4

# Bowtie2 index (if not built yet)
if [ ! -f "${BASE_DIR}/mm10.1.bt2" ]; then
    echo "Building Bowtie2 index..."
    bowtie2-build $GENOME $BASE_DIR/mm10
fi

# ---------------------------
# 1️⃣ Define histone marks and samples
# ---------------------------
declare -A MARKS
# key = histone mark, value = FASTQ files (control first, treatment second)
MARKS["AcH2B"]="SRR647745.fastq SRR647746.fastq"
MARKS["AcH3K9,14"]="SRR647752.fastq SRR647753.fastq"
MARKS["AcH4K12"]="SRR647755.fastq SRR647756.fastq"

# ---------------------------
# 2️⃣ Loop over each histone mark
# ---------------------------
for mark in "${!MARKS[@]}"; do
    echo "================== Processing $mark =================="
    cd $BASE_DIR/$mark

    # Create output folders
    mkdir -p bam
    mkdir -p peaks

    # Split FASTQ files (control treatment)
    readarray -t fastqs <<< "${MARKS[$mark]}"

    CONTROL_FASTQ=${fastqs[0]}
    TREAT_FASTQ=${fastqs[1]}

    # 2a️⃣ Alignment with Bowtie2
    echo "Aligning $CONTROL_FASTQ and $TREAT_FASTQ..."
    bowtie2 -x $BASE_DIR/mm10 -U $CONTROL_FASTQ -S ${CONTROL_FASTQ%.fastq}.sam --very-sensitive -p $THREADS
    bowtie2 -x $BASE_DIR/mm10 -U $TREAT_FASTQ -S ${TREAT_FASTQ%.fastq}.sam --very-sensitive -p $THREADS

    # 2b️⃣ Convert SAM → sorted BAM
    for f in ${CONTROL_FASTQ} ${TREAT_FASTQ}; do
        samtools view -bS ${f%.fastq}.sam | samtools sort -o bam/${f%.fastq}_sorted.bam
        samtools index bam/${f%.fastq}_sorted.bam
    done

    # Remove SAM to save space
    rm *.sam

    # 2c️⃣ Remove PCR duplicates
    for f in ${CONTROL_FASTQ} ${TREAT_FASTQ}; do
        samtools markdup -r bam/${f%.fastq}_sorted.bam bam/${f%.fastq}_rmdup.bam
        samtools index bam/${f%.fastq}_rmdup.bam
    done

    # ---------------------------
    # 3️⃣ MACS2 Peak Calling
    # ---------------------------
    CONTROL_BAM=bam/${CONTROL_FASTQ%.fastq}_rmdup.bam
    TREAT_BAM=bam/${TREAT_FASTQ%.fastq}_rmdup.bam
    OUT_NAME=${mark}_TSA_vs_Veh

    echo "Running MACS2 peak calling..."
    macs2 callpeak \
        -t $TREAT_BAM \
        -c $CONTROL_BAM \
        -f BAM \
        -g mm \
        -n $OUT_NAME \
        --broad \
        --broad-cutoff 0.1 \
        --outdir peaks

    echo "Finished $mark"
    echo "------------------------------------------------------"
done

echo "================== ALL MARKS DONE ===================="
