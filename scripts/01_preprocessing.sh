#!/bin/bash

# Set paths
DATA_DIR=../data
RAW_DIR=$DATA_DIR/raw
FASTQ_DIR=$DATA_DIR/fastq
BAM_DIR=$DATA_DIR/bam
PEAK_DIR=$DATA_DIR/peaks

mkdir -p $RAW_DIR $FASTQ_DIR $BAM_DIR $PEAK_DIR

# ------------------------------
# 1. Download SRR files
# ------------------------------
# TSA and Vehicle samples for AcH2B, AcH3K9,14, AcH4K12
SRRS=("SRR647745" "SRR647746" "SRR647747" "SRR647752" "SRR647753" "SRR647755" "SRR647756")

for SRR in "${SRRS[@]}"
do
  echo "Downloading $SRR..."
  prefetch $SRR -O $RAW_DIR
  fastq-dump --split-3 --outdir $FASTQ_DIR $RAW_DIR/$SRR.sra
done

# ------------------------------
# 2. Align reads to mm10 genome
# ------------------------------
for FASTQ in $FASTQ_DIR/*.fastq
do
  SAMPLE=$(basename $FASTQ .fastq)
  echo "Aligning $SAMPLE..."
  bowtie2 -x /path/to/mm10_index -U $FASTQ -S $BAM_DIR/$SAMPLE.sam
  samtools view -bS $BAM_DIR/$SAMPLE.sam > $BAM_DIR/$SAMPLE.bam
  samtools sort $BAM_DIR/$SAMPLE.bam -o $BAM_DIR/${SAMPLE}_sorted.bam
  samtools index $BAM_DIR/${SAMPLE}_sorted.bam
done

# ------------------------------
# 3. Peak calling using MACS2
# ------------------------------
# Example: AcH2B TSA vs Vehicle
macs2 callpeak \
  -t $BAM_DIR/SRR647747_sorted.bam \
  -c $BAM_DIR/SRR647745_sorted.bam,$BAM_DIR/SRR647746_sorted.bam \
  -f BAM \
  -g mm \
  -n AcH2B_TSA_vs_Veh \
  --outdir $PEAK_DIR \
  --broad \
  --broad-cutoff 0.1

# Repeat MACS2 calls for AcH3K9,14 and AcH4K12 with their TSA vs Vehicle BAMs
