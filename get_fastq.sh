#!/bin/bash

# Usage: sbatch get_fastq.sh <SampleList>
    # <SampleList> is a csv with 3 columns (no header): SRR, Sample Name, ENA link to fastq

#SBATCH --job-name=get_fastq
#SBATCH --partition=short
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem 64G
#SBATCH -t 8:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --out=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.err

# Directories
OUTDIR=/scratch/${USER}/data
mkdir -p $OUTDIR

SAMPLES=$(realpath "$1")

# Iterate through each SRR of each sample
for sample in $(awk -F, '{ print $2 }' "$SAMPLES" | uniq) # Fetal1, Fetal2, ..., Adult2, Adult3
do
    echo "$(date) Processing sample: $sample"
    mkdir -p $OUTDIR/$sample
    pushd $OUTDIR/$sample
    awk -v sname="$sample" -F, '$2 == sname { print $3 }' "$SAMPLES" | \
    xargs curl --parallel -O
    popd
done