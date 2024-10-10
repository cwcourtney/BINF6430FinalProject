#!/bin/bash

# Usage: sbatch get_fastq <SampleList>
    # <SampleList> is a csv with two columns (no header): SRR, Sample Name

#SBATCH --job-name=get_fastq
#SBATCH --partition=express
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem 32G
#SBATCH -t 1:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --out=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.err

# Modules
module load sratoolkit

# Directories
OUTDIR=/scratch/${USER}/data
mkdir -p $OUTDIR
TMP=/scratch/${USER}/tmp
mkdir -p $TMP

# Iterate through each SRR of each sample
for sample in $(awk -F, '{ print $2 }' $1 | uniq) # Fetal1, Fetal2, ..., Adult2, Adult3
do
    echo "$(date) Processing sample: $sample"
    mkdir -p $OUTDIR/$sample
    for srr in $(awk -v sname="$sample" -F, '$2 == sname { print $1 }' $1) # SRR...
    do
        echo "$(date) Getting FASTQ for run: $srr"
        fasterq-dump $srr \
        --outdir $OUTDIR/$sample \
        --temp $TMP \
        --threads $SLURM_CPUS_PER_TASK \
        --split-files \
        --details
    done
done