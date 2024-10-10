#!/bin/bash
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

# Extract SR Run IDs (col 5) where Title (col 2) contains 'heart', omit header row
for sample in $(awk -F, 'NR >1 && $2 ~ /heart/ { print $5  }' SampleList.csv)
do
    echo "$(date) Getting FASTQ from: $sample"
    fasterq-dump $sample \
    --outdir $OUTDIR \
    --temp $TMP \
    --threads $SLURM_CPUS_PER_TASK \
    --progress \
    --split-files
done