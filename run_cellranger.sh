#!/bin/bash

#SBATCH --job-name=run_cellranger
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
SAMPLE="BALF-C141" # testing with just 1 sample at the moment
FASTQDIR=$OUTDIR/$SAMPLE/
REF=$OUTDIR/refdata-gex-GRCh38-2024-A

# Add CellRanger to path
export PATH=/scratch/${USER}/cellranger-8.0.1:$PATH

cd $FASTQDIR

# Run CellRanger
cellranger count \
--id=$SAMPLE \
--transcriptome=$REF \
--fastqs=$FASTQDIR \
--create-bam=true \
--localcores=$SLURM_CPUS_PER_TASK \
--localmem=$SLURM_MEM_PER_NODE

# Copy key outputs
mkdir /home/${USER}/propeller
cp $SAMPLE/outs/filtered_feature_bc_matrix.h5 /home/${USER}/propeller/
cp -r $SAMPLE/outs/filtered_feature_bc_matrix /home/${USER}/propeller/
cp $SAMPLE/outs/web_summary.html /home/${USER}/propeller/