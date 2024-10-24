#!/bin/bash

# Directories
OUTDIR=/scratch/${USER}/data
mkdir -p $OUTDIR
$SAMPLE="BALF-C141" # testing with just 1 sample at the moment
FASTQDIR=$OUTDIR/$SAMPLE/

# Suffixes
OLD_R1="_1.fastq.gz"
OLD_R2="_2.fastq.gz"
NEW_R1="_S1_R1_001.fastq.gz"
NEW_R2="_S1_R2_001.fastq.gz"

# Install CellRanger
cd /scratch/${USER}/
curl -o cellranger-8.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.1.tar.gz?Expires=1729755313&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=i0mJn-cHGbBWyMrYPqvJcm-gxcAL6s0IzMaurS6uxVIFCNhppV3f68RlechIK7sa6NvYQLH8NdF8vYZJuNEGNUWHo2oujJubOnqixUA43vCpoTC94xGIeQ6mEs~z4HZN7YdYdwqtBFZH77bRY2BwFJGKHf1g6ZAEFW4vsmYfxwQXhMx7-8sygMefwYMnwTu8mH0aHoTTH8XUryisTZ4RIZAxVtyGsn-SqRVMb8mEnJPEMjnE7Gnf6oPBaahQW7oKV1L5v482M8dHFmUOnw1qOuUzllXuuKibSkKxMmxEMOwYP4uhgOunZWNB7JE90TJEJzI65Y-UGKJuqdnQLUVlOw__"
tar -xzvf cellranger-8.0.1.tar.gz

# Download and unzip 10x human reference dataset
cd $OUTDIR
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz

# Rename FASTQ to follow cellranger expected input
cd $FASTQDIR
for ff in *$OLD_R1
do
    baseName="${ff/$OLD_R1/}"
    mv -- "$baseName$OLD_R1" "$baseName$NEW_R1"
    mv -- "$baseName$OLD_R2" "$baseName$NEW_R2"
done
