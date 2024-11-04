# Overview

**Project Title**: Bioinformatics Pipeline Recreation for BINF6430 Transcriptomics in Bioinformatics

**Contributors**: Dawn Mitchell, Courtney Wong, Jessica Kelly, Ethan Littlestone, Harshini Muthukumar

**Date**: Fall 2024

The purpose of this project is to identify a recent peer-reviewed journal article on an RNA-specific bioinformatics analysis technique, then explore, reproduce, verify, and discuss the presented pipeline methods.

# Selected Paper

We have chosen to perform a reproducibility study of the paper “propeller: testing for differences in cell type proportions in single cell data” (Phipson et al., 2022) published in the journal Bioinformatics. The main purpose of this study was to develop a robust method to utilize single-cell RNA-sequencing data for detecting significant differences in cell type composition between samples under different conditions. The challenge of using scRNA-seq data to research such effects is determining if observed differences in cell type proportions between samples of different conditions is biologically significant. This type of analysis is prone to false positive results due to technical bias and the inherent heterogeneous nature of cell type proportions between samples.

# Methods

The steps below, individually described in the following sections, are used to reproduce the analysis pipeline described in the paper.

1. `get_fastq.sh`: Obtains raw FASTQ zipped files for local use
2. `setup_cellranger.sh`: First-time setup for cellranger
3. `run_cellranger.sh`: Runs cellranger on FASTQ files in a directory

## get_fastq.sh

**Overview**: Obtains raw FASTQ zipped files for local use. Iterates through .csv containing samples and their links and feeds them into `curl` to download the file.

**Input**: A .csv file with 3 columns (no header): SRR, Sample Name, ENA link to fastq

**Output**: Downloaded FASTQ files directed into `/scratch/${USER}/data/XYZ/` where `XYZ` are the Sample Name identified in the .csv file. Each SRR sample will therefore become its own subdirectory and contain its associated FASTQ files.

**Example Usage**:

```
sbatch get_fastq.sh SampleList_covid.csv
```

**Example Input**: Note that these SRR accession IDs and associated FASTQ file link(s) can be easily obtained in bulk by searching for the project ID at https://sra-explorer.info/ . For example, the project ID of the COVID samples (bronchoalveolar immune cells) mentioned in the paper is SRP250732, as indicated at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145926 , and the first 3 paired samples are shown below.

```
SRR11181956,BALF-C143,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR111/056/SRR11181956/SRR11181956_1.fastq.gz
SRR11181956,BALF-C143,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR111/056/SRR11181956/SRR11181956_2.fastq.gz
SRR11181957,BALF-C144,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR111/057/SRR11181957/SRR11181957_1.fastq.gz
SRR11181957,BALF-C144,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR111/057/SRR11181957/SRR11181957_2.fastq.gz
SRR11537950,C149,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/050/SRR11537950/SRR11537950_1.fastq.gz
SRR11537950,C149,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/050/SRR11537950/SRR11537950_2.fastq.gz
...
```

**Example Output**: For the COVID samples, the following file structure will appear under `/scratch/${USER}/data/`

```
.
|-- BALF-C141
|   |-- SRR11181954_1.fastq.gz
|   `-- SRR11181954_1.fastq.gz
|-- BALF-C142
|   |-- SRR11181955_1.fastq.gz
|   `-- SRR11181955_2.fastq.gz
|-- BALF-C143
|   |-- SRR11181956_1.fastq.gz
|   `-- SRR11181956_2.fastq.gz
|-- BALF-C144
|   |-- SRR11181957_1.fastq.gz
|   `-- SRR11181957_2.fastq.gz
|-- BALF-C145
|   |-- SRR11181958_1.fastq.gz
|   `-- SRR11181958_2.fastq.gz
|-- BALF-C146
|   |-- SRR11181959_1.fastq.gz
|   `-- SRR11181959_2.fastq.gz
|-- C100
|   |-- SRR11537948_1.fastq.gz
|   `-- SRR11537948_2.fastq.gz
|-- C148
|   |-- SRR11537949_1.fastq.gz
|   `-- SRR11537949_2.fastq.gz
|-- C149
|   |-- SRR11537950_1.fastq.gz
|   `-- SRR11537950_2.fastq.gz
|-- C152
|   |-- SRR11537951_1.fastq.gz
|   `-- SRR11537951_2.fastq.gz
|-- C51
|   |-- SRR11537946_1.fastq.gz
|   `-- SRR11537946_2.fastq.gz
`-- C52
    |-- SRR11537947_1.fastq.gz
    `-- SRR11537947_2.fastq.gz
```

Note that Read 1 is the 16nt barcode and 10nt UMI:

```
(base) [wong.co@login-00 data]$ cd BALF-C141
(base) [wong.co@login-00 BALF-C141]$ zcat SRR11181954_S1_R1_001.fastq.gz | head
@SRR11181954.1 /1
GAACATCGTAACTGTATATTATTCTC
+
DCCC@;8C@7936;<->440.):8->
@SRR11181954.2 /1
AAGCCGCAGTGCGATGGTACACAGCG
+
B:DCDDB?CBC=<::BA6<>.@.9(A
@SRR11181954.3 /1
ACTGCTCCAGCCTGTGCCTGATGGTG
```

And Read 2 is the actual transcript:

```
(base) [wong.co@login-00 BALF-C141]$ zcat SRR11181954_S1_R2_001.fastq.gz | head
@SRR11181954.1 /2
GTCTGGCCAGCTGGTGAACTGAATGTGAGTCACCTCTCTTCCAGTTGCTTTTTCTTTTTTATTTACAATGTTCAATTTCTGAATGATGTAAGCTGGACAT
+
BCBBB6B4BCBBBB;BBAA@CBCB:ACD@>9D95C5?BC@@A@B:CB?A6::B@?>?@@@B<;BA=BA7E3?ACB7AA@B>A@9?<CB+@C1AACC(:A>
@SRR11181954.2 /2
ATTGCGCCAGGTTTCAATTTCTATCGCCTATACTTTATTTGGGTAAATGGTTTGGCTAAGGTTGTCTGGTAGTAAGGTGGAGTGGGTTTGGGGCTAGGCT
+
<A@DDCCDCCB1@A<AD@=BBCCCAC;DCADCCCCBCDBDDDC=CACCACC6CD@BDACEC:BEAC>DB@CEBACCDDDCEE=CDC@4>CDDC(DCCD%=
@SRR11181954.3 /2
CCTTTCCTGTTCACTCTACCCTTTGACTCTAAATCTCAAAGCCAGTGTTGGGGCCCAGTGGCTCCATTCGATTGAAACATGGCCAATGATATCCAAGAGC
```

## setup_cellranger.sh

**Purpose**: First-time setup for cellranger.

**Input**: None *(in work)*

**Output**: Downloads and unzips CellRanger from 10X Genomics into `/scratch/${USER}/cellranger-8.0.1`. Downloads and unzips the 10X Genomics GRCh38 human reference dataset into `/scratch/${USER}/data/refdata-gex-GRCh38-2024-A`. Renames all paired FASTQ files within a hard-coded directory *(in work)* to conform with cellranger's expected format: `[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz` or `[Sample Name]_S1_[Read Type]_001.fastq.gz` per https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-specifying-fastqs

## run_cellranger.sh

**Purpose**: Runs cellranger in `count` mode to perform sample de-multiplexing, barcode processing, and single-cell 5’ UMI counting with human GRCh38 as the reference genome.

**Input**: FASTQs for the sample must be located under `/scratch/${USER}/data/<SampleName>/`

**Output**: Creates default CellRanger outputs inside the same folder as the FASTQs: `/scratch/${USER}/data/<SampleName>/<SampleName>/`. Key files are moved to `/home/${USER}/propeller`.

**Example Usage**: The script is currently hardcoded to create counts for the FASTQs located under `/scratch/${USER}/data/BALF-C141` only but may be modified to loop over multiple samples.

**Example Output**: The key outputs reside in the subfolder `outs/`, with the following file structure. The files we are primarily interested in are the summary report `web_summary.html` and the feature barcode matrix `filtered_feature_bc_matrix.h5`.

```
.
|-- analysis
|   |-- clustering
|   |-- diffexp
|   |-- pca
|   |-- tsne
|   `-- umap
|-- cloupe.cloupe
|-- filtered_feature_bc_matrix
|   |-- barcodes.tsv.gz
|   |-- features.tsv.gz
|   `-- matrix.mtx.gz
|-- filtered_feature_bc_matrix.h5
|-- metrics_summary.csv
|-- molecule_info.h5
|-- possorted_genome_bam.bam
|-- possorted_genome_bam.bam.bai
|-- raw_feature_bc_matrix
|   |-- barcodes.tsv.gz
|   |-- features.tsv.gz
|   `-- matrix.mtx.gz
|-- raw_feature_bc_matrix.h5
`-- web_summary.html
```

# BALF_C141_seurat.Rmd

**Purpose**: Conduct post-processing analysis on the filtered feature barcode matrix from a single selected sample, BALF-C141, using the R package `Seurat`.

# References

Liao, M., Liu, Y., Yuan, J. et al. Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19. Nat Med 26, 842–844 (2020). https://doi.org/10.1038/s41591-020-0901-9

Phipson, B., Sim, C. B., Porrello, E. R., Hewitt, A. W., Powell, J., & Oshlack, A. (2022). propeller: testing for differences in cell type proportions in single cell data. Bioinformatics, 38(20), 4720–4726. https://doi.org/10.1093/bioinformatics/btac582

Zheng, G. X. Y. et al. (2017). Massively parallel digital transcriptional profiling of single cells. Nature Communications 8: 1-12, doi:10.1038/ncomms14049
