# Overview

**Project Title**: Bioinformatics Pipeline Recreation for BINF6430 Transcriptomics in Bioinformatics

**Contributors**: Dawn Mitchell, Courtney Wong, Jessica Kelly, Ethan Littlestone, Harshini Muthukumar

**Date**: Fall 2024

The purpose of this project is to identify a recent peer-reviewed journal article on an RNA-specific bioinformatics analysis technique, then explore, reproduce, verify, and discuss the presented pipeline methods.

# Selected Paper

We have chosen to perform a reproducibility study of the paper “propeller: testing for differences in cell type proportions in single cell data” (Phipson et al., 2022) published in the journal Bioinformatics. The main purpose of this study was to develop a robust method to utilize single-cell RNA-sequencing data for detecting significant differences in cell type composition between samples under different conditions. The challenge of using scRNA-seq data to research such effects is determining if observed differences in cell type proportions between samples of different conditions is biologically significant. This type of analysis is prone to false positive results due to technical bias and the inherent heterogeneous nature of cell type proportions between samples.

# Methods

The steps below, individually described in the following sections, are used to reproduce the analysis pipeline described in the paper.

1. `get_fastq.sh`: Obtain raw FASTQ zipped files for local use
2. `align_reads.sh`: Align reads to a reference genome
3. ...

## get_fastq.sh

...

## align_reads.sh

...

# References

Phipson, B., Sim, C. B., Porrello, E. R., Hewitt, A. W., Powell, J., & Oshlack, A. (2022). propeller: testing for differences in cell type proportions in single cell data. Bioinformatics, 38(20), 4720–4726. https://doi.org/10.1093/bioinformatics/btac582
