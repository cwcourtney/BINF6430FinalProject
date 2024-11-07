# Overview

**Project Title**: Bioinformatics Pipeline Recreation for BINF6430 Transcriptomics in Bioinformatics

**Contributors**: Dawn Mitchell, Courtney Wong, Jessica Kelly, Ethan Littlestone, Harshini Muthukumar

**Date**: Fall 2024

The purpose of this project is to identify a recent peer-reviewed journal article on an RNA-specific bioinformatics analysis technique, then explore, reproduce, verify, and discuss the presented pipeline methods.

# Selected Paper

We have chosen to perform a reproducibility study of the paper “propeller: testing for differences in cell type proportions in single cell data” (Phipson et al., 2022) published in the journal Bioinformatics. The main purpose of this study was to develop a robust method to utilize single-cell RNA-sequencing data for detecting significant differences in cell type composition between samples under different conditions. The challenge of using scRNA-seq data to research such effects is determining if observed differences in cell type proportions between samples of different conditions is biologically significant. This type of analysis is prone to false positive results due to technical bias and the inherent heterogeneous nature of cell type proportions between samples.

# Methods

## scRNA-seq Workflow Exploration on COVID lung samples

See `scripts/scRNAseq_workflow.html`

# References

Liao, M., Liu, Y., Yuan, J. et al. Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19. Nat Med 26, 842–844 (2020). https://doi.org/10.1038/s41591-020-0901-9

Phipson, B., Sim, C. B., Porrello, E. R., Hewitt, A. W., Powell, J., & Oshlack, A. (2022). propeller: testing for differences in cell type proportions in single cell data. Bioinformatics, 38(20), 4720–4726. https://doi.org/10.1093/bioinformatics/btac582

Zheng, G. X. Y. et al. (2017). Massively parallel digital transcriptional profiling of single cells. Nature Communications 8: 1-12, doi:10.1038/ncomms14049
