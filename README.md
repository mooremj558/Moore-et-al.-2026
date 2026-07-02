# Moore-et-al.-2026
Code for manuscript Moore et al. 2026

## Compute 
Two compute systems were used in this study.
The WashU HTCF and the Wang Lab systems. Both are slurm based system.

## Workflows/Code
Handling of ATAC, WGBS, and RNA-seq data was described previously (https://github.com/twlab/epitherapy_induced_antigen_GBM)

Scripts for processing of raw Hi-C and CTCF CUT&RUN files are in the "Data Processing" directory.
Scripts for downstream data processing (i.e. TAD/loop calling for Hi-C) are included in the "Analysis" directory in separate sub-directories.

R and Jupyter notebooks detailing all data analysis for each figure are included in the "Analysis" directory.
In the "Analysis" directory is also a "paper_final_master.sh" file which contains sections for each analysis beginning with running preprocessing scripts as well as commands for manipulating input files for use with R/Jupyter notebook code.
