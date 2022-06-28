# Estimating the potential to prevent locally acquired HIV infections in a UNAIDS Fast-Track City, Amsterdam

This repository includes code and partial data for the analyses.

## Data
The data folder contains
1)  Estimated time-to-diagnosis data by risk group and migrant group
2)  Reconstructed phylogenetic trees labelled with risk group and year of sequence sample

## Code
The analysis is run in 3 stages:

### Phylogenetic analysis
1) scripts/pre-process-sequences.R - Pre-processes the sequence data
2) scripts/phylo-analysis.R - Runs phylogenetic analysis

### Estimating the undiagnosed population
1) scripts/undiagnosed_bytrmgroup.R - Runs the stan model to estimate the proportion undiagnosed


### Estimating locally acquired infections
1) submit-job.R - Writes a shell script to run the analysis on a computer cluster including preparing the stan data (stan-make-data.R), sampling using cmdstan, and post-processing

[![DOI](https://zenodo.org/badge/504520616.svg)](https://zenodo.org/badge/latestdoi/504520616)
