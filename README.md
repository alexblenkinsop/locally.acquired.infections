# Estimating the potential to prevent locally acquired HIV infections in a UNAIDS Fast-Track City, Amsterdam

This repository includes code and partial data for the analyses in [Blenkinsop, Monod, van Sighem et al. (2022)](https://doi.org/10.7554/eLife.76487) [![DOI](https://zenodo.org/badge/504520616.svg)](https://zenodo.org/badge/latestdoi/504520616)

## Data
The data folder contains subfolders with the following input files:
* trees
    * reconstructed phylogenetic trees labelled with risk group and year of sequence sample
* subgraphs
    * Amsterdam subgraphs extracted from trees
* subgraph_metatadata
    * classification of subgraphs as pre-existing by 2014 or emergent since 2014
* patient data
    * file containing a flag to indicate whether patient ID was virally suppressed by 2014
    * file containing a flag to indicate whether a patient had an estimated infection date after 2014
    * number of diagnosed individuals estimated to have been infected since 2014 by transmission risk group and place of birth
    * number of sequenced individuals estimated to have been infected since 2014 by transmission risk group, place of birth and HIV subtype
* infection times
    * estimated time-to-diagnosis data by risk group and migrant group
    * estimated infected individuals by year for MSM/non-MSM in Amsterdam from the European Centres for Disease Control (ECDC) HIV modelling tool


## Code
The analysis is run in 3 stages. The first stage requires sequence data and patient meta-data. The second two stages can be run using aggregated patient and phylogenetic data in the repository to replicate the results of the manuscript. 

### Phylogenetic analysis
1) scripts/pre-process-sequences.R - Pre-processes sequence data
2) scripts/phylo-analysis.R - Runs phylogenetic analysis

### Estimating the undiagnosed population
1) run-undiagnosed.R - Writes a shell script to estimate the proportion undiagnosed.

### Estimating locally acquired infections
1) submit-job-MSM.R & submit-job-HSX.R - Writes shell script to run the analysis on a computer cluster including preparing the stan data (stan-make-data.R), sampling using cmdstan, and post-processing. Requires the job name of the undiagnosed model as input.
2) scripts/post-processing-combine-MSM-HSX-results.R - Writes a shell script to combine the results of the independent MSM and heterosexual model. Both jobs must have completed before running.
