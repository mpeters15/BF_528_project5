## Project Description

This repository contains the code for BF528 Project 5 - Individual Project.

This project revisits Project 3: Concordance of microarray and RNA-Seq differential gene expression. It seeks to reproduce a subset of analyses that were performed in Wang et al. 2014 (PMID: 4243706). The analyses done in this paper had several goals:
1. To characterize the concordance of differential gene expression across platforms
2. Test and compare how effective each platform is at detecting expected pathway-level effects based on each treatment’s mechanism of action
3. Assess the mechanism of action prediction accuracy of each platform using a test set.

The Data Curator and Programmer roles are performed in this project. NextFlow [DSL 2](https://nextflow.io/docs/latest/dsl2.html) is used to run the data processing and analyses from these roles.

## Contributors

 - Michael Peters

## Repository Contents

### How to use this repository:

1. Log into SCC and navigate to a working directory.
2. Run:
    ```
    module load nextflow
    nice nextflow run http://github.com/mpeters15/BF_528_project5 -r <latest commit hash>
    ```
3. Outputs will be stored in `nf_out`
