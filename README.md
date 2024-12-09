# 1. Introduction
This script implements a two-step permutation test for cis protein quantitative trait loci (cis-pQTLs) analysis, described in Fang et al. (submitted). It requires plink v1.90 (https://www.cog-genomics.org/plink/). Please contact huatang@stanford.edu or hyfang@cnu.edu.cn for questions or bug report.

# 2. Usage
## 2.1 Prepare genotype data file, protein abundance file, covariate information file and protein meta information file

The genotype data file, protein abundance file and covariates information file should be in plink format. More details can be found in plink's manual and the example data.

The protein meta information file includes 6 columns without header and these 6 columns are id, chr, start, end, strand and non-missing count across samples for each gene. The column "strand" in protein meta information file is currently not be used. More details can be found in the example data.

## 2.2 Run two-step permutation test

Change variable names in Line 16-19 of "run_cispqtl_perm2step.sh" according to the genotype file, the protein abundance file, the covariate information file and the protein meta information file.

Use "source run_cispqtl_perm2step.sh;" in terminal for running the two-step permutation test.

## 2.3 Description of the results
The file "out_cispqtl.csv" under the folder "output/" is the final output that includes protein information (Columns 1-5), number of observations (Column 6), is or is not involved in 2 steps (Column 7), index SNP information (Column 8-14), permuation p-value (Column 15), q-value from permutation p-values (Column 16).

# 3. Example
The example input data is under the folder "data/". Use "source run_cispqtl_perm2step.sh;" in terminal and then the output file "out_cispqtl.csv" can be found under the folder "output/".


