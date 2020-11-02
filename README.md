# ssCTPR

## Description

ssCTPR is a polygenic risk score method based on elastic net regression with an additional cross-trait panelty, using summary statistics from trait of interest and trait(s) that are genetically correlated with the trait of interest, and accounting for Linkage Disequilibrium (LD) by a reference panel from target population.

## Installation

    install.packages("devtools")
    devtools::install_github("yingxi-kaylee/ssCTPR")
    
## References

Wonil Chung, Jun Chen, Constance Turman, Sara Lindstrom, Zhaozhong Zhu, Po-Ru Loh, Peter Kraft and Liming Liang (2019), Efficient cross-trait penalized regression increases prediction accuracy in large cohorts using secondary phenotypes. *Nature Communications*, 10(1), 569.

Mak, TSH, Porsch, RM, Choi, SW, Zhou, X, Sham, PC. Polygenic scores via penalized regression on summary statistics. *Genet. Epidemiol*. 2017; 41: 469– 480.
[DOI:10.1002/gepi.22050](https://doi.org/10.1002/gepi.22050).