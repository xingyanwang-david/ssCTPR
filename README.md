# ssCTPR

## Description

ssCTPR is a Polygenic Risk Score (PRS) method based on elastic net regression with an additional cross-trait penalty, using summary statistics from trait of interest and trait(s) that are genetically correlated with the trait of interest while accounting for Linkage Disequilibrium (LD) via a reference panel from the target population.

## Installation

    install.packages("devtools")
    devtools::install_github("yingxi-kaylee/ssCTPR")
    
## Reference

Chung, W., Chen, J., Turman, C., Lindstrom, S., Zhu, Z., Loh, P. R., . . . Liang, L. (2019). Efficient cross-trait penalized regression increases prediction accuracy in large cohorts using secondary phenotypes. *Nat Commun*, 10(1), 569. [doi:10.1038/s41467-019-08535-0](https://doi:10.1038/s41467-019-08535-0)

Mak, T. S. H., Porsch, R. M., Choi, S. W., Zhou, X., & Sham, P. C. (2017). Polygenic scores via penalized regression on summary statistics. *Genet Epidemiol*, 41(6), 469-480. [doi:10.1002/gepi.22050](https://doi.org/10.1002/gepi.22050)