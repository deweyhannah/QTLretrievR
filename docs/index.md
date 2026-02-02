# QTLretrievR

## Overview

QTLretrievR is an integration of **qtl2** and **intermediate** to:

- map genotypes to probes (genoprobs)
- calculate Quantitative Trait Loci (QTL) peaks
- identify mediators of QTL peaks
- identify phenotypic effects within experimental data

The package is designed for use with experimental data from **Diversity
Outbred (DO) mice**.

------------------------------------------------------------------------

## Getting started

### Installation

QTLretrievR is not available on CRAN. To install, please use the
instructions below:

    remotes::install_github("deweyhannah/QTLretrievR")

Please note that QTLretrievR has 2 Biocunductor packages as
dependencies, they may not install if `BiocManager` is not available.

### Submitting to a Job Scheduler

Currently, `QTLretrievR` is set up to identify available cores in a
linux, mac, or windows environment. It is also able to identify the
number of cores if you are submitting a job through a SLURM scheduler.
Just make sure to include both of the following lines in your header:

    #SBATCH --ntasks=N
    #SBATCH --cpus-per-task=X

When the number of cores is calculated for parallelization purposes, the
total number of cores that are identified are `N*X`. If you don’t
include both, and are submitting with a SLURM scheduler, it may not
identify the correct number of cores.

If you are using any other scheduler, make sure that you specify the
total number of cpus for the whole job. If applicable the number of
processes should be 1.

## Getting Help

If you find a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/deweyhannah/QTLretrievR/issues).
