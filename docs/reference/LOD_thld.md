# Determine significant and suggestive LOD thresholds

Determine significant and suggestive LOD thresholds for your QTLs, using
`scan1perm` from `qtl2`. Select the number of genes to permute, the
number of permutations, and the number of genes to run in each batch.

## Usage

``` r
LOD_thld(
  mapping,
  tissue,
  annots = NULL,
  n.gene = 75,
  n.perm = 750,
  batch.size = 5
)
```

## Arguments

- mapping:

  Mapping list from `mapQTL`, or full path to `.rds` containing one.

- tissue:

  Tissue to determine thresholds for.

- annots:

  Data frame of phenotype annotations for filtering to autosomal
  phenotypes. Columns must include "id", "symbol", "start", "end".

- n.gene:

  Number of phenotypes to run permutations on. Default is 75.

- n.perm:

  Number of permutations to run per phenotype. Default is 750.

- batch.size:

  Number of genes in each parallelized batch.

## Value

A list containing:

- perms:

  A `scan1perm` object of permutations for all selected genes.

- thresholds:

  A list with two elements:

  significant

  :   Median LOD threshold at alpha = 0.05.

  suggestive

  :   Median LOD threshold at alpha = 0.4.
