# Prepare and run mediation for a set of QTL peaks

Prepare and run mediation for a set of QTL peaks

## Usage

``` r
modiFinder(
  peaks,
  mapping,
  sigLOD = 7.5,
  annots,
  outdir = NULL,
  med_out = "mediation.rds",
  total_cores = NULL,
  save = "sr",
  distOnly = TRUE,
  hsOnly = FALSE
)
```

## Arguments

- peaks:

  List of dataframes containing QTL peaks for each tissue, or full path
  to `.rds` containing one.

- mapping:

  Mapping list from `mapQTL`, or full path to `.rds` containing one.

- sigLOD:

  Significant LOD threshold to use to filter phenotypes for mediation.
  Default is 7.5

- annots:

  Annotations file. Contains mapping information for phenotypes.
  Dataframe, or tsv. Columns must include "id", "symbol", "start",
  "end".

- outdir:

  Directory to save output files. Default is `NULL`.

- med_out:

  String indicating the name of the output file containing mediation
  results for mediation within a phenotype. This file will be saved in
  `.rds` format and used for downstream analysis and visualization.
  Should end in `.rds`. Default is "`mediation.rds`"

- total_cores:

  Number of available cores to use for parallelization. Default is
  `NULL`.

- save:

  Indicates object return/save behavior. One of `c("sr", "so", "ro")`;
  save & return, save only, return only. Default is "sr".

- distOnly:

  Logical. Mediate only the distal peaks? Default is `TRUE`.

- hsOnly:

  Logical. Mediate only on peaks identified within a hotspot. Default is
  `FALSE`.

## Value

A list containing mediation results for each tissue
