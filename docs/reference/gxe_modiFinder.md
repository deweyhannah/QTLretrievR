# Prepare and run mediation for delta G x E peaks

Prepare and run mediation for delta G x E peaks

## Usage

``` r
gxe_modiFinder(
  peaks,
  mapping,
  by = "delta",
  sigLOD = 7.5,
  annots,
  outdir = NULL,
  med_out = "gxe_mediation.rds",
  total_cores = NULL,
  save = "sr",
  distOnly = TRUE,
  hsOnly = FALSE,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- peaks:

  List of dataframes containing delta G x E peaks (name 'delta')

- mapping:

  Mapping list from `gxeQTL`. Needs to contain mapping information for
  at minimum the delta condition, and the condition being used for
  mediation

- by:

  String indicating the name of the condition to be used for mediation.
  Default 'delta'.

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
  Should end in `.rds`. Default is "`delta_mediation.rds`"

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

- BPPARAM:

  BiocParallel Parameter

## Value

A list containing mediation results for the peaks called from the delta
G x E method. Mediated by 'delta', '', or ''.
