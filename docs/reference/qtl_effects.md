# Calculate founder effects for significant peaks

Calculate founder effects for significant peaks

## Usage

``` r
qtl_effects(
  peaks,
  mapping,
  suggLOD = 6,
  outdir = NULL,
  effects_out = "effects.rds",
  total_cores = NULL,
  save = "sr"
)
```

## Arguments

- peaks:

  List of dataframes containing annotated QTL peaks for each tissue, or
  full path to `.rds` containing one.

- mapping:

  Mapping list from `mapQTL`, or full path to `.rds` containing one.

- suggLOD:

  Suggestive LOD threshold to use for filtering phenotypes. Default is
  6.

- outdir:

  Directory to save effects output files. Default is `NULL`.

- effects_out:

  String indicating the name of the output file containing founder
  haplotype effects results. This file will be saved in `.rds` format
  and used for downstream analysis and visualization. Should end in
  `.rds`. Default is "`effects.rds`"

- total_cores:

  Number of available cores to use for parallelization. Default is
  `NULL`.

- save:

  Indicates object return/save behavior. One of `c("sr", "so", "ro")`;
  save & return, save only, return only. Default is "sr".

## Value

A list containing:

- effects_blupQTL effect BLUPs from scan along one chromosome. Output
  from
  [`qtl2::scan1blup()`](https://rdrr.io/pkg/qtl2/man/scan1blup.html)

- peaksannotated peaks with LOD scores above suggestive threshold
