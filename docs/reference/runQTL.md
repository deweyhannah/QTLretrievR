# Wrapper function to generate mapping, peaks, mediation, and effects data

Wrapper function to generate mapping, peaks, mediation, and effects data

## Usage

``` r
runQTL(
  geno_out = "gbrs_interpolated_genoprobs.rds",
  peaks_out = "peaks.rds",
  map_out = "mapping.rds",
  med_out = "mediation_res.rds",
  effects_out = "effects.rds",
  outdir = NULL,
  gbrs_fileLoc,
  metadata,
  expr_mats,
  covar_factors,
  annots,
  tissues = c(),
  gridFile = gridfile,
  total_cores = NULL,
  save_t = "sr",
  ...
)
```

## Arguments

- geno_out:

  String indicating the name of the output file containing interpolated
  genotype probabilities. This file will be saved in `.rds` format and
  used for QTL mapping. Should end in `.rds`. Default is
  "`gbrs_interpolated_genoprobs.rds`".

- peaks_out:

  String indicating the name for output peaks file. This file will be
  saved in `.rds` format and will be used as an input for downstream
  analysis. Should end in `.rds`. Default is "`peaks.rds`"

- map_out:

  String indicating the name of the output file containing QTL mapping
  results. This file will be saved in `.rds` format and will be used in
  downstream analyses. Should end in `.rds`. Default is "`map.rds`"

- med_out:

  String indicating the name of the output file containing mediation
  results for mediation within a phenotype. This file will be saved in
  `.rds` format and used for downstream analysis and visualization.
  Should end in `.rds`. Default is "`mediation.rds`"

- effects_out:

  String indicating the name of the output file containing founder
  haplotype effects results. This file will be saved in `.rds` format
  and used for downstream analysis and visualization. Should end in
  `.rds`. Default is "`effects.rds`"

- outdir:

  Directory to save output files. Default is `NULL`.

- gbrs_fileLoc:

  Path to GBRS interpolated tsv files.

- metadata:

  Sample metadata. Either a string pointing to the file, or the object
  itself.

- expr_mats:

  List of normalized count matrices (objects), or character paths to the
  file. One matrix per tissue. The order *must match* the tissue order
  in `genoprobs`.

- covar_factors:

  Additive covariate factors. These need to be columns in the sample
  metadata.

- annots:

  Annotations file. Contains mapping information for phenotypes.
  Dataframe, or tsv. Columns must include "id", "symbol", "start",
  "end".

- tissues:

  Vector of strings indicating tissues or conditions in project.

- gridFile:

  Genome Grid. Path to location or object. Defaults to 75k grid loaded
  with package.

- total_cores:

  Number of available cores to use for parallelization. Default is
  `NULL`.

- save_t:

  Indicates object return/save behavior. One of `c("sr", "so", "ro")`;
  save & return, save only, return only. Default is "sr".

## Value

A list containing

- peaks_listUnfiltered peaks for each tissue.

- maps_listList of objects associated with mapping. See
  [mapQTL](https://deweyhannah.github.io/QTLretrievR/reference/mapQTL.md)
  help for details.

- res_listList containing mediation results for each tissue.

- effects_resList of objects associated with effects. See
  [qtl_effects](https://deweyhannah.github.io/QTLretrievR/reference/qtl_effects.md)
  help for details.
