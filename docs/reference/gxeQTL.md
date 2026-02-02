# Generate mapping data and peaks for G x E QTL analysis

This function performs G x E QTL mapping across multiple tissues,
generating genome-wide LOD scores and identifying significant peaks. It
supports parallel processing and flexible input formats for expression
data and covariates.

## Usage

``` r
gxeQTL(
  genoprobs,
  samp_meta,
  expr_mats,
  covar_factors,
  thrA = 5,
  thrX = 5,
  gridFile = gridfile,
  localRange = 2e+06,
  outdir = NULL,
  peaks_out = "gxe_peaks.rds",
  map_out = "gxe_map.rds",
  annots = NULL,
  total_cores = NULL,
  save = "sr",
  delta = FALSE,
  ctrl,
  env,
  rz = FALSE,
  ...
)
```

## Arguments

- genoprobs:

  Genotype probabilities object (qtl2-formatted), or path to .rds file
  containing one.

- samp_meta:

  Sample metadata. Either a string pointing to the file, or the object
  itself.

- expr_mats:

  List of normalized count matrices (objects), or character paths to the
  file. One matrix per tissue. The order *must match* the tissue order
  in `genoprobs`.

- covar_factors:

  Additive covariate factors. These need to be columns in the factor
  metadata.

- thrA:

  Minimum reported LOD threshold for autosomes. Default is 5.

- thrX:

  Minimum reported LOD threshold for X chromosome. Default is 5.

- gridFile:

  Genome Grid. Path to location or object. Defaults to 75k grid loaded
  with package.

- localRange:

  Definition of "local" in bp. Default is 2e6.

- outdir:

  Directory to save output files. Default is `NULL`.

- peaks_out:

  String indicating the name for output peaks file. This file will be
  saved in `.rds` format and and be used as an input for downstream
  analysis. Should end in `.rds`. Default is "`gxe_peaks.rds`"

- map_out:

  String indicating the name of the output file containing G x E QTL
  mapping results. This file will be saved in `.rds` format and will be
  used in downstream analyses. Should end in `.rds`. Default is
  "`gxe_map.rds`"

- annots:

  Annotations file. Contains mapping information for phenotypes.
  Dataframe, or tsv. Columns must include "id", "symbol", "start",
  "end".

- total_cores:

  Number of available cores to use for parallelization. Default is
  `NULL`.

- save:

  Indicates object return/save behavior. One of `c("sr", "so", "ro")`;
  save & return, save only, return only. Default is "sr".

- delta:

  Logical. The delta method for G x E analysis (env - ctrl) should be
  used. Default is `FALSE`.

- ctrl:

  String indicating your control or background gene name (ex: "ctrl" or
  "CTRL").

- env:

  String indicating your exposed/treated samples (ex: "trt" or "treated"
  or "\<your_treatment_here\>").

- rz:

  Logical. Set to `TRUE` if expression data is already
  rankZ-transformed. Default is `FALSE`.

## Value

A list containing:

- maps_listA list of dataframes and lists to that can be used for future
  analyses and in other functions

  - qtlprobsGenome probabilities in qtl2format

  - covar_listlist of covariate matrices for each tissue

  - expr_listOriginal normalized expression data for each tissue

  - exprZ_listRank Z normalized expression data for each tissue

  - kinship_locoKisnhip Matrix calculated using the "loco" option in
    [`qtl2::calc_kinship`](https://rdrr.io/pkg/qtl2/man/calc_kinship.html)

  - gmapGenomic map of markers

  - map_dat2Combined genomic and physical map of markers

  - pmapPhysical map of markers

  - tissue_sampMetadata broken down for each tissue

- peaks_listA list of peaks list for each tissue.
