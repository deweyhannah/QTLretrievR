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
  phys = TRUE,
  min_cores = 4,
  BPPARAM = BiocParallel::SerialParam()
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
  in `genoprobs`. These *cannot* be rankZ transformed phenotypes, but
  they should be normalized and not raw.

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

- phys:

  Logical. if `TRUE`, use the physical map for peak calling; otherwise
  use the genomic map. Default is `TRUE`.

- min_cores:

  Numeric. The minimum number of cores to be used in each process for
  mapping.

- BPPARAM:

  BiocParallel Parameter

## Value

A list containing:

- maps_list:

  A list of dataframes and lists that can be used for future analyses
  and in other functions, containing:

  qtlprobs

  :   Genome probabilities in qtl2 format

  covar_list

  :   List of covariate matrices for each tissue

  expr_list

  :   Original normalized expression data for each tissue

  exprZ_list

  :   Rank Z-normalized expression data for each tissue

  kinship_loco

  :   Kinship matrix calculated using the "loco" option in
      [`qtl2::calc_kinship`](https://rdrr.io/pkg/qtl2/man/calc_kinship.html)

  gmap

  :   Genomic map of markers

  map_dat2

  :   Combined genomic and physical map of markers

  pmap

  :   Physical map of markers

  tissue_samp

  :   Metadata split out by tissue

- peaks_list:

  A list of peak objects for each tissue
