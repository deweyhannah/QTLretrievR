# Prepare and run between phenotype mediation for a set of QTL peaks with provided expression data.

Mediation scan between two different phenotypes. Provide peaks and
mapping for the phenotype to mediate, and the rank Z-transformed
phenotype quantification to use as mediators.

## Usage

``` r
multi_modiFinder(
  peaks,
  mapping,
  exprZ,
  suggLOD = 7,
  annots,
  outdir = NULL,
  med_out = "multi_pheno_mediation.rds",
  total_cores = NULL,
  save = "sr"
)
```

## Arguments

- peaks:

  List of dataframes containing QTL peaks for each tissue, or full path
  pointing to saved peaks object.

- mapping:

  Mapping list from `mapQTL`, or full path to `.rds` containing one.

- exprZ:

  rankZ transformed expression to use for the mediation.

- suggLOD:

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
  results for mediation between phenotypes. This file will be saved in
  `.rds` format and used for downstream analysis and visualization.
  Should end in `.rds`. Default is "`multi_pheno_mediation.rds`"

- total_cores:

  Number of available cores to use for parallelization. Default is
  `NULL`.

- save:

  Indicates object return/save behavior. One of `c("sr", "so", "ro")`;
  save & return, save only, return only. Default is "sr".

## Value

A list containing mediation results for each tissue
