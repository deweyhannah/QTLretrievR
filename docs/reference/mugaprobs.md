# Convert and process MUGA probabilities to qtl2 format

Converts MUGA genotyping final reports into a 3D array format compatible
with `qtl2`, organized by tissue and sample, optionally saves the result
as an RDS file

## Usage

``` r
mugaprobs(
  type = "GM",
  covarLoc,
  covar_file,
  i.files,
  genoPrefix = "gm4qtl2",
  probsOut = "muga_interpolated_genoprobs.rds",
  saveDir = getwd(),
  tissues = c()
)
```

## Arguments

- type:

  One of `c("GM", or "MM")`; GigaMUGA or MegaMUGA. Default is "GM".

- covarLoc:

  Location of covariate file.

- covar_file:

  Covariate file including at minimum the sex (sex) and generation
  (ngen) of each sample, this needs to be a .csv file

- i.files:

  Either a string of the directory where the chromosome specific
  genotype files are or a list of final report files to process - if
  passing the final report files they need to be either unzipped or in
  .gz format

- genoPrefix:

  Prefix for the chromosome specific genotype files (excluding "\_geno")

- probsOut:

  File name to save probabilities, default is
  "muga_interpolated_genoprobs.rds".

- saveDir:

  Directory to save genome probability object. Default is the current
  working directory.

- tissues:

  List of tissues included in analysis. If left blank tissue will be set
  to "a".

## Value

none
