# Convert GBRS tsv genome probabilities to genoprobs format.

Converts GBRS-formatted TSV genotype probability files into a 3D array
format compatible with `qtl2`, organized by tissue and sample,
optionally saves the result as an RDS file

## Usage

``` r
genoprobably(
  outfile = "./gbrs_interpolated_genoprobs.rds",
  gbrsFileLoc,
  tissues = c(),
  gridFile = gridfile,
  save = "sr"
)
```

## Arguments

- outfile:

  File path to save the resulting genopbrobs list. Defaults to
  "gbrs_interpolated_genoprobs.rds".

- gbrsFileLoc:

  File path to where the GBRS files are located.

- tissues:

  Character vector of tissue names. These should match substrings in the
  GBRS file names. if empty, defaults fo "a".

- gridFile:

  Either a data frame or file path to the genome grid used to assign
  marker names. Default is 75K grid loaded with the package.

- save:

  Character. Determines output behavior: "sr" to save and return, "so"
  to save only, "ro" to return only. Default is "sr".

## Value

A named list of 3D genotype probability arrays (one per tissue),
formatted for use with `qtl2`.
