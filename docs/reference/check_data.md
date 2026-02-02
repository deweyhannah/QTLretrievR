# Checking the inputted data for QTLretrievR

Checking the inputted data for QTLretrievR

## Usage

``` r
check_data(x, type = "")
```

## Arguments

- x:

  can be a path to an Rds object or a list of R objects.

- type:

  The name of the object if a single one is passed. Use "genoprobs" for
  genotype probabilities, "mediation" for mediation results and "peaks"
  for a table of eQTL peaks.

## Value

List of properly names R objects to be used in downstream analysis of
QTLretrievR.
