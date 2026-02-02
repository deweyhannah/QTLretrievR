# Plot peaks associated with specific genes.

Plot peaks associated with specific genes.

## Usage

``` r
peak_plot(
  mapping,
  tissue,
  pheno,
  pop = "do",
  outdir = NULL,
  pname = NULL,
  psave = TRUE,
  chromosome = NULL,
  effects = FALSE,
  founders = NULL,
  palette = NULL
)
```

## Arguments

- mapping:

  Mapping list from `mapQTL`, or full path to `.rds` containing one.

- tissue:

  Tissue to derive plots from

- pheno:

  Phenotype to plot

- pop:

  One of `c("do", "cc", "other")` to indicate founder population.
  Default is "`do`".

- outdir:

  Directory to save plots. Default is `NULL`.

- pname:

  File name to save plot (needs to end in `.png`).

- psave:

  Logical. Save the plot as `.png`. Default `TRUE`.

- chromosome:

  Chromosome that the peaks is present on.

- effects:

  Logical. Plotting effects on chromosome. Default is `FALSE`.

- founders:

  If `pop == "other"`, list of founders in haplotype order.

- palette:

  Founder color map (\> 8 founders).

## Value

Plot of phenotype specific peaks, or phenotype specific peak with
effects
