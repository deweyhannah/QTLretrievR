# Create a plot with founder effects for a given phenotype

Create a plot with founder effects for a given phenotype

## Usage

``` r
tailWag(
  effects,
  tissue,
  feat,
  chromosome,
  color = "#0073C2FF",
  psave = TRUE,
  pname = NULL,
  outdir = NULL,
  pop = "do",
  founders = NULL,
  symbol = TRUE
)
```

## Arguments

- effects:

  Effects object (output from `qtl_effects` function).

- tissue:

  String indicating which tissue to plot the founder effects for.

- feat:

  Individual phenotype to plot effects of.

- chromosome:

  Chromosome that the peak is present on.

- color:

  Plot color. Default is "#0073C2FF".

- psave:

  Logical. Save the plot as `.png`. Default `TRUE`.

- pname:

  File name to save plot (needs to end in `.png`). Default is
  `wag_effects_<tissue>_<feat>.png`

- outdir:

  Directory to save plots. Default is `NULL`.

- pop:

  One of `c("do", "cc", "other")` to indicate founder population.
  Default is "`do`".

- founders:

  If `pop == "other"`, list of founders in haplotype order.

- symbol:

  Logical. The phenotype to be plotted is passed as a symbol instead of
  phenotype id. Default `TRUE`.

## Value

ggplot object of the effects plot
