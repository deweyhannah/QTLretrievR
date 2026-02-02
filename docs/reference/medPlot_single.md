# Identify mediators for provided features +/- x Mb from a given position on a given chromosome

Identify mediators for provided features +/- x Mb from a given position
on a given chromosome

## Usage

``` r
medPlot_single(
  meds,
  range = 2,
  position,
  feats,
  chromosome,
  top_n = 5,
  psave = TRUE,
  pname = NULL,
  outdir = NULL,
  plot = "padj",
  showTarg = TRUE
)
```

## Arguments

- meds:

  Dataframe of mediation to pull potential mediators from (output from
  `modiFinder`, then select `$tissue`).

- range:

  The range from the position which mediatiors should be pulled. In Mbp.

- position:

  The position of the peak(s) to be investigated (single position only).
  In bp.

- feats:

  Vector of targets to identify mediators for.

- chromosome:

  Chromosome that the peaks(s) is present on.

- top_n:

  The number of top mediators per target to show. Default 5.

- psave:

  Logical. Save the plot as `.png`. Default `TRUE`.

- pname:

  String indicating the name of plot to save as a .png. Default
  "mediation_plot\_*chrtop*\<top_n\>.png".

- outdir:

  Directory to save plots. Default is `NULL`.

- plot:

  One of c("padj", "pval", "per_drop", "ranks") depending on what
  statistic should be plotted in the heatmap. Default "padj".

- showTarg:

  Logical. Show target names in resulting plot. Default `TRUE`.

## Value

A list containing a dataframe representing the ranking of each mediator
within the region for each of the provided targets and the heatmap
object.
