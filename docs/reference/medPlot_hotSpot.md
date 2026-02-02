# Identify top mediators within a hotspot (molecular QTL)

Identify top mediators within a hotspot (molecular QTL)

## Usage

``` r
medPlot_hotSpot(
  peaks,
  meds,
  tbands,
  chromosome,
  sigLOD,
  hsNum = 1,
  top_n = 5,
  plot = "padj",
  topFeats = NULL,
  psave = TRUE,
  pname = NULL,
  outdir = NULL,
  showTarg = T
)
```

## Arguments

- peaks:

  Dataframe of peaks to pull targets from (output of `mapQTL`, then
  select `peaks_list$tissue`).

- meds:

  Dataframe of mediation to pull potential mediators from (output from
  `modiFinder`, then select `$tissue`).

- tbands:

  Dataframe of all transbands (output of `transbands`, then select
  `bands.rna$tissue`).

- chromosome:

  Chromosome that the transband (hotspot) is present on.

- sigLOD:

  Significant LOD threshold to use for filtering phenotypes.

- hsNum:

  If there are multiple hotspots on a chromosome, indicate which one.
  Default is 1.

- top_n:

  The number of top mediators per target to show, determined by LOD
  drop. Default is 5. *NB* If you want to see all the mediators, set
  this to the number of mediators in your hotspot.

- plot:

  One of c("padj", "pval", "per_drop", "ranks") depending on what
  statistic should be plotted in the heatmap. Default "padj".

- topFeats:

  Optional. Top hotspot features based on PC analysis.

- psave:

  Logical. Save the plot as `.png`. Default `TRUE`.

- pname:

  String indicating the name of plot to save as a .png. Default
  "mediation_plot_chr*top*\<top_n\>.png".

- outdir:

  Directory to save plots. Default is `NULL`.

- showTarg:

  Logical. Show target names in resulting plot. Default `TRUE`.

## Value

A list containing a dataframe representing the ranking of each mediator
within the hotspot for each target of the hotspot and the heatmap
object.
