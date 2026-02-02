# Map founder haplotype effects for a given hotspot.

Map founder haplotype effects for a given hotspot.

## Usage

``` r
hsHapEffects(
  effects,
  tbands,
  chromosome,
  tissue,
  sigLOD,
  hsNum = 1,
  pop = "do",
  founders = NULL,
  palette = NULL,
  topFeats = NULL,
  psave = TRUE,
  pname = NULL,
  outdir = NULL,
  vert = FALSE,
  showTarg = TRUE,
  ...
)
```

## Arguments

- effects:

  Effects object `qtl_effects`.

- tbands:

  List of transband locations (`bands.rna` list from `transbands`
  function)

- chromosome:

  Chromosome that the transband (hotspot) is present on.

- tissue:

  Tissue that the transband (hotspot) is present in

- sigLOD:

  Significant LOD threshold to use for filtering phenotypes.

- hsNum:

  If there are multiple hotspots on a chromosome, indicate which one.
  Default is 1.

- pop:

  One of `c("do", "cc", "other")` to indicate founder population.
  Default is "`do`".

- founders:

  If `pop == "other"`, list of founders in haplotype order.

- palette:

  Founder color map.

- topFeats:

  Optional. Top hotspot features based on PCA analysis.

- psave:

  Logical. Save the plot as `.png`. Default `TRUE`.

- pname:

  File name to save plot (needs to end in `.png`). Default is
  `haplotype_effects_<tissue>_transband_<hsNum>_chromosome_<chromosome>.png`

- outdir:

  Directory to save plots. Default is `NULL`.

- vert:

  Logical. Rotate plot to vertical orientation. Default `FALSE`.

- showTarg:

  Logical. Show target names in resulting plot. Default `TRUE`.

- ...:

  Additional arguments to pass to ComplexHeatmap

## Value

List of plots. `ht_plot`: Plot only, formatted with legend on the left.
`ht`: Original ComplexHeatmap object.
