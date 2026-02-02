# Plot an individual peak or Principal Component of a hotspot.

Plot an individual peak or Principal Component of a hotspot.

## Usage

``` r
hsPeakPlot(
  mapping,
  feats,
  tbands,
  chromosome,
  tissue,
  candidateMed = NULL,
  hsNum = 1,
  color = "#0073C2FF",
  pc = FALSE,
  wag = FALSE,
  sigLOD,
  topPC = 0.5,
  psave = TRUE,
  pname = NULL,
  outdir = NULL,
  ...
)
```

## Arguments

- mapping:

  Mapping list from `mapQTL`.

- feats:

  List of phenotypes to include. Individual/Candidate Mediator Only:
  Local to hotspot. PC: Targets of hotspot.

- tbands:

  List of transband locations (`bands.rna` list from `transbands`
  function).

- chromosome:

  Chromosome that the transband (hotspot) is present on.

- tissue:

  Tissue that the transband (hotspot) is present in.

- candidateMed:

  Dataframe of candidate mediator information (can be from annotations).
  Columns must include "id", "symbol", "start", "end".

- hsNum:

  If there are multiple hotspots on a chromosome, indicate which one.
  Default is 1.

- color:

  Plot color. Default is "#0073C2FF".

- pc:

  Logical. Plot the first principal component of the targets. Default
  `FALSE`.

- wag:

  Logical. Include `tailWag` effects plot in final. Default is `FALSE`.

- sigLOD:

  Significant LOD threshold to use for filtering phenotypes.

- topPC:

  Proportion of top hotspot contributers to be identified with Principal
  Component Analysis. Default is 0.5.

- psave:

  Logical. Save the plot as `.png`. Default `TRUE`.

- pname:

  File name to save plot (needs to end in `.png`).

- outdir:

  Directory to save plots. Default is `NULL`.

- ...:

  Additional arguments

## Value

ggplot/patchwork object: Peak plot (with or without phenotype effects).
If plotting the Principal Component of the hotspot, also includes a list
of the top proportion of phenotypes contributing to the hotspot.
