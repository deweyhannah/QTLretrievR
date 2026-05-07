# Identify distal hotspots above a suggestive and significant LOD score.

Identify distal hotspots above a suggestive and significant LOD score.

## Usage

``` r
transbands(
  map_dat,
  peaks,
  sigLOD = 7.5,
  suggLOD = 6,
  psave = TRUE,
  pname = NULL,
  outdir = NULL,
  color = "#0073C2FF"
)
```

## Arguments

- map_dat:

  `map_dat2` from `mapQTL` mapping list.

- peaks:

  List of dataframes containing QTL peaks for each tissue (annotated).

- sigLOD:

  Significant LOD threshold to use for filtering phenotypes. Default is
  7.5

- suggLOD:

  Suggestive LOD threshold to use for filtering phenotypes. Default is
  6.

- psave:

  Logical. Save the plot as `.png`. Default `TRUE`.

- pname:

  File name to save plot (needs to end in `.png`). Default is
  `hotspots_<tissue>.png`

- outdir:

  Directory to save plots. Default is `NULL`.

- color:

  Plot color. Default is "#0073C2FF".

## Value

A list containing:

- bands.rna:

  A list of tibbles containing hotspot information for each tissue. Each
  tibble includes:

  Chromosome

  :   Chromosome of the hotspot

  Start/Stop

  :   Start and stop positions (bp) of the hotspot

  Significant peaks

  :   Number of significant peaks in the hotspot

  Suggestive peaks

  :   Number of suggestive peaks in the hotspot

- trans_band_plots:

  A list of `ggplot2` objects (one per tissue) showing the number of
  significant distal hotspot peaks.
