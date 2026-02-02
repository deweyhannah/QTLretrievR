# Plot QTL maps (peak vs gene)

Plot QTL maps (peak vs gene)

## Usage

``` r
plot_eqtlmap(
  map_dat,
  peaks,
  sigLOD = 7.5,
  outdir = NULL,
  pname = NULL,
  psave = TRUE,
  unit = "bp",
  map_col = "#0073C2FF"
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

- outdir:

  Directory to save plots. Default is `NULL`.

- pname:

  File name to save plot (needs to end in `.png`). Default is
  "eqtl_map_LOD\_.png".

- psave:

  Logical. Save the plot as `.png`. Default `TRUE`.

- unit:

  One of `c("bp", "mbp")`; annotation position units. Default is "bp".

- map_col:

  Plot color. Default is "#0073C2FF"

## Value

A list of peak maps (ggplot objects) for each tissue
