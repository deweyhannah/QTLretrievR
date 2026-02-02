# interp_bp

Interpolate bp peak positions to bp based on marker maps

## Usage

``` r
interp_bp(df, genmap, physmap)
```

## Arguments

- df:

  A peak dataframe with the peak position in cM.

- genmap:

  Genetic map - found in `maps_list` (gmap).

- physmap:

  Physical map - found in `maps_list` (pmap).

## Value

Peak dataframe with interpolated bp column appended
