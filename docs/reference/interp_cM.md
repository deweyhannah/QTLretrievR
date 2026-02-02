# interp_cM

Interpolate bp peak positions to cM based on marker maps.

## Usage

``` r
interp_cM(df, genmap, physmap)
```

## Arguments

- df:

  A peak dataframe with the peak position in bp.

- genmap:

  Genetic map - found in `maps_list` (gmap).

- physmap:

  Physical map - found in `maps_list` (pmap).

## Value

Peak dataframe with interpolated bp column appended
