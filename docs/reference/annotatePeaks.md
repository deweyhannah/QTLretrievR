# Interpolate peak positions and attach annotations

Interpolate peak positions and attach annotations

## Usage

``` r
annotatePeaks(mapping, peaks, annots, localRange = 2e+06)
```

## Arguments

- mapping:

  Mapping list from `mapQTL`.

- peaks:

  List of un-annotated peaks.

- annots:

  Annotations file. Contains mapping information for phenotypes.
  Dataframe, or tsv. Columns must include "id", "symbol", "start",
  "end".

- localRange:

  Definition of "local" in bp. Default is 2e6.

## Value

List of annotated peaks for each tissue.
