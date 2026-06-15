# Tidy an integration trajectory

Long data frame with one row per focal-species-by-window cell: the
window label and center, the plots in the window, the focal support, the
count of native-bearing plots, and the focal-to-native `distance`.

## Usage

``` r
# S3 method for class 'specvec_integration'
as.data.frame(x, ..., na.rm = FALSE)
```

## Arguments

- x:

  A `specvec_integration`.

- ...:

  Unused.

- na.rm:

  Drop cells with no measured distance.

## Value

A data frame.
