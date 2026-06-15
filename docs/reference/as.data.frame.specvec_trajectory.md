# Tidy a species trajectory

Long-to-wide data frame with one row per focal-species-by-window cell:
the window label and center, the support (plots backing the cell), and
the `dim` coordinate columns `d1..dD`.

## Usage

``` r
# S3 method for class 'specvec_trajectory'
as.data.frame(x, ..., na.rm = FALSE)
```

## Arguments

- x:

  A `specvec_trajectory`.

- ...:

  Unused.

- na.rm:

  Drop cells with no co-occurrence (all-`NA` coordinates).

## Value

A data frame.
