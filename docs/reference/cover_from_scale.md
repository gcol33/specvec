# Convert an abundance column to cover-proportion

Convert an abundance column to cover-proportion

## Usage

``` r
cover_from_scale(
  x,
  scale = c("percent", "proportion", "braun_blanquet"),
  mapping = NULL
)
```

## Arguments

- x:

  Abundance values: numeric percent, numeric proportion, or ordinal
  cover-abundance codes.

- scale:

  One of `"percent"`, `"proportion"`, `"braun_blanquet"`, or any
  registered scale.

- mapping:

  Optional named numeric vector overriding the default lookup (used by
  ordinal scales such as Braun-Blanquet).

## Value

Numeric cover-proportion in `[0,1]`, `NA` where missing or unrecognized.

## Examples

``` r
cover_from_scale(c(0, 50, 100), scale = "percent")
cover_from_scale(c("r", "+", "2", "5"), scale = "braun_blanquet")
```
