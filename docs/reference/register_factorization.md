# Register a factorization

Register a factorization

## Usage

``` r
register_factorization(
  name,
  fn,
  kind = c("matrix", "objective"),
  accepts = c("sym", "implicit", "counts")
)
```

## Arguments

- name:

  Factorization name.

- fn:

  Factorization function `(operator, dim, ...)` returning a species x
  dim matrix with species row names.

- kind:

  `"matrix"` (svd/eigen; composes with a symmetric weighting) or
  `"objective"` (brings its own loss and subsumes its weighting).

- accepts:

  Operator kind this factorization consumes (`"sym"`, `"implicit"`, or
  `"counts"`).

## Value

Invisibly `NULL`; called for the side effect of registration.
