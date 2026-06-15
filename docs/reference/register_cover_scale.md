# Register a cover scale

Register a cover scale

## Usage

``` r
register_cover_scale(name, fn)
```

## Arguments

- name:

  Scale name.

- fn:

  Function `(x, mapping)` returning cover-proportion in `[0,1]`, with
  `NA` for missing or unrecognized values (imputed downstream).

## Value

Invisibly `NULL`; called for the side effect of registration.
