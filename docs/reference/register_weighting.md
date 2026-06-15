# Register a weighting

A weighting builds the operator a factorization consumes. `fn` takes
`(data, kept_species, n_plots, min_cooccurrence)` and returns an
operator list with a `kind` field (`"sym"`, `"counts"`, or
`"implicit"`).

## Usage

``` r
register_weighting(
  name,
  fn,
  input = c("species_species", "plot_species"),
  supports_abundance = FALSE
)
```

## Arguments

- name:

  Weighting name.

- fn:

  Operator-building function.

- input:

  `"species_species"` or `"plot_species"`.

- supports_abundance:

  Logical; whether the weighting uses cover.

## Value

Invisibly `NULL`; called for the side effect of registration.
