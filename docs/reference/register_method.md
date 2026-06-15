# Register a method preset

Maps a friendly method name to a `(weighting, factorization)` pair plus
a capability descriptor used for dispatch and validation.

## Usage

``` r
register_method(
  name,
  weighting,
  factorization,
  input = c("species_species", "plot_species"),
  native_output = c("species", "species_and_community"),
  supports_abundance = FALSE,
  default_dim = 64L
)
```

## Arguments

- name:

  Method name (e.g. `"abund_pmi"`).

- weighting:

  Registered weighting name.

- factorization:

  Registered factorization name.

- input:

  `"species_species"` or `"plot_species"`.

- native_output:

  `"species"` or `"species_and_community"`.

- supports_abundance:

  Logical.

- default_dim:

  Default embedding dimension.

## Value

Invisibly `NULL`; called for the side effect of registration.
