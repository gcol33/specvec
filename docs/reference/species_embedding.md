# Species embedding

Learn species vectors from co-occurrence. `method` is sugar for a
registered `(weighting, factorization)` pair; pass
`weighting`/`factorization` to override for power use. The default
`"abund_pmi"` is abundance-weighted PMI.

## Usage

``` r
species_embedding(
  x,
  method = "abund_pmi",
  dim = 64L,
  weighting = NULL,
  factorization = NULL,
  time = NULL,
  min_occurrence = 5L,
  min_cooccurrence = 1L,
  glove_iter = 20L
)
```

## Arguments

- x:

  A `specvec_data` object.

- method:

  One of
  [`specvec_methods()`](https://gcol33.github.io/specvec/reference/specvec_methods.md);
  default `"abund_pmi"`. Pass `NULL` to drive the engine purely from
  `weighting` + `factorization`.

- dim:

  Embedding dimension.

- weighting:

  Optional weighting name overriding the method's.

- factorization:

  Optional factorization name overriding the method's.

- time:

  Optional time window: fit on plots whose stored `time` falls in the
  window. Numeric time selects the closed interval
  `[min(time), max(time)]`; otherwise set membership. `NULL` (default)
  uses all plots. For trajectories comparable across windows, use
  [`species_trajectory()`](https://gcol33.github.io/specvec/reference/species_trajectory.md).

- min_occurrence:

  Drop species occurring in fewer than this many plots.

- min_cooccurrence:

  Keep a species pair only if it co-occurs in at least this many plots.

- glove_iter:

  Iterations for the GloVe factorizer.

## Value

A `specvec_embedding`: species x dim matrix `V` (species row names) plus
method, capability, and preprocessing metadata.

## Examples

``` r
df <- data.frame(plot = rep(paste0("p", 1:6), each = 2),
                 species = c("A","B","A","B","A","C","B","C","A","B","B","C"))
x <- specvec(df, "plot", "species")
species_embedding(x, method = "pmi", dim = 2, min_occurrence = 1)
```
