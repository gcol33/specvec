# Community embedding

Place each plot in the species embedding space by pooling the vectors of
the species it contains. The readout is uniform across methods (pooled
species vectors), so plot embeddings stay comparable.
`weights = "cover"` pools by cover (the abundance moat at plot level);
`weights = "presence"` is the plain mean of present species vectors.

## Usage

``` r
community_embedding(
  x,
  method = "abund_pmi",
  dim = 64L,
  weights = c("cover", "presence"),
  embedding = NULL,
  normalize = FALSE,
  time = NULL,
  min_occurrence = 5L,
  min_cooccurrence = 1L,
  ...
)
```

## Arguments

- x:

  A `specvec_data` object.

- method:

  Embedding method (see
  [`specvec_methods()`](https://gcol33.github.io/specvec/reference/specvec_methods.md));
  ignored if `embedding` is supplied.

- dim:

  Embedding dimension; ignored if `embedding` is supplied.

- weights:

  `"cover"` (default when cover is present) or `"presence"`.

- embedding:

  Optional pre-fitted `specvec_embedding` to pool, instead of fitting
  one here.

- normalize:

  If `TRUE`, L2-normalize each plot vector after pooling.

- time:

  Optional time window (see
  [`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)):
  pool (and, when fitting, train on) only plots in the window. `NULL`
  (default) uses all plots. For per-window communities in one shared
  frame, use
  [`community_trajectory()`](https://gcol33.github.io/specvec/reference/community_trajectory.md).

- min_occurrence, min_cooccurrence:

  Passed to
  [`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)
  when fitting.

- ...:

  Passed to
  [`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md).

## Value

A `specvec_community`: plot x dim matrix `U` (plot row names) plus
pooling and provenance metadata.

## Examples

``` r
df <- data.frame(plot = rep(paste0("p", 1:8), each = 2),
  species = c("A","B","A","B","A","C","B","C","A","B","B","C","A","C","A","B"),
  cover = c(80,20,50,50,60,40,50,50,70,30,40,60,90,10,55,45))
x <- specvec(df, "plot", "species", abundance = "cover")
community_embedding(x, method = "abund_pmi", dim = 3, min_occurrence = 1)
```
