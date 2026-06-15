# Community trajectory and novelty over time

Embed each time window's communities in one fixed frame and score how
novel they are relative to a reference window. The species frame is
fitted once on all plots, so the per-window community embeddings live in
the same space and novelty is comparable across time (the temporal
counterpart of
[`community_novelty()`](https://gcol33.github.io/specvec/reference/community_novelty.md)).

## Usage

``` r
community_trajectory(
  x,
  by = NULL,
  reference = NULL,
  method = "abund_pmi",
  dim = 64L,
  weights = c("cover", "presence"),
  normalize = FALSE,
  k = 5L,
  min_occurrence = 5L,
  min_cooccurrence = 1L,
  frame_embedding = NULL,
  ...
)
```

## Arguments

- x:

  A `specvec_data` object built with a `time` column.

- by:

  Window definition (see
  [`species_trajectory()`](https://gcol33.github.io/specvec/reference/species_trajectory.md)).

- reference:

  Reference window the novelty is measured against: `NULL` (default)
  uses the first window, a window label or index picks one, or a
  `specvec_community`/matrix supplies an external baseline.

- method, dim, min_occurrence, min_cooccurrence:

  Frame embedding controls.

- weights:

  `"cover"` (default when present) or `"presence"` pooling.

- normalize:

  If `TRUE`, L2-normalize each community vector before scoring.

- k:

  Neighbours averaged in the novelty distance.

- frame_embedding:

  Optional pre-fitted frame embedding to reuse.

- ...:

  Passed to
  [`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)
  when fitting the frame.

## Value

A `specvec_community_trajectory`: `communities` (per-window plot x dim
matrices in the shared frame), a `novelty` table (window, center,
n_plots, mean and median per-plot novelty vs the reference), the fixed
`frame`, and the resolved `reference` label.

## See also

[`species_trajectory()`](https://gcol33.github.io/specvec/reference/species_trajectory.md),
[`community_novelty()`](https://gcol33.github.io/specvec/reference/community_novelty.md)

## Examples

``` r
set.seed(1)
df <- data.frame(
  plot = rep(paste0("p", 1:90), each = 3),
  species = sample(paste0("s", 1:15), 270, replace = TRUE),
  decade = rep(c(1990, 2000, 2010), each = 90)
)
x <- specvec(df, "plot", "species", time = "decade")
community_trajectory(x, dim = 4, min_occurrence = 1, k = 3)
```
