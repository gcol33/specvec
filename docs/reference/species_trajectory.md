# Species trajectory through a fixed embedding frame

Track focal species through time without the cross-time alignment
problem. A single *frame* embedding is fitted once (the stable
coordinate system), and each focal species is then placed, per time
window, at the cover-weighted centroid of the frame species it co-occurs
with in that window. Because every window is read out in one fixed
frame, the trajectory points are directly comparable across time – no
per-window rotation, no Procrustes alignment, no circularity. This is
the engine behind alien-integration trajectories (v0.3).

## Usage

``` r
species_trajectory(
  x,
  species,
  frame = NULL,
  by = NULL,
  method = "abund_pmi",
  dim = 64L,
  weights = c("cover", "presence"),
  min_occurrence = 5L,
  min_cooccurrence = 1L,
  frame_embedding = NULL,
  ...
)
```

## Arguments

- x:

  A `specvec_data` object built with a `time` column.

- species:

  Character vector of focal species ids to trace.

- frame:

  Optional character vector of species defining the fixed frame; default
  is every species except `species`.

- by:

  Window definition passed to the time splitter: `NULL` (default) for
  one window per distinct time value, or a numeric break vector to bin
  time.

- method:

  Frame embedding method (see
  [`specvec_methods()`](https://gcol33.github.io/specvec/reference/specvec_methods.md)).

- dim:

  Frame embedding dimension.

- weights:

  `"cover"` (default when cover is present) weights co-occurrence by the
  geometric mean of covers, matching AbundPMI; `"presence"` counts
  co-occurring plots.

- min_occurrence, min_cooccurrence:

  Frame species/pair filters.

- frame_embedding:

  Optional pre-fitted `specvec_embedding` to reuse as the frame instead
  of fitting one (its species become the frame).

- ...:

  Passed to
  [`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)
  when fitting the frame.

## Value

A `specvec_trajectory`: `U` (focal x window x dim array, `NA` where a
focal species has no co-occurrence in a window), `support` (focal x
window count of window plots containing the focal species), the fixed
`frame` embedding, the `windows` table, and pooling provenance.

## Details

The frame defaults to every species except the focal ones, so the focal
species move through the background community rather than helping define
it. Pass `frame` to fix the frame to a reference set (e.g. native
species only).

## See also

[`community_trajectory()`](https://gcol33.github.io/specvec/reference/community_trajectory.md),
[`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)

## Examples

``` r
set.seed(1)
df <- data.frame(
  plot = rep(paste0("p", 1:60), each = 3),
  species = sample(c("focal", paste0("s", 1:12)), 180, replace = TRUE),
  decade = rep(c(1990, 2000, 2010), each = 60)
)
x <- specvec(df, "plot", "species", time = "decade")
tr <- species_trajectory(x, species = "focal", dim = 4, min_occurrence = 1)
as.data.frame(tr)
```
