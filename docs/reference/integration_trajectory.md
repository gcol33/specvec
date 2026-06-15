# Alien integration trajectory

Track how a neophyte moves toward (or away from) the native community
over time, in one fixed embedding frame. The focal species is placed per
time window at the cover-weighted centroid of the species it co-occurs
with (the
[`species_trajectory()`](https://gcol33.github.io/specvec/reference/species_trajectory.md)
projection), the native community is placed per window at the centroid
of its pooled plot embeddings
([`community_embedding()`](https://gcol33.github.io/specvec/reference/community_embedding.md)'s
readout), and the trajectory is the focal-to-native distance per window.
A falling distance is integration: the neophyte's associates shift from
where it arrived toward the resident native community.

## Usage

``` r
integration_trajectory(
  x,
  species,
  native = NULL,
  frame = NULL,
  by = NULL,
  method = "abund_pmi",
  dim = 64L,
  weights = c("cover", "presence"),
  metric = c("euclidean", "cosine"),
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

  Character vector of focal neophyte ids to trace.

- native:

  Character vector of native species defining the integration target;
  default is every frame species (the resident community at large).

- frame:

  Optional character vector fixing the coordinate system; default is
  every species except `species`.

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

  `"cover"` (default when cover is present) weights co-occurrence and
  pooling by cover; `"presence"` counts plots.

- metric:

  Distance from the focal species to the native centroid: `"euclidean"`
  (default) or `"cosine"` (returned as `1 - cosine similarity`, so a
  smaller value still means closer for both).

- min_occurrence, min_cooccurrence:

  Frame species/pair filters.

- frame_embedding:

  Optional pre-fitted `specvec_embedding` to reuse as the frame instead
  of fitting one.

- ...:

  Passed to
  [`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)
  when fitting the frame.

## Value

A `specvec_integration`: `distance` (focal x window distance to the
native centroid, `NA` where the focal has no co-occurrence or a window
holds no native community), `support` (focal x window plot counts), the
per-window `native_centroid` and `native_support`, the `windows` table,
the underlying `specvec_trajectory`, and the fixed `frame`.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) tidies
it.

## Details

Two species sets shape the measurement.

- `frame` is the fixed coordinate system, fitted once. It defaults to
  every species except the focal ones, so the neophyte moves through the
  resident pool rather than helping define it. Keeping the whole
  resident pool in the frame (not the native subset alone) means the
  pole the neophyte starts from – the disturbed or ruderal flora it
  arrives with – is represented, so an early position is well defined
  and the motion toward the natives is real.

- `native` is the subset of the frame whose centroid is the integration
  target. It defaults to the whole frame (the resident community at
  large); for a meaningful measurement supply the species that are
  native at the study region (e.g. the rows flagged native in an EVA
  `STATUS` column).

**ReSurvey anchoring.** Raw multi-decade plot data confounds integration
with where and what was sampled each decade. Restrict `x` to resampled
plots (the `ReSurvey plot (Y/N)` flag in EVA header exports) before
building the `specvec_data`, so the same locations are tracked through
time. The function takes whatever plots it is given and does not
hard-code a survey-design column; the restriction is a data-preparation
step, demonstrated in the science vignette.

## See also

[`species_trajectory()`](https://gcol33.github.io/specvec/reference/species_trajectory.md),
[`community_trajectory()`](https://gcol33.github.io/specvec/reference/community_trajectory.md)

## Examples

``` r
set.seed(1)
native  <- paste0("nat", 1:6)
ruderal <- paste0("rud", 1:6)
rows <- list(); pid <- 0L
for (dec in c(1990, 2000, 2010)) {
  share_nat <- (dec - 1990) / 20             # the neophyte shifts ruderal -> native
  for (i in 1:50) {
    pid <- pid + 1L
    pool <- if (stats::runif(1) < share_nat) native else ruderal
    rows[[pid]] <- data.frame(plot = paste0("p", pid),
                              species = c(sample(pool, 3), "alien"), decade = dec)
  }
  for (i in 1:30) {                          # resident native backbone, every decade
    pid <- pid + 1L
    rows[[pid]] <- data.frame(plot = paste0("p", pid),
                              species = sample(native, 3), decade = dec)
  }
}
df <- do.call(rbind, rows)
x <- specvec(df, "plot", "species", time = "decade")
it <- integration_trajectory(x, species = "alien", native = native,
                             dim = 4, weights = "presence", min_occurrence = 1)
as.data.frame(it)
```
