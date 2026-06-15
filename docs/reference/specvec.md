# Build a specvec data object

Turn a long `plot, species (, abundance)` table into the sparse matrices
the embedding engine needs. Column roles are passed by name. Abundance,
time, and labels are optional. A stored `time` column powers the
windowed embeddings and trajectories
([`species_trajectory()`](https://gcol33.github.io/specvec/reference/species_trajectory.md),
[`community_trajectory()`](https://gcol33.github.io/specvec/reference/community_trajectory.md)).

## Usage

``` r
specvec(
  data,
  plot,
  species,
  abundance = NULL,
  time = NULL,
  labels = NULL,
  cover_scale = c("percent", "proportion", "braun_blanquet"),
  duplicates = c("max", "sum", "first", "error"),
  cover_mapping = NULL
)

as_specvec(
  data,
  plot,
  species,
  abundance = NULL,
  time = NULL,
  labels = NULL,
  cover_scale = c("percent", "proportion", "braun_blanquet"),
  duplicates = c("max", "sum", "first", "error"),
  cover_mapping = NULL
)
```

## Arguments

- data:

  A data frame in long form (one row per plot-species record).

- plot:

  Name of the plot/site id column.

- species:

  Name of the species column.

- abundance:

  Optional name of the cover/abundance column.

- time:

  Optional name of a plot-level time column (e.g. decade or year),
  stored for windowed embeddings and trajectories.

- labels:

  Optional character vector of plot-level label column names (carried
  for the benchmark; the core embedding never needs them).

- cover_scale:

  How to read `abundance`: `"percent"` (default), `"proportion"`, or
  `"braun_blanquet"`. See
  [`cover_from_scale()`](https://gcol33.github.io/specvec/reference/cover_from_scale.md).

- duplicates:

  How to aggregate duplicated `plot x species` rows: `"max"` (default),
  `"sum"`, `"first"`, or `"error"`.

- cover_mapping:

  Optional named numeric vector overriding an ordinal cover-scale
  lookup.

## Value

A `specvec_data` object: sparse `P` (presence), `COV` (cover or `NULL`),
sorted `species`/`plots` id maps, optional `time`/`labels`, `meta`.

## Examples

``` r
df <- data.frame(
  plot = c("p1","p1","p2","p2","p3","p3"),
  species = c("A","B","A","B","A","C"),
  cover = c(40, 10, 60, 5, 30, 80)
)
specvec(df, plot = "plot", species = "species", abundance = "cover")
```
