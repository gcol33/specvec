# Species similarity

Similarity of one species to another, or to all species.

## Usage

``` r
species_similarity(x, species, to = NULL, metric = c("cosine", "euclidean"))
```

## Arguments

- x:

  A `specvec_embedding`.

- species:

  Focal species id.

- to:

  Optional second species id; if `NULL`, returns the named similarity
  (cosine) or distance (euclidean) of `species` to every species.

- metric:

  `"cosine"` (default) or `"euclidean"`.

## Value

A scalar when `to` is supplied, otherwise a named numeric vector.
