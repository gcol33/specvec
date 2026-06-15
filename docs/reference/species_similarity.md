# Species similarity

Similarity of one species to another, or to all species.

## Usage

``` r
species_similarity(emb, a, b = NULL, metric = c("cosine", "euclidean"))
```

## Arguments

- emb:

  A `specvec_embedding`.

- a:

  Focal species id.

- b:

  Optional second species id; if `NULL`, returns the named similarity
  (cosine) or distance (euclidean) of `a` to every species.

- metric:

  `"cosine"` (default) or `"euclidean"`.

## Value

A scalar when `b` is supplied, otherwise a named numeric vector.
