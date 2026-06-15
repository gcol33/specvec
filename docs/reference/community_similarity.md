# Community similarity

Pairwise similarity among plot embeddings, or of each plot to a
reference set.

## Usage

``` r
community_similarity(
  object,
  reference = NULL,
  metric = c("cosine", "euclidean")
)
```

## Arguments

- object:

  A `specvec_community`.

- reference:

  Optional `specvec_community` or matrix to compare against; if `NULL`,
  compares the object's plots to themselves.

- metric:

  `"cosine"` (default) returns a similarity matrix; `"euclidean"`
  returns a distance matrix.

## Value

A matrix of `object` rows by `reference` rows.
