# Community novelty

Per-plot novelty: the mean distance to the `k` nearest reference
communities. A plot far from everything in the reference set scores
high.

## Usage

``` r
community_novelty(object, reference, k = 5L)
```

## Arguments

- object:

  A `specvec_community` (the plots to score).

- reference:

  A `specvec_community` or matrix of reference communities.

- k:

  Number of nearest reference communities to average over.

## Value

A named numeric vector of novelty per plot.
