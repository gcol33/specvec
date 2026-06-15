# Co-occurrence operator (engine primitive)

Build the weighted species x species operator a factorization would
consume (for `chi_square`, return the plot x species contingency, whose
residual is applied implicitly at factorization). Exported for power
users; most callers use
[`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md).

## Usage

``` r
cooc_matrix(
  x,
  weighting = c("counts", "ppmi", "abundance_pmi", "chi_square"),
  min_occurrence = 5L,
  min_cooccurrence = 1L,
  time = NULL
)
```

## Arguments

- x:

  A `specvec_data` object.

- weighting:

  One of `"counts"`, `"ppmi"`, `"abundance_pmi"`, `"chi_square"`.

- min_occurrence:

  Drop species occurring in fewer than this many plots.

- min_cooccurrence:

  Keep a species pair only if it co-occurs in at least this many plots
  (species x species weightings only).

- time:

  Optional time window: restrict to plots whose stored `time` falls in
  the window before building the operator. Numeric time selects the
  closed interval `[min(time), max(time)]`; otherwise set membership.
  `NULL` uses all plots.

## Value

A sparse matrix: species x species for `counts`/`ppmi`/ `abundance_pmi`,
or plot x species for `chi_square`.

## Examples

``` r
df <- data.frame(plot = c("p1","p1","p2","p2","p3","p3"),
                 species = c("A","B","A","B","A","C"))
x <- specvec(df, "plot", "species")
cooc_matrix(x, "ppmi", min_occurrence = 1)
```
