# Benchmark embedding methods

Score registered embedding methods on neutral ecological tasks under one
fair protocol: a plot-level train/test split, one scored species-pair
set per seed, never-seen negatives for link prediction, and a single
presence-mean community definition shared across methods. Co-occurrence
metrics need only `x`; the EUNIS and trait metrics activate when labels
or a trait table are supplied. The species set is filtered once
(`min_occurrence`) before splitting, so every method sees the same
species. Runs on any user dataset, not just the reference data that
tuned the package default.

## Usage

``` r
compare_embeddings(
  x,
  methods = c("ca", "pmi", "abund_pmi", "glove"),
  metrics = c("cooc_ppmi", "cooc_raw", "link_auc", "eunis", "trait"),
  dim = 64L,
  seeds = 1:3,
  test_frac = 0.2,
  min_occurrence = 5L,
  min_cooccurrence = 1L,
  labels = NULL,
  traits = NULL,
  reference = "ca",
  glove_iter = 20L,
  control = .benchmark_control()
)

specvec_benchmark(
  x,
  methods = c("ca", "pmi", "abund_pmi", "glove"),
  metrics = c("cooc_ppmi", "cooc_raw", "link_auc", "eunis", "trait"),
  dim = 64L,
  seeds = 1:3,
  test_frac = 0.2,
  min_occurrence = 5L,
  min_cooccurrence = 1L,
  labels = NULL,
  traits = NULL,
  reference = "ca",
  glove_iter = 20L,
  control = .benchmark_control()
)
```

## Arguments

- x:

  A `specvec_data` object.

- methods:

  Method names to compare (see
  [`specvec_methods()`](https://gcol33.github.io/specvec/reference/specvec_methods.md)).
  `"glove"` is skipped with a note if `text2vec` is not installed.

- metrics:

  Which metrics to report: any of `"cooc_ppmi"`, `"cooc_raw"`,
  `"link_auc"`, `"eunis"`, `"trait"`. EUNIS/trait skip cleanly when
  their inputs are absent.

- dim:

  Embedding dimension.

- seeds:

  Integer seeds; each gives an independent split. Reported as
  across-seed mean and SD.

- test_frac:

  Held-out plot fraction (and species fraction for traits).

- min_occurrence, min_cooccurrence:

  Species/pair filters, applied globally before splitting.

- labels:

  Plot-level label for the EUNIS metric: a column name in `x$labels`, a
  plot-ordered vector, or `NULL` to use the first label column.

- traits:

  Species-level numeric attributes for the trait metric: a matrix or
  data frame with species row names, or `NULL` to skip.

- reference:

  Method the verdict compares against (default `"ca"`).

- glove_iter:

  Iterations for the GloVe factorizer.

- control:

  Scale guards and kNN neighbourhood; see the defaults in the source.
  Override only for very large data.

## Value

A `specvec_benchmark`: a tidy `method x metric` summary (mean, sd, n),
the per-seed raw rows, the run config, and the fixed-criterion verdict.
A method "beats" the reference only if it exceeds it on both neutral
metrics (`cooc_raw` and `link_auc`) by more than two pooled SDs.

## Examples

``` r
df <- data.frame(
  plot = rep(paste0("p", 1:40), each = 3),
  species = sample(paste0("s", 1:15), 120, replace = TRUE),
  cover = round(runif(120, 1, 100), 1)
)
x <- specvec(df, "plot", "species", abundance = "cover")
# \donttest{
compare_embeddings(x, methods = c("ca", "pmi", "abund_pmi"),
                   dim = 4, seeds = 1:2, min_occurrence = 1)
# }
```
