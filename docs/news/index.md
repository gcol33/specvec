# Changelog

## specvec 0.2.0

Temporal embeddings: track how species and communities move through
time.

### Trajectories

- [`species_trajectory()`](https://gcol33.github.io/specvec/reference/species_trajectory.md)
  traces focal species through a fixed embedding frame. One frame
  embedding is fitted once, then each focal species is placed per time
  window at the cover-weighted centroid of the species it co-occurs
  with. Every window shares one frame, so the points compare directly
  across time without alignment or rotation. The output carries a
  per-cell support count, so thin windows are visible rather than
  hidden.
- [`community_trajectory()`](https://gcol33.github.io/specvec/reference/community_trajectory.md)
  embeds each window’s communities in the shared frame and reports
  per-window novelty against a reference window, the temporal
  counterpart of
  [`community_novelty()`](https://gcol33.github.io/specvec/reference/community_novelty.md).

### Windowed embeddings

- [`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md),
  [`community_embedding()`](https://gcol33.github.io/specvec/reference/community_embedding.md),
  and
  [`cooc_matrix()`](https://gcol33.github.io/specvec/reference/cooc_matrix.md)
  gain a `time =` argument that fits on the plots in a time window. The
  reserved `time` seam from 0.1 is now live.

## specvec 0.1.0

First public release.

### Embeddings

- [`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md)
  builds a specvec object from a long plot-species table. Abundance is
  optional and read on the percent, proportion, or Braun-Blanquet scale
  (`cover_scale =`).
- [`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)
  learns species vectors from co-occurrence. The default `abund_pmi`
  method is abundance-weighted PMI; `ca`, `pmi`, and `glove` are
  available through the method registry.
- [`community_embedding()`](https://gcol33.github.io/specvec/reference/community_embedding.md)
  places each plot in the species space by pooling its species vectors,
  weighted by cover or presence.
- [`nearest_species()`](https://gcol33.github.io/specvec/reference/nearest_species.md),
  [`species_similarity()`](https://gcol33.github.io/specvec/reference/species_similarity.md),
  [`community_similarity()`](https://gcol33.github.io/specvec/reference/community_similarity.md),
  and
  [`community_novelty()`](https://gcol33.github.io/specvec/reference/community_novelty.md)
  query the learned space.

### Benchmark

- [`compare_embeddings()`](https://gcol33.github.io/specvec/reference/compare_embeddings.md)
  scores methods on your own data under one fair protocol: held-out
  co-occurrence recovery, link prediction, and optional habitat (EUNIS)
  and trait recovery. It reports which method wins by more than two
  pooled standard deviations.

### Extending

- A method is a registered `(weighting, factorization)` pair.
  [`register_weighting()`](https://gcol33.github.io/specvec/reference/register_weighting.md),
  [`register_factorization()`](https://gcol33.github.io/specvec/reference/register_factorization.md),
  and
  [`register_method()`](https://gcol33.github.io/specvec/reference/register_method.md)
  add new methods, and
  [`specvec_methods()`](https://gcol33.github.io/specvec/reference/specvec_methods.md)
  lists what is registered.

### Documentation

- [`vignette("specvec-methods")`](https://gcol33.github.io/specvec/articles/specvec-methods.md)
  derives the co-occurrence operator, the AbundPMI formula, and the
  baselines.
- [`vignette("specvec-benchmark")`](https://gcol33.github.io/specvec/articles/specvec-benchmark.md)
  runs
  [`compare_embeddings()`](https://gcol33.github.io/specvec/reference/compare_embeddings.md)
  on simulated data and on an EVA export.
