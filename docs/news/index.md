# Changelog

## specvec 0.3.0

Alien integration, and a compositional method for the abundance fork.

### Integration trajectories

- [`integration_trajectory()`](https://gcol33.github.io/specvec/reference/integration_trajectory.md)
  reports a focal neophyte’s distance to the native-community centroid
  per time window, in one fixed embedding frame. The neophyte is placed
  by
  [`species_trajectory()`](https://gcol33.github.io/specvec/reference/species_trajectory.md)
  at the cover-weighted centroid of the species it co-occurs with; the
  native community is placed by pooling its plots in the same frame; the
  readout is the distance between them. A falling distance is
  integration. `native =` names the target community (defaults to the
  resident pool); `frame =` fixes the coordinate system so the pole the
  neophyte arrives from stays represented. Restrict the input to
  resampled (ReSurvey) plots before building the object so the same
  locations are tracked through time. See
  [`vignette("specvec-integration")`](https://gcol33.github.io/specvec/articles/specvec-integration.md).

### Compositional method

- `clr` adds the compositional route to abundance (the research fork of
  the methods grid): each plot is a composition, and the species
  operator is the centered-log-ratio covariance estimated from the
  variation matrix of pairwise log-ratios, double-centered. It registers
  through the same `(weighting, factorization)` machinery as the other
  methods. On the neutral co-occurrence metrics the cover-weighted PMI
  operand (`abund_pmi`) recovers more held-out signal, so it stays the
  default; `clr` is available for users who want the compositional
  treatment.
  [`vignette("specvec-methods")`](https://gcol33.github.io/specvec/articles/specvec-methods.md)
  describes it.

### API

- The query verbs now share one calling convention. The object queried
  is `x` in all of them: `nearest_species(x, ...)`,
  `species_similarity(x, ...)`, `community_similarity(x, ...)`,
  `community_novelty(x, ...)` (was `emb` / `object`). In
  [`species_similarity()`](https://gcol33.github.io/specvec/reference/species_similarity.md)
  the focal and comparison species are `species` and `to` (was `a` and
  `b`), matching `nearest_species(x, species)`. Positional calls are
  unaffected.
- `as_specvec()` is removed;
  [`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md)
  is the single constructor.

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
