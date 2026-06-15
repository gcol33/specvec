# Package index

## Data

Build a specvec object from a long plot-species table and convert cover
scales.

- [`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md) :
  Build a specvec data object
- [`cover_from_scale()`](https://gcol33.github.io/specvec/reference/cover_from_scale.md)
  : Convert an abundance column to cover-proportion
- [`register_cover_scale()`](https://gcol33.github.io/specvec/reference/register_cover_scale.md)
  : Register a cover scale

## Species embeddings

Learn species vectors and query their neighbours.

- [`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)
  : Species embedding
- [`nearest_species()`](https://gcol33.github.io/specvec/reference/nearest_species.md)
  : Nearest species
- [`species_similarity()`](https://gcol33.github.io/specvec/reference/species_similarity.md)
  : Species similarity
- [`cooc_matrix()`](https://gcol33.github.io/specvec/reference/cooc_matrix.md)
  : Co-occurrence operator (engine primitive)

## Community embeddings

Place plots in the species space and score novelty.

- [`community_embedding()`](https://gcol33.github.io/specvec/reference/community_embedding.md)
  : Community embedding
- [`community_similarity()`](https://gcol33.github.io/specvec/reference/community_similarity.md)
  : Community similarity
- [`community_novelty()`](https://gcol33.github.io/specvec/reference/community_novelty.md)
  : Community novelty

## Temporal

Window embeddings by time and track species and communities through a
fixed frame.

- [`species_trajectory()`](https://gcol33.github.io/specvec/reference/species_trajectory.md)
  : Species trajectory through a fixed embedding frame
- [`community_trajectory()`](https://gcol33.github.io/specvec/reference/community_trajectory.md)
  : Community trajectory and novelty over time
- [`as.data.frame(`*`<specvec_trajectory>`*`)`](https://gcol33.github.io/specvec/reference/as.data.frame.specvec_trajectory.md)
  : Tidy a species trajectory

## Alien integration

Track a neophyte’s distance to the native community over time.

- [`integration_trajectory()`](https://gcol33.github.io/specvec/reference/integration_trajectory.md)
  : Alien integration trajectory
- [`as.data.frame(`*`<specvec_integration>`*`)`](https://gcol33.github.io/specvec/reference/as.data.frame.specvec_integration.md)
  : Tidy an integration trajectory

## Benchmark

Score embedding methods on your own data.

- [`compare_embeddings()`](https://gcol33.github.io/specvec/reference/compare_embeddings.md)
  [`specvec_benchmark()`](https://gcol33.github.io/specvec/reference/compare_embeddings.md)
  : Benchmark embedding methods

## Extend

Register new weightings, factorizations, and method presets.

- [`specvec_methods()`](https://gcol33.github.io/specvec/reference/specvec_methods.md)
  : List registered methods
- [`register_weighting()`](https://gcol33.github.io/specvec/reference/register_weighting.md)
  : Register a weighting
- [`register_factorization()`](https://gcol33.github.io/specvec/reference/register_factorization.md)
  : Register a factorization
- [`register_method()`](https://gcol33.github.io/specvec/reference/register_method.md)
  : Register a method preset
