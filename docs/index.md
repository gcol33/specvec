# specvec

*species and communities as points in one space*

[![R-CMD-check](https://github.com/gcol33/specvec/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gcol33/specvec/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Abundance-weighted PMI embeddings of species co-occurrence, with
first-class community vectors and a method bake-off.**

Feed it a long `plot x species (x cover)` table. `specvec` turns the
co-occurrence into dense vectors through a sparse eigendecomposition
(`Matrix` + `RSpectra`): every species becomes a point, every plot
becomes a point in the same space, and species that keep the same
company sit close together. Cover enters the weighting, so a dominant
species counts for more than a trace occurrence. A plot-by-species table
records who occurs where; `specvec` reads off which species behave alike
and places communities and species in one comparable space.

``` r

library(specvec)

# long table: one row per plot-species record, cover optional
data <- data.frame(
  plot    = c("p1","p1","p2","p2","p3","p3","p4","p4"),
  species = c("Acer campestre","Sambucus nigra","Acer campestre","Prunus avium",
              "Sambucus nigra","Prunus avium","Acer campestre","Sambucus nigra"),
  cover   = c(40, 10, 55, 20, 30, 60, 25, 15)
)

x    <- specvec(data, plot = "plot", species = "species", abundance = "cover")
emb  <- species_embedding(x, dim = 2, min_occurrence = 1)   # abund_pmi by default
nearest_species(emb, "Acer campestre")
comm <- community_embedding(x, embedding = emb)             # one vector per plot
```

## What the neighbours look like

Fitted on a 10,000-plot sample of European vegetation, the nearest
neighbours of three established neophytes recover the communities each
one grows in.

Robinia pseudoacacia is a nitrogen-fixing tree of dry, nutrient-rich
woodland edges, and its neighbours are the hedgerow and ruderal-forest
species that follow it:

    nearest_species(emb, "Robinia pseudoacacia")
                    species  similarity
         Euonymus europaeus       0.770
             Sambucus nigra       0.715
               Prunus avium       0.710
             Acer campestre       0.710
              Rubus caesius       0.705

Solidago gigantea and Impatiens glandulifera are two riverbank invaders,
and they land in nearly the same floodplain neighbourhood of native
willows, co-invading neophytes, and tall riparian herbs:

     Solidago gigantea                Impatiens glandulifera
       Humulus lupulus     0.838        Carduus crispus     0.857
       Acer negundo        0.810        Salix euxina        0.852
       Salix alba          0.801        Salix alba          0.839
       Echinocystis lobata 0.778        Humulus lupulus     0.823
       Solidago canadensis 0.772        Echinocystis lobata 0.807

The dry-edge tree and the two riparian herbs occupy distinct regions of
the space, each beside the native and co-invading species it actually
grows with. These results come from European Vegetation Archive (EVA)
data, an access-controlled database available on request rather than as
an open download. `specvec` ships no vegetation data;
[`vignette("specvec-benchmark")`](https://gcol33.github.io/specvec/articles/specvec-benchmark.md)
reproduces the workflow on simulated plots and shows how to map an EVA
export onto
[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md)
once you have access.

## Abundance-aware, and benchmarked

The default method, `abund_pmi`, is pointwise mutual information on a
cover-weighted co-occurrence operator: each co-occurring plot
contributes the geometric mean of the two species’ covers, so an
abundant pair weighs more than a pair where one species is a trace.
[`compare_embeddings()`](https://gcol33.github.io/specvec/reference/compare_embeddings.md)
scores the registered methods on neutral ecological tasks under one
fixed protocol (held-out co-occurrence recovery, link prediction, and
habitat recovery) and reports which recovers the most held-out signal.

``` r

compare_embeddings(x, methods = c("ca", "pmi", "abund_pmi", "glove"))
```

The win criterion is fixed before the numbers are seen: a method beats
the reference only if it exceeds it on both neutral co-occurrence
metrics by more than two pooled standard deviations. On continental
vegetation data this selects `abund_pmi`, and its lead over
correspondence analysis widens with the number of plots.

## What’s in the box

- **[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md)**:
  build the data object from a long table, with percent, proportion, and
  Braun-Blanquet cover scales.
- **[`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)**:
  species vectors through a registered `(weighting, factorization)`
  method (`abund_pmi`, `pmi`, `ca`, `glove`, `clr`).
- **[`community_embedding()`](https://gcol33.github.io/specvec/reference/community_embedding.md)**:
  first-class plot vectors, cover- or presence-pooled, in the same space
  as the species.
- **[`nearest_species()`](https://gcol33.github.io/specvec/reference/nearest_species.md)
  /
  [`species_similarity()`](https://gcol33.github.io/specvec/reference/species_similarity.md)**:
  neighbour and similarity queries on the species vectors.
- **[`community_similarity()`](https://gcol33.github.io/specvec/reference/community_similarity.md)
  /
  [`community_novelty()`](https://gcol33.github.io/specvec/reference/community_novelty.md)**:
  plot-to-plot similarity and distance-to-reference novelty.
- **[`species_trajectory()`](https://gcol33.github.io/specvec/reference/species_trajectory.md)
  /
  [`community_trajectory()`](https://gcol33.github.io/specvec/reference/community_trajectory.md)
  /
  [`integration_trajectory()`](https://gcol33.github.io/specvec/reference/integration_trajectory.md)**:
  change through time in one fixed frame, including a neophyte’s
  distance to the native community.
- **[`compare_embeddings()`](https://gcol33.github.io/specvec/reference/compare_embeddings.md)**:
  the shipped method bake-off with a fixed verdict.

A method is a registered `(weighting, factorization)` pair, so a new
method is added by registering a weighting or a factorization, not by
copying a function.

## Installation

``` r

install.packages("pak")
pak::pak("gcol33/specvec")
```

The core methods need only `Matrix` and `RSpectra`. `text2vec` is
required only for the `glove` method, and `FNN` only to speed up
nearest-neighbour queries.

## Documentation

- [Getting
  started](https://gcol33.github.io/specvec/articles/specvec-quickstart.html)
- [Preparing data and cover
  scales](https://gcol33.github.io/specvec/articles/specvec-data.html)
- [Methods: from plots to
  embeddings](https://gcol33.github.io/specvec/articles/specvec-methods.html)
- [Community embeddings and
  novelty](https://gcol33.github.io/specvec/articles/specvec-communities.html)
- [Tracking species and communities through
  time](https://gcol33.github.io/specvec/articles/specvec-temporal.html)
- [Alien integration
  trajectories](https://gcol33.github.io/specvec/articles/specvec-integration.html)
- [Comparing embedding
  methods](https://gcol33.github.io/specvec/articles/specvec-benchmark.html)
- [Extending
  specvec](https://gcol33.github.io/specvec/articles/specvec-extending.html)
- [Function reference](https://gcol33.github.io/specvec/reference/)

## Support

> “Software is like sex: it’s better when it’s free.” — Linus Torvalds

I’m a PhD student who builds R packages in my free time because I
believe good tools should be free and open. I started these projects for
my own work and figured others might find them useful too.

If this package saved you some time, buying me a coffee is a nice way to
say thanks. It helps with my coffee addiction.

[![Buy Me A
Coffee](https://img.shields.io/badge/-Buy%20me%20a%20coffee-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/gcol33)

## License

MIT (c) Gilles Colling. Department of Botany and Biodiversity Research,
University of Vienna.

## Citation

``` bibtex
@software{specvec,
  author = {Colling, Gilles},
  title  = {specvec: Abundance-Aware Species and Community Embeddings},
  year   = {2026},
  url    = {https://github.com/gcol33/specvec}
}
```
