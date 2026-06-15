# specvec

SpecVec learns abundance-aware species and community embeddings from ecological plot data.

A vegetation plot records which species grow together and how much cover each one has. SpecVec turns that co-occurrence into dense vectors: species that share habitat sit close together, and every plot becomes a point in the same space. Cover enters the weighting, so a dominant species carries more weight than a trace occurrence.

## Installation

```r
# install.packages("remotes")
remotes::install_github("gcol33/specvec")
```

## Quick start

```r
library(specvec)

## a long plot x species table, one row per record
data <- data.frame(
  plot    = c("p1","p1","p2","p2","p3","p3","p4","p4"),
  species = c("Acer campestre","Sambucus nigra",
              "Acer campestre","Prunus avium",
              "Sambucus nigra","Prunus avium",
              "Acer campestre","Sambucus nigra"),
  cover   = c(40, 10, 55, 20, 30, 60, 25, 15)
)

x    <- specvec(data, plot = "plot", species = "species", abundance = "cover")
emb  <- species_embedding(x, dim = 2, min_occurrence = 1)  # abund_pmi by default
nearest_species(emb, "Acer campestre")
comm <- community_embedding(x, embedding = emb)            # one vector per plot
```

`specvec()` takes a long table with one row per plot-species record. Cover is optional and read as percent by default; the proportion and Braun-Blanquet scales are built in (`cover_scale =`).

## What the neighbours look like

Fitted on a 10,000-plot European vegetation sample (EVA), the nearest neighbours of three established neophytes recover the communities each one grows in.

Robinia pseudoacacia is a nitrogen-fixing tree of dry, nutrient-rich woodland edges, and its neighbours are the hedgerow and ruderal-forest species that follow it:

```
nearest_species(emb, "Robinia pseudoacacia")
                species  similarity
     Euonymus europaeus       0.770
         Sambucus nigra       0.715
           Prunus avium       0.710
         Acer campestre       0.710
          Rubus caesius       0.705
```

Solidago gigantea and Impatiens glandulifera are two riverbank invaders, and they land in nearly the same floodplain neighbourhood of native willows, co-invading neophytes, and tall riparian herbs:

```
 Solidago gigantea                Impatiens glandulifera
   Humulus lupulus     0.838        Carduus crispus     0.857
   Acer negundo        0.810        Salix euxina        0.852
   Salix alba          0.801        Salix alba          0.839
   Echinocystis lobata 0.778        Humulus lupulus     0.823
   Solidago canadensis 0.772        Echinocystis lobata 0.807
```

The dry-edge tree and the two riparian herbs occupy distinct regions of the space, each beside the native and co-invading species it actually grows with. `vignette("specvec-benchmark")` shows how to build a `specvec` object from an EVA export and reproduce this.

## Methods and benchmark

`species_embedding()` and `community_embedding()` share one engine: a co-occurrence operator paired with a factorization. The default `abund_pmi` is abundance-weighted PMI. Correspondence analysis (`ca`), plain PMI (`pmi`), and `glove` reach the same call through a method registry, so adding a method is a registration.

`compare_embeddings()` scores the methods on your own data under one fair protocol (held-out co-occurrence, link prediction, and habitat recovery) and reports which method wins:

```r
compare_embeddings(x, methods = c("ca", "pmi", "abund_pmi", "glove"))
```

See `vignette("specvec-methods")` for the definitions and `vignette("specvec-benchmark")` to run the comparison.

## Status

v0.1 ships species and community embeddings, neighbour and novelty queries, and the benchmark. Temporal embeddings and alien-integration trajectories are on the roadmap. The core methods install from `Matrix` and `RSpectra` alone; `text2vec` is needed only for the GloVe method.

## License

MIT (c) Gilles Colling. Department of Botany and Biodiversity Research, University of Vienna.
