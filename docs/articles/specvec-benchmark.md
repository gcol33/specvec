# Compare embedding methods on your own data

SpecVec ships its bake-off as a function so you can score the embedding
methods on your community data, not only on the data that tuned the
package default.
[`compare_embeddings()`](https://gcol33.github.io/specvec/reference/compare_embeddings.md)
fits each method on training plots, scores it on held-out plots, and
reports which method recovers the most held-out signal.

``` r

library(specvec)
```

## A worked example on simulated plots

The simulation below draws species presence from smooth latent niches,
so species with nearby niches co-occur and carry cover that tracks how
central the plot is to each niche.

``` r

sim_plots <- function(M = 500, S = 60, seed = 1) {
  set.seed(seed)
  mu  <- matrix(rnorm(S * 2), S, 2)
  env <- matrix(rnorm(M * 2), M, 2)
  rows <- lapply(seq_len(M), function(p) {
    d2   <- rowSums((mu - matrix(env[p, ], S, 2, byrow = TRUE))^2)
    prob <- exp(-d2 / 2)
    present <- which(runif(S) < prob)
    if (length(present) < 2) present <- order(prob, decreasing = TRUE)[1:2]
    data.frame(plot = paste0("p", p), species = paste0("sp", present),
               cover = round(100 * prob[present] / max(prob[present]), 1))
  })
  do.call(rbind, rows)
}

df <- sim_plots(M = 500, S = 60)
x  <- specvec(df, "plot", "species", abundance = "cover")
x
#> <specvec_data> plots=500  species=60
#>   presence: nnz=10352  density=34.5067%
#>   abundance: yes (cover_scale=percent)  duplicates=max
```

Run the comparison. The co-occurrence metrics need only `x`:

``` r

b <- compare_embeddings(
  x,
  methods = c("ca", "pmi", "abund_pmi"),
  metrics = c("cooc_ppmi", "cooc_raw", "link_auc"),
  dim = 16, seeds = 1:3, min_occurrence = 3
)
b
#> <specvec_benchmark> plots=500  species=60  dim=16  seeds=1,2,3
#>   methods: ca, pmi, abund_pmi
#> 
#>   method            cooc_ppmi         cooc_raw
#>   abund_pmi      0.827+-0.028     0.652+-0.022
#>   pmi            0.820+-0.028     0.589+-0.029
#>   ca             0.784+-0.031     0.502+-0.023
#> 
#>   verdict vs 'ca' (beats on cooc_raw & link_auc by >2 pooled SDs):
#>     reference 'ca' absent; ordering only.
```

## Reading the output

Each method is fit on the same training plots and scored on the same
pairs:

- `cooc_raw`: Spearman correlation between the embedding’s pair scores
  and held-out raw co-occurrence counts. This is the neutral
  co-occurrence metric.
- `link_auc`: how well the pair scores rank co-occurring pairs above
  pairs that never co-occur in either split. An AUC of 0.5 is chance.
- `cooc_ppmi`: the same Spearman against held-out PPMI. PMI methods
  optimise a PPMI-style objective, so this is reported for reference.

The verdict applies a fixed rule decided before looking: a method beats
the reference (`ca` by default) only if it exceeds it on both neutral
metrics (`cooc_raw` and `link_auc`) by more than two pooled standard
deviations.

The latent-niche simulation is the unimodal-gradient regime
correspondence analysis is built for, so CA scores well here. The
abundance-weighted advantage that sets `abund_pmi` as the package
default was measured on continental vegetation data at 10,000 and
200,000 plots, where co-occurrence is sparse and the species set is
large.

## Habitat and trait recovery

Two more metrics activate when you supply labels. EUNIS habitat is
plot-level: pass a label column carried on the `specvec` object and add
`"eunis"` to `metrics`. Traits are species-level: pass a
species-by-trait matrix to `traits` and add `"trait"`. Both score
through k-nearest-neighbour recovery and skip cleanly when their inputs
are absent.

``` r

b <- compare_embeddings(
  x,
  methods = c("ca", "pmi", "abund_pmi", "glove"),
  metrics = c("cooc_raw", "link_auc", "eunis", "trait"),
  labels  = "habitat",     # a label column on `x`
  traits  = trait_matrix,  # rows named by species
  dim = 64, seeds = 1:3
)
```

## Using your own data

[`compare_embeddings()`](https://gcol33.github.io/specvec/reference/compare_embeddings.md)
takes any `specvec` object, so the bake-off runs on a vegetation
database export. A European Vegetation Archive (EVA) species table with
a plot id, a taxon name, and a cover column becomes a `specvec` like
this:

``` r

library(data.table)
sp <- fread("species_export.csv")
sp <- sp[rank == "species"]                       # taxonomic scope is the caller's

x <- specvec(sp, plot = "PlotObservationID", species = "taxon",
             abundance = "cover", labels = "Eunis_lvl2")

compare_embeddings(x, dim = 64, seeds = 1:3, min_occurrence = 5)
```

The package carries no vegetation data of its own; the comparison is
meant to run on yours.
