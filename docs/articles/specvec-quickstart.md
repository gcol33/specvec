# Getting started with specvec

``` r

library(specvec)
```

## Introduction

specvec learns abundance-aware vectors for species and for whole
communities from ecological plot data. The input is the long table that
nearly every vegetation survey already produces: one row per
plot-species record, optionally carrying a cover value. From that table
specvec places every species at a point in a shared coordinate space,
where species that keep the same company across plots land near each
other. Two species that repeatedly share habitat end up with similar
vectors; two species that never co-occur sit far apart.

Cover enters the geometry through the weighting, not as an afterthought.
In the default method, a co-occurring plot contributes the geometric
mean of the two species’ covers, so a plot where both species dominate
counts for more than a plot where one is a trace. A species that is
consistently abundant in the company it keeps gets pulled harder toward
that company than a species that merely brushes past it. This is what
“abundance-aware” means here: the same word (“co-occurrence”) gets a
quantitative weight instead of a yes/no.

The second idea is that communities are first-class citizens of the same
space. Once species have vectors, a plot becomes the cover-weighted
average of the vectors of the species it holds. Plots therefore live in
the same coordinate system as species, and the distance between two
plots reflects how similar their weighted species composition is in
embedding terms. That single readout supports similarity matrices
between plots, novelty scores for a plot against a reference set, and
trajectories over time. A plot vector and a species vector can sit side
by side and be compared with the same cosine, because they were built
from the same coordinates.

It helps to contrast this with the table that usually starts a
vegetation analysis. A plot-by-species matrix records who is present and
how much, but it says nothing directly about which species behave alike:
that has to be inferred through an ordination or a distance measure
layered on top. An embedding bakes the similarity into the coordinates
themselves. Two species sit close exactly when their cover-weighted
co-occurrence is higher than chance would predict, so reading off
neighbours, clustering, and projecting all reduce to ordinary geometry
in a modest number of dimensions, typically a few dozen rather than the
thousands of raw species columns.

The shift in object is what makes the downstream work simple. A
plot-by-species matrix is wide and sparse: a regional survey easily
reaches several thousand species columns, most of them zero in any given
plot, and the rare species that fill those columns add noise more than
signal. An embedding folds that width into a few dozen dense coordinates
per species, learned from how species share plots. The coordinates are
continuous, so two species that never literally co-occur can still sit
close when they keep similar third-party company, which a raw
co-occurrence count would record as a flat zero. The same compression
carries to plots: a community becomes a short dense vector instead of a
long binary or cover row, and the cosine between two such vectors is a
direct read on compositional similarity that needs no separate distance
matrix built on top of the table.

We use specvec in three settings, and the same machinery covers all
three. The first is descriptive: which species the data treat as
ecological partners, read off as nearest neighbours. The second is
comparative across plots: which communities resemble each other, and
which ones stand apart from a chosen reference, which is the novelty
score. The third is diachronic: how a species or a community moves
through a fixed coordinate frame as the data span decades, which is the
trajectory and integration layer. Each setting reuses the species
vectors as its building block, so a single fit feeds neighbour queries,
community pooling, novelty, and trajectories without refitting.

Where the embedding pays off is in the operations it makes cheap. Once
species and plots are points, clustering is a clustering of vectors,
ordination is a projection of vectors, and a query for the most similar
plot is a nearest-neighbour search in a few dozen coordinates. None of
those steps needs the original plot-by-species table again, and none
needs a bespoke ecological distance: the geometry already encodes
association, so ordinary Euclidean and cosine measures apply. The same
vectors can also feed a downstream model that expects fixed-length
numeric inputs, which is awkward to do with a ragged species list and
trivial with a dense embedding row. specvec stops at the embeddings and
the few readouts built directly on them; what to do with the coordinates
afterward is left to the analysis at hand.

This vignette walks the full path on small simulated data: build a
`specvec_data` object, fit a species embedding, read off neighbours and
pairwise similarities, pool plots into community vectors, compare
cover-weighted against presence-weighted pooling, score novelty against
a reference block, and run the shipped method bake-off. Two short
sections then show how to inspect the return objects through their print
methods and `$` components, and why two fits of the same data are
byte-identical. Every step uses base R graphics so the example knits
anywhere, and the data are small enough (a few hundred plots) that each
fit returns in a moment. The closing section is a set of short, runnable
pointers into the deeper vignettes.

## A worked dataset

A clean way to see the geometry work is to simulate it from a known
cause. Give each of `S` species a latent niche, a point `mu` in a
two-dimensional environment. Give each of `M` plots an environmental
position `env`. A species occurs in a plot with probability that falls
off with the squared distance between the species’ niche and the plot’s
environment, and where it does occur its cover scales with that same
probability. Species with nearby niches will tend to show up together in
plots whose environment sits between them, which is the co-occurrence
structure specvec is built to recover.

The simulation is deliberately transparent so the embedding has a known
answer to match. Two species whose niches `mu` are close share many of
the same plots and both reach high cover there, which is exactly the
cover-weighted co-occurrence the default operator rewards; two species
whose niches sit far apart rarely meet. If the method works, the species
map it produces should mirror the latent niche map that generated the
data, up to rotation and scale.

A known-truth simulation buys three things that real data cannot. It
supplies an answer key: because we drew the niches `mu`, we know which
species ought to land together, so a recovered neighbour list can be
checked rather than just admired. It isolates the mechanism: the only
structure in the data is the Gaussian niche falloff, with no confounds
from sampling design, detection bias, or unmeasured environment, so when
the embedding recovers the niche layout the credit belongs to the method
and not to a lucky correlate. And it gives a controllable knob: the
falloff width, the number of plots, and the cover scaling can be turned
up or down to see when recovery holds and when it breaks. The cost is
that simulated data are cleaner than any survey, so a method that
recovers the truth here has cleared a necessary bar, not a sufficient
one. The benchmark vignette closes that gap by scoring methods on
held-out tasks rather than on a known generator.

``` r

sim_plots <- function(M = 300, S = 40, seed = 1) {
  set.seed(seed)
  mu  <- matrix(rnorm(S * 2), S, 2)
  env <- matrix(rnorm(M * 2), M, 2)
  rows <- lapply(seq_len(M), function(p) {
    d2 <- rowSums((mu - matrix(env[p, ], S, 2, byrow = TRUE))^2)
    prob <- exp(-d2 / 2)
    present <- which(runif(S) < prob)
    if (length(present) < 2) present <- order(prob, decreasing = TRUE)[1:2]
    data.frame(plot = paste0("p", p), species = paste0("sp", present),
               cover = round(100 * prob[present] / max(prob[present]), 1))
  })
  df <- do.call(rbind, rows)
  attr(df, "mu") <- mu
  df
}
df <- sim_plots()
head(df, 4)
#>   plot species cover
#> 1   p1     sp3 100.0
#> 2   p1     sp4   9.8
#> 3   p1     sp6  34.1
#> 4   p1     sp8  44.2
```

The result is a long `plot, species, cover` table, exactly the shape
[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md)
expects. The default run draws 300 plots over 40 species, and the
falloff term `exp(-d2 / 2)` keeps each plot to the handful of species
whose niches lie near its environment. A short guard promotes the two
highest-probability species into any plot that would otherwise hold
fewer than two, so no plot is empty. The cover column is the same
probability rescaled so the most likely species in each plot reads 100,
which makes the cover numbers comparable across plots of different
richness.

We pass the column roles by name and read cover on the percent scale
(values run from near 0 up to 100). Internally
[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md)
builds two sparse matrices over a sorted species index: a presence
matrix `P` of zeros and ones, and a cover matrix `COV` rescaled into
`[0, 1]`. The print method reports how many plots and species survived,
whether cover is present, and how duplicate plot-species rows were
aggregated; the `duplicates = "max"` default keeps the largest cover
when a species appears twice in one plot.

``` r

x <- specvec(df, plot = "plot", species = "species", abundance = "cover",
             cover_scale = "percent")
x
#> <specvec_data> plots=300  species=40
#>   presence: nnz=4187  density=34.8917%
#>   abundance: yes (cover_scale=percent)  duplicates=max
```

The object exposes its pieces directly. `x$species` and `x$plots` are
the sorted id maps, `x$P` and `x$COV` are the sparse matrices, and
`x$meta` records the choices above. Sorting the ids up front is what
makes a refit reproducible: the row and column order of every downstream
matrix is fixed by the sort, not by the order rows happened to arrive in
the data frame. The latent niches `mu` stay attached to `df` as an
attribute, so after fitting we can check whether species that were
assigned nearby niches really did end up as embedding neighbours, which
is the property the whole method is meant to recover.

## Species embedding

[`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)
turns the data object into species vectors. The `method` argument is
shorthand for a registered `(weighting, factorization)` pair; the
default `"abund_pmi"` is abundance-weighted pointwise mutual
information. It builds a cover-weighted co-occurrence operator, keeps
only the positive PMI entries, and takes the top eigenvectors scaled by
the square roots of their eigenvalues. The practical effect is that two
species sharing a lot of cover-heavy plots get a large positive PMI and
are pulled together, while pairs that co-occur no more than chance
contribute nothing. When a dataset carries no cover column, `abund_pmi`
falls back to the presence version and records that on the object.

Pointwise mutual information measures how much more two species co-occur
than they would if they were placed in plots independently. Raw counts
favour common species, which dominate any pair they touch simply by
being everywhere; PMI divides the observed co-occurrence by the product
of the two species’ marginals, so the score reflects association rather
than abundance alone. The “abund” part replaces the binary presence with
the square root of cover, which means a plot contributes the geometric
mean of the two covers to the pair’s score, and each species’ marginal
becomes its total cover mass instead of its plot count. Negative PMI
values, the pairs that co-occur less than chance, are clipped to zero
before the factorization, which keeps the operator sparse and the
geometry stable.

The factorization step turns that association matrix into coordinates.
The positive PMI operator is symmetric, so an eigendecomposition gives
an orthogonal set of axes ordered by how much association each explains.
Keeping the top `dim` eigenvectors and scaling each by the square root
of its eigenvalue places species in a space where the dot product
between two species vectors approximates their PMI: a pair with high
mutual information lands with a large positive inner product, a pair
with none lands near orthogonal. That is why the cosine in
[`species_similarity()`](https://gcol33.github.io/specvec/reference/species_similarity.md)
is a direct read on association, and why the first few columns, which
carry the largest eigenvalues, dominate the picture. The eigenvalues
themselves trail off, so the later columns add detail at diminishing
returns, which is the reason a modest `dim` captures most of the
recoverable structure on small data.

``` r

emb <- species_embedding(x, method = "abund_pmi", dim = 8, min_occurrence = 3)
emb
#> <specvec_embedding> method=abund_pmi  dim=8  species=40
#>   weighting=abundance_pmi  factorization=eigen
#>   kept 40/40 species (min_occurrence=3, min_cooccurrence=1)  plots=300
dim(emb$V)
#> [1] 40  8
```

`emb$V` is a species-by-dimension matrix with species as row names; here
it is eight-dimensional and holds one row per species that occurred in
at least three plots. `min_occurrence = 3` drops the rarest species
before the operator is built, since species seen once or twice give
unstable vectors, and `emb$preprocessing` records how many were kept
against the total. The full list of fitted species is in `emb$species`.
The fit is deterministic: eigenvectors are unique only up to sign, so
specvec fixes each column’s sign by its largest-magnitude entry, and two
fits of the same data return identical matrices.

The dimension is a modelling choice with a clear trade-off. A larger
`dim` keeps more of the operator’s spectrum, so finer distinctions
between species survive, at the cost of more coordinates to estimate
from the same plots and a longer tail of low-variance dimensions that
mostly carry noise on small data. A smaller `dim` compresses harder and
denoises, which helps when plots are scarce, and it speeds up every
downstream cosine. We use `dim = 8` here because 40 species over 300
plots gives little signal beyond the first several eigenvectors; on
continental data the package default of 64 is a reasonable starting
point, and the benchmark can score a few candidate dimensions side by
side.

The sticky readout is “which species sit next to this one”.
[`nearest_species()`](https://gcol33.github.io/specvec/reference/nearest_species.md)
ranks the closest species to a focal species by cosine similarity (the
default) or Euclidean distance, returning a tidy data frame.

``` r

focal <- emb$species[1]
nearest_species(emb, focal, n = 5)
#>   species similarity
#> 1    sp29  0.9200359
#> 2    sp13  0.9052377
#> 3    sp37  0.8903815
#> 4     sp6  0.8830617
#> 5    sp17  0.8536480
```

[`species_similarity()`](https://gcol33.github.io/specvec/reference/species_similarity.md)
answers the pairwise question directly. Pass two species and it returns
one cosine in `[-1, 1]`; pass one species and it returns the full named
vector of similarities to every species, handy for sorting or
thresholding.

``` r

others <- emb$species[2:3]
species_similarity(emb, focal, others[1])
#> [1] 0.4175535
sim_all <- species_similarity(emb, focal)
round(sort(sim_all, decreasing = TRUE)[1:4], 3)
#>   sp1  sp29  sp13  sp37 
#> 1.000 0.920 0.905 0.890
```

The cosine ranges over `[-1, 1]`: values near 1 mark species the
embedding treats as habitat partners, values near 0 mark species with no
shared signal, and negative values mark species that avoid each other.
The named-vector form is the basis for any downstream sort or threshold,
for instance pulling every species within a chosen similarity of a focal
taxon.

Because we know the latent niches, we can check the recovery directly
rather than trusting the neighbour list on faith. For the focal species
we line up its embedding neighbour ranking against the true distances
between latent niches `mu`, and ask whether the species the embedding
calls closest are also the ones whose niches sit nearest the focal
niche. A positive rank correlation between the two orderings means the
embedding has reconstructed the generating geometry.

``` r

mu <- attr(df, "mu")
rownames(mu) <- paste0("sp", seq_len(nrow(mu)))
mu <- mu[emb$species, ]                       # align truth to fitted species
fi <- match(focal, emb$species)
true_d <- sqrt(rowSums((mu - matrix(mu[fi, ], nrow(mu), 2, byrow = TRUE))^2))
emb_sim <- species_similarity(emb, focal)     # higher = closer in embedding
others <- setdiff(emb$species, focal)
cor(emb_sim[others], -true_d[others], method = "spearman")
#> [1] 0.8963563
```

A clearly positive correlation means the embedding neighbours of the
focal species are the species whose latent niches really are closest,
which is the recovery the method promises. The negative sign on `true_d`
lines the directions up: small latent distance should pair with high
embedding similarity. On this small run the correlation is strong rather
than perfect, because eight dimensions over a few hundred plots cannot
recover every nuance of a 40-species layout, and the cover falloff adds
sampling noise that no amount of dimension can remove. The check is
worth running on any simulation before trusting a method on data where
the truth is hidden.

The first two embedding columns carry the leading eigenvalues, so a
scatter of column 1 against column 2 gives a readable map of the
dominant gradient. The full geometry lives in all eight dimensions, and
this projection flattens the rest, but it is enough to see clusters
form. We label a handful of species so the picture stays uncluttered.

``` r

V <- emb$V
show <- emb$species[1:12]
plot(V[show, 1], V[show, 2], pch = 19, col = "#3366aa",
     xlab = "embedding dim 1", ylab = "embedding dim 2",
     main = "Species in the first two embedding dimensions")
text(V[show, 1], V[show, 2], labels = show, pos = 3, cex = 0.7)
```

![Two-dimensional scatter of the first two species-embedding columns for
a handful of species, each point labelled with its species
id.](specvec-quickstart_files/figure-html/plot-species-1.svg)

Points that cluster on this map are species the embedding judges to
share habitat, which traces back to the latent niches that generated the
data. The two-dimensional view is a projection of the eight-dimensional
fit, so two species that look close here can still differ on a later
dimension; the cosine in
[`species_similarity()`](https://gcol33.github.io/specvec/reference/species_similarity.md)
reads the full space and is the measure to trust when the scatter and
the numbers disagree.

## Community embeddings

Species vectors are the building block; community vectors are the
payoff.
[`community_embedding()`](https://gcol33.github.io/specvec/reference/community_embedding.md)
places each plot in the same space by pooling the vectors of the species
the plot contains. Passing the already-fitted embedding through
`embedding =` reuses those exact vectors instead of refitting, which
keeps the species and community spaces aligned and saves the fit. With
`weights = "cover"` (the default when cover is present) the pooling is a
cover-weighted average, so a dominant species steers the plot vector
more than a trace; with `weights = "presence"` it is the plain mean of
the present species’ vectors.

``` r

comm <- community_embedding(x, embedding = emb, weights = "cover")
comm
#> <specvec_community> plots=300  dim=8  pooling=cover
#>   from: method=abund_pmi  weighting=abundance_pmi  factorization=eigen
dim(comm$U)
#> [1] 300   8
```

`comm$U` is a plot-by-dimension matrix with plots as row names, one row
per plot that contributed species in the fitted set. The pooled-vector
readout is uniform across embedding methods, so plot vectors stay
comparable no matter which method produced the species vectors
underneath; the species space can change recipe while the community
readout stays a weighted average over it. A plot lands at the weighted
centre of its species, so a plot dominated by one cover-heavy species
sits close to that species’ vector, and a plot spread evenly across
several sits near their average.

The pooling is a weighted mean of rows of the species matrix, with the
weights taken from the plot’s cover values under `weights = "cover"` and
set to one for every present species under `weights = "presence"`.
Because the operation is linear, the plot vector inherits the geometry
of the species space: a plot whose species cluster tightly lands inside
that cluster, and a plot that mixes two distinct groups lands between
them at a point set by their relative cover. This is also why community
vectors are comparable across methods. Whatever recipe produced the
species coordinates, the plot is the same weighted average over them, so
two analyses that differ only in the species fit still pool their plots
the same way and the community distances stay interpretable on the same
footing.

Reusing the fitted embedding is the recommended pattern, and it matters
for more than speed. When the species vectors are fixed once and pooled
for every plot, the species and the communities live in literally the
same coordinates, so a species vector and a plot vector are directly
comparable and so are any two plots. If instead
[`community_embedding()`](https://gcol33.github.io/specvec/reference/community_embedding.md)
were allowed to refit (which it does when no `embedding =` is passed),
the plot vectors would still be internally consistent, but they would
belong to a fresh fit and would not align with an embedding fitted
elsewhere. The rule of thumb is to fit
[`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)
once and thread that object through every community call that should
share its frame.

The presence-weighted version answers a slightly different question. It
asks which species are there at all rather than how much of each there
is, and the two can disagree for plots dominated by one species. Cover
pooling lets a single abundant species pull the plot toward its own
corner of the space; presence pooling gives every present species an
equal vote regardless of how much ground it covers. Comparing the two
placements is a quick way to see how much abundance, rather than mere
species list, is shaping a given plot’s position.

``` r

comm_pres <- community_embedding(x, embedding = emb, weights = "presence")
comm_pres$pooling
#> [1] "presence"
```

[`community_similarity()`](https://gcol33.github.io/specvec/reference/community_similarity.md)
returns a plot-by-plot matrix of cosines. Restricting it to a few plots
makes the structure easy to read: values near 1 mark plots with
near-identical weighted composition, values near 0 mark plots with
little shared signal.

``` r

plots5 <- comm$plots[1:5]
sub <- structure(list(U = comm$U[plots5, , drop = FALSE], pooling = comm$pooling,
                      normalized = comm$normalized, from = comm$from, plots = plots5),
                 class = "specvec_community")
round(community_similarity(sub, metric = "cosine"), 2)
#>        p1  p10 p100 p101 p102
#> p1   1.00 0.69 0.67 0.85 0.62
#> p10  0.69 1.00 0.85 0.74 0.98
#> p100 0.67 0.85 1.00 0.49 0.90
#> p101 0.85 0.74 0.49 1.00 0.62
#> p102 0.62 0.98 0.90 0.62 1.00
```

[`community_novelty()`](https://gcol33.github.io/specvec/reference/community_novelty.md)
scores how unusual a plot is against a reference set: for each scored
plot it averages the distance to its `k` nearest reference communities,
so a plot that sits far from everything in the reference scores high. We
split the plots into a reference block and a query block and score the
queries against the reference.

``` r

ref_plots <- comm$plots[1:250]
qry_plots <- comm$plots[251:length(comm$plots)]
ref <- structure(list(U = comm$U[ref_plots, , drop = FALSE], pooling = comm$pooling,
                      normalized = comm$normalized, from = comm$from, plots = ref_plots),
                 class = "specvec_community")
qry <- structure(list(U = comm$U[qry_plots, , drop = FALSE], pooling = comm$pooling,
                      normalized = comm$normalized, from = comm$from, plots = qry_plots),
                 class = "specvec_community")
nov <- community_novelty(qry, reference = ref, k = 5)
round(sort(nov, decreasing = TRUE)[1:4], 3)
#>   p66   p70   p61   p88 
#> 0.277 0.250 0.226 0.207
```

The plots at the top of that list are the ones whose weighted
composition is least like anything in the reference block. Averaging
over the `k` nearest references rather than the single closest one keeps
the score from hinging on one quirky reference plot; raising `k` smooths
it further and lowers the scores of plots that sit near a small
reference cluster. On real data the reference is often a baseline period
or a reference region, and high-novelty plots flag compositions worth a
closer look, whether that is a community drifting away from its
historical state or a survey reaching habitat the reference never
sampled.

The two community sub-objects above are built by hand so the vignette
can slice `comm$U` into a reference and a query block without a second
fit. In practice the two blocks usually come from two
[`community_embedding()`](https://gcol33.github.io/specvec/reference/community_embedding.md)
calls that share one fitted `embedding =`, for example a baseline period
and a later period, which is exactly the pattern the novelty taste test
below and the communities vignette use. The hand-built version is shown
here only to make the split explicit on one simulated dataset.

## Comparing methods

specvec registers several embedding recipes, and
[`compare_embeddings()`](https://gcol33.github.io/specvec/reference/compare_embeddings.md)
scores them side by side under one protocol so the choice is grounded in
held-out performance rather than taste. It filters the species set once,
splits plots into train and test, fits each method on the training plots
through the same engine path, and scores the methods on neutral tasks.
Here we compare correspondence analysis (`ca`), presence PMI (`pmi`),
and the abundance-weighted default (`abund_pmi`) on three co-occurrence
metrics across two seeds.

``` r

bench <- compare_embeddings(x, methods = c("ca", "pmi", "abund_pmi"),
                            metrics = c("cooc_ppmi", "cooc_raw", "link_auc"),
                            dim = 8, seeds = 1:2, min_occurrence = 3)
bench
#> <specvec_benchmark> plots=300  species=40  dim=8  seeds=1,2
#>   methods: ca, pmi, abund_pmi
#> 
#>   method            cooc_ppmi         cooc_raw
#>   abund_pmi      0.714+-0.071     0.602+-0.036
#>   pmi            0.688+-0.085     0.519+-0.027
#>   ca             0.652+-0.076     0.454+-0.058
#> 
#>   verdict vs 'ca' (beats on cooc_raw & link_auc by >2 pooled SDs):
#>     reference 'ca' absent; ordering only.
```

The protocol matters as much as the scores. The species set is filtered
once, globally, before any split, so every method sees the same species
and the comparison is not confounded by one method quietly dropping
rarer taxa. The plots are then split into train and test, each method is
fit on the training plots through the same engine path that
[`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)
uses, and the held-out plots supply the scoring signal. Each seed is an
independent split, and the table reports the across-seed mean and
standard deviation so a method that wins by luck on one split is visible
as a wide spread.

The metrics read as follows. `cooc_raw` is the Spearman correlation
between embedding-space similarity and held-out raw co-occurrence
counts: it asks whether the geometry recovers how often pairs actually
appear together in plots the model never saw. `cooc_ppmi` is the same
idea scored against held-out positive PMI rather than raw counts, which
downweights the common species that dominate raw counts. `link_auc` is
link-prediction AUC: given pairs that do co-occur and never-seen
negatives that do not, it measures how well embedding similarity
separates them, with 0.5 being a coin flip and 1.0 perfect separation.
The three together cover both the strength of association the geometry
encodes and its ability to call a pair present or absent.

The verdict line applies a rule fixed before the numbers are seen: a
method beats the reference (`ca` by default) only if it exceeds it on
both neutral metrics (`cooc_raw` and `link_auc`) by more than two pooled
standard deviations. Fixing the criterion in advance keeps the
comparison honest: there is no scanning for the metric on which a
favoured method happens to look best. The same criterion is reported
whatever the outcome, so a method that merely ties, or wins on one
metric while tying the other, is labelled `partial` or `no advantage`
rather than promoted. On this small run the `link_auc` column reads as
missing, because link prediction needs more never-seen negatives than 40
species over 300 plots reliably supply, so the printed verdict falls
back to `ordering only`: the table is ranked but the two-sigma test
cannot run without both neutral metrics. The verdict sharpens as the
data grow, both because the standard deviations narrow with more seeds
and plots and because the link-prediction metric populates once the
species pool is large enough. On continental vegetation data this
protocol selects `abund_pmi`, and its lead over correspondence analysis
widens with the number of plots. The companion benchmark vignette runs
the full comparison, including the optional `glove` method, which
factorizes the same co-occurrence operator through a weighted
least-squares objective and activates when the `text2vec` package is
installed.

The numeric scores live on the object as well as in the print method.
`bench$summary` is a tidy `method x metric` frame with the mean, the
standard deviation, and the count of seeds behind each cell, ready for a
table or a plot; `bench$raw` holds the per-seed rows before aggregation;
and `bench$config` records the dimension, seeds, filter thresholds, and
reference method, so a result carries the settings that produced it.

``` r

head(bench$summary)
#>      method    metric      mean         sd n
#> 1        ca cooc_ppmi 0.6515608 0.07588473 2
#> 2        ca  cooc_raw 0.4538273 0.05835988 2
#> 3       pmi cooc_ppmi 0.6883899 0.08490871 2
#> 4       pmi  cooc_raw 0.5193675 0.02728614 2
#> 5 abund_pmi cooc_ppmi 0.7140779 0.07069192 2
#> 6 abund_pmi  cooc_raw 0.6019236 0.03561142 2
```

## Inspecting objects

Every specvec function returns an S3 object that prints a short summary
and exposes its contents through `$`. The print methods are for a quick
look; the `$` components are what downstream code reads. There is no
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html), or
[`confint()`](https://rdrr.io/r/stats/confint.html) method, because
these objects are fitted matrices with metadata rather than fitted
models with a coefficient table.

The objects are designed to be self-describing, so a fitted result
carries the information needed to interpret it without the call that
produced it. An embedding records its weighting, its factorization, and
the filter thresholds; a community records its pooling rule and the
embedding it was pooled from; a benchmark records its seeds, dimension,
and reference. This matters when a result is saved and read back later,
or passed to a colleague: the metadata travels with the matrix, so the
provenance is never lost to a forgotten script. The print methods
surface a few of these fields for a glance, and the `$` access below
reaches the rest.

A `specvec_data` object carries the sparse matrices and the id maps.
`x$P` is the presence matrix, `x$COV` is the cover matrix or `NULL`,
`x$species` and `x$plots` are the sorted ids, and `x$meta` records the
cover scale, the duplicate rule, and whether abundance was supplied.

``` r

dim(x$P)
#> [1] 300  40
x$meta$cover_scale
#> [1] "percent"
x$meta$has_abundance
#> [1] TRUE
x$species[1:5]
#> [1] "sp1"  "sp10" "sp11" "sp12" "sp13"
```

A `specvec_embedding` carries the fitted vectors and the recipe that
produced them. `emb$V` is the species-by-dimension matrix, `emb$species`
the row ids in order, `emb$method`, `emb$weighting`, and
`emb$factorization` the recipe, and `emb$preprocessing` the filter
thresholds and the kept-versus-total species counts.

``` r

dim(emb$V)
#> [1] 40  8
emb$method
#> [1] "abund_pmi"
c(emb$weighting, emb$factorization)
#> [1] "abundance_pmi" "eigen"
emb$preprocessing$n_species_kept
#> [1] 40
```

Trajectory and integration objects add an
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) method
that tidies the per-window cells into a long frame, with one row per
focal-species-by-window combination and the coordinate or distance
columns alongside the window label, center, and support count. The data
prep here is tiny so the call runs inline; the temporal and integration
vignettes use it on the real diachronic features.

``` r

set.seed(1)
dft <- data.frame(
  plot = rep(paste0("p", 1:60), each = 3),
  species = sample(c("focal", paste0("s", 1:12)), 180, replace = TRUE),
  decade = rep(c(1990, 2000, 2010), each = 60)
)
xt <- specvec(dft, "plot", "species", time = "decade")
tr <- species_trajectory(xt, species = "focal", dim = 4, min_occurrence = 1)
#> specvec: no cover in data; trajectory uses presence co-occurrence.
#> specvec: no abundance/cover in data; 'abund_pmi' falls back to presence PMI.
head(as.data.frame(tr, na.rm = TRUE), 3)
#>   species window center support        d1         d2         d3          d4
#> 1   focal   1990   1990       3 0.2453864 0.11941532 0.02738079  0.09443296
#> 2   focal   2000   2000       5 0.3259786 0.01402241 0.07018302 -0.03418499
#> 3   focal   2010   2010       7 0.2239117 0.10045860 0.12684386  0.08267929
```

## Reproducibility

Two fits of the same data return the same matrix to the bit, which is
what lets downstream comparisons be stable across sessions and machines.
Two ingredients make this hold. The species and plot ids are sorted when
the `specvec_data` object is built, so the row and column order of every
matrix is fixed by the sort rather than by the order rows arrived in the
data frame. And the eigenvectors, which are unique only up to a sign
flip, have each column’s sign fixed by its largest-magnitude entry, so a
solver that returns the opposite sign is corrected before the matrix is
returned.

``` r

emb1 <- species_embedding(x, method = "abund_pmi", dim = 8, min_occurrence = 3)
emb2 <- species_embedding(x, method = "abund_pmi", dim = 8, min_occurrence = 3)
identical(emb1$V, emb2$V)
#> [1] TRUE
```

The same determinism carries through pooling, so two community
embeddings built from the same fixed species embedding match as well.
This is why the species and community spaces stay aligned across calls,
and why a novelty score or a trajectory computed today reproduces
tomorrow without storing the fitted object.

``` r

c1 <- community_embedding(x, embedding = emb1, weights = "cover")
c2 <- community_embedding(x, embedding = emb2, weights = "cover")
identical(c1$U, c2$U)
#> [1] TRUE
```

## Where to go next

The five steps above are the whole core loop: build a data object, fit
species vectors, pool them into community vectors, score novelty,
compare methods. Each downstream topic has a dedicated vignette that
picks up where this one stops. The taste tests below are runnable
sketches, each a few lines plus a pointer to where the detail lives.

**Data prep and cover scales.**
[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md)
reads percent cover by default, and also converts proportions and
ordinal Braun-Blanquet codes, with a `duplicates` rule for repeated
plot-species rows. Ordinal scales are first-class: the Braun-Blanquet
codes below are mapped to cover-class midpoints and rescaled into
`[0, 1]` internally, and
[`cover_from_scale()`](https://gcol33.github.io/specvec/reference/cover_from_scale.md)
exposes that conversion on its own so a mapping can be checked before a
full build.

``` r

cover_from_scale(c("r", "+", "2", "5"), scale = "braun_blanquet")
#> [1] 0.001 0.005 0.150 0.875
bb <- data.frame(plot = c("a","a","b","b"), species = c("X","Y","X","Z"),
                 bb = c("3", "+", "5", "2"))
xb <- specvec(bb, "plot", "species", abundance = "bb",
              cover_scale = "braun_blanquet", duplicates = "max")
xb$meta$cover_scale
#> [1] "braun_blanquet"
```

Two preparation choices on the data object are worth knowing up front.
The `cover_scale` argument decides how the abundance column is read, and
a custom `cover_mapping` overrides an ordinal lookup when a survey uses
split classes or a local convention. The `duplicates` argument decides
what happens when the same species appears twice in one plot, which
surveys produce when layers or visits are stored as separate rows:
`"max"` keeps the largest cover, `"sum"` adds the covers, `"first"`
takes the first row, and `"error"` stops so the duplication can be
inspected. The print method reports how many duplicate pairs were
aggregated, so a surprising count flags a data-shape problem before it
reaches the embedding. The cover-scale lookups, custom `cover_mapping`,
and aggregation rules are covered in
[`vignette("specvec-data")`](https://gcol33.github.io/specvec/articles/specvec-data.md).

**Community novelty in depth.** Novelty against a reference set scales
to large reference blocks, and the interpretation shifts with the choice
of reference and `k`. The runnable sketch below scores the simulation’s
query plots against its reference plots at two values of `k` and plots
the distribution, so the smoothing effect of a larger neighbourhood is
visible.

``` r

nov5  <- community_novelty(qry, reference = ref, k = 5)
nov15 <- community_novelty(qry, reference = ref, k = 15)
hist(nov5, breaks = 12, col = "#88aacc", border = "white",
     main = "Novelty against the reference block", xlab = "novelty (k = 5)")
```

![Histogram of per-plot novelty scores against a reference
block.](specvec-quickstart_files/figure-html/taste-novelty-1.svg)

``` r

round(c(k5 = mean(nov5), k15 = mean(nov15)), 3)
#>    k5   k15 
#> 0.118 0.181
```

Raising `k` from 5 to 15 lowers the mean novelty, because averaging over
more reference neighbours pulls each score toward the bulk of the
reference and damps the influence of a single nearby plot. The right `k`
depends on how dense the reference is and how local a notion of novelty
the question wants: a small `k` flags plots that sit far from their
immediate reference neighbours, a larger `k` flags plots far from the
reference as a whole. How to choose a reference and read the novelty
scores is covered in
[`vignette("specvec-communities")`](https://gcol33.github.io/specvec/articles/specvec-communities.md).

**Tracking change over time.** Store a plot-level time column on the
data object and specvec fits one frame embedding, then places a focal
species per window at the cover-weighted centroid of the species it
co-occurs with in that fixed frame. The sketch below simulates three
decades in which a focal species drifts from one species group to
another, fits its trajectory, and plots the first frame coordinate
across windows.

``` r

set.seed(2)
grpA <- paste0("a", 1:6); grpB <- paste0("b", 1:6)
rows <- list(); pid <- 0L
for (dec in c(1990, 2000, 2010)) {
  pB <- (dec - 1990) / 20                  # focal drifts from group A to group B
  for (i in 1:40) { pid <- pid + 1L
    pool <- if (runif(1) < pB) grpB else grpA
    rows[[pid]] <- data.frame(plot = paste0("p", pid),
                              species = c(sample(pool, 3), "focal"), decade = dec) }
}
xt2 <- specvec(do.call(rbind, rows), "plot", "species", time = "decade")
tr2 <- species_trajectory(xt2, species = "focal", dim = 4,
                          weights = "presence", min_occurrence = 1)
#> specvec: no abundance/cover in data; 'abund_pmi' falls back to presence PMI.
d <- as.data.frame(tr2)
plot(d$center, d$d1, type = "b", pch = 19, xlab = "decade", ylab = "frame dim 1",
     main = "Focal species moving through a fixed frame")
```

![First frame coordinate of a focal species across three
decades.](specvec-quickstart_files/figure-html/taste-time-1.svg)

The first frame coordinate of the focal species shifts across the three
decades, tracing its move from the company of group A toward group B.
Because the frame is fitted once and every window is read out in it, the
points are directly comparable across time, with no per-window rotation
to undo. Windowed fits and trajectories for species and communities are
in
[`vignette("specvec-temporal")`](https://gcol33.github.io/specvec/articles/specvec-temporal.md).

**Alien integration.**
[`integration_trajectory()`](https://gcol33.github.io/specvec/reference/integration_trajectory.md)
reads a neophyte’s distance to the native-community centroid across
windows, in one fixed frame. A falling distance is integration: the
neophyte’s associates shift from where it arrived toward the resident
native flora. The sketch reuses the drift simulation idea, with the
focal alien moving from a ruderal pool into the native pool over three
decades, and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the distance per window.

``` r

set.seed(3)
native <- paste0("nat", 1:6); ruderal <- paste0("rud", 1:6)
rows <- list(); pid <- 0L
for (dec in c(1990, 2000, 2010)) {
  share_nat <- (dec - 1990) / 20                 # alien shifts ruderal -> native
  for (i in 1:50) { pid <- pid + 1L
    pool <- if (runif(1) < share_nat) native else ruderal
    rows[[pid]] <- data.frame(plot = paste0("p", pid),
                              species = c(sample(pool, 3), "alien"), decade = dec) }
  for (i in 1:30) { pid <- pid + 1L
    rows[[pid]] <- data.frame(plot = paste0("p", pid),
                              species = sample(native, 3), decade = dec) }
}
xi <- specvec(do.call(rbind, rows), "plot", "species", time = "decade")
it <- integration_trajectory(xi, species = "alien", native = native,
                             dim = 4, weights = "presence", min_occurrence = 1)
#> specvec: no abundance/cover in data; 'abund_pmi' falls back to presence PMI.
as.data.frame(it)[, c("window", "support", "distance")]
#>   window support    distance
#> 1   1990      50 0.949255908
#> 2   2000      50 0.493990698
#> 3   2010      50 0.002885728
```

The distance column falls across the three decades, which is the
integration signal: as the alien is recorded more often alongside native
species, its co-occurrence centroid slides toward the native-community
centroid in the fixed frame. Two choices shape the measurement. The
frame defaults to every species except the focal alien, so the alien
moves through the resident pool rather than helping define the
coordinates it is measured in, and the `native` set picks which subset
of that pool the centroid is built from. On real surveys the native flag
usually comes from a status column in the header data, and the plots are
restricted to resampled locations first so the same ground is tracked
through time. The native-set choice, the ReSurvey anchoring step, and a
worked EVA example are in
[`vignette("specvec-integration")`](https://gcol33.github.io/specvec/articles/specvec-integration.md).

**Comparing methods and the bake-off verdict.** The bake-off scores
registered methods under the fixed protocol and applies the two-sigma
verdict. The compact call below runs the default and presence PMI
against the `ca` reference on the co-occurrence metrics, and reads the
verdict structure straight off the object.

``` r

b <- compare_embeddings(x, methods = c("ca", "pmi", "abund_pmi"),
                        metrics = c("cooc_raw", "link_auc"),
                        dim = 8, seeds = 1:2, min_occurrence = 3)
b$verdict$reference
#> [1] "ca"
b$config$methods
#> [1] "ca"        "pmi"       "abund_pmi"
```

The full comparison, the EUNIS and trait metrics, and how to read each
verdict tag are in
[`vignette("specvec-benchmark")`](https://gcol33.github.io/specvec/articles/specvec-benchmark.md).

**Extending the registry.** A new method enters by registration, never
by copying a function. The sketch registers a trivial cover scale that
reads tenths (`"3"` becomes 0.3) and confirms
[`cover_from_scale()`](https://gcol33.github.io/specvec/reference/cover_from_scale.md)
picks it up, then registers a Jaccard-style weighting and a method
preset built on it, which
[`specvec_methods()`](https://gcol33.github.io/specvec/reference/specvec_methods.md)
then lists alongside the shipped five.

``` r

register_cover_scale("tenths", function(x, mapping = NULL)
  pmin(pmax(suppressWarnings(as.numeric(x)), 0), 10) / 10)
cover_from_scale(c("0", "3", "10"), scale = "tenths")
#> [1] 0.0 0.3 1.0

register_weighting("jaccard", function(data, ks, n_plots, min_cooccurrence) {
  P <- data$P[, ks, drop = FALSE]
  C <- as.matrix(Matrix::crossprod(P))               # species x species co-occurrence
  occ <- Matrix::diag(C)
  J <- C / (outer(occ, occ, "+") - C); diag(J) <- 0  # Jaccard index
  list(kind = "sym", M = Matrix::Matrix(J, sparse = TRUE), species = data$species[ks])
}, input = "species_species")
register_method("jaccard", weighting = "jaccard", factorization = "eigen")
specvec_methods()
#> [1] "abund_pmi" "ca"        "clr"       "glove"     "jaccard"   "pmi"
```

The full registry model, the operator-kind contract, and worked
weightings and factorizations are in
[`vignette("specvec-extending")`](https://gcol33.github.io/specvec/articles/specvec-extending.md).

**The maths.** The operator-to-factorization pipeline, the PMI and
AbundPMI derivations, and the compositional `clr` route are in
[`vignette("specvec-methods")`](https://gcol33.github.io/specvec/articles/specvec-methods.md).

``` r

specvec_methods()
#> [1] "abund_pmi" "ca"        "clr"       "glove"     "jaccard"   "pmi"
```
