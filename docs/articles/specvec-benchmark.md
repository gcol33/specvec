# Compare embedding methods on your own data

``` r

library(specvec)
```

specvec ships its bake-off as a function so you can score the embedding
methods on your own community data, not only on the data that tuned the
package default.
[`compare_embeddings()`](https://gcol33.github.io/specvec/reference/compare_embeddings.md)
fits each method on training plots, scores it on held-out plots, and
reports which method recovers the most held-out signal. The same call
also drives an alias,
[`specvec_benchmark()`](https://gcol33.github.io/specvec/reference/compare_embeddings.md),
for readers who prefer the longer name; the two are the same function.

The point of a shipped comparison is that the choice of method becomes
an empirical question on your data rather than a matter of habit.
Picking a method by taste is easy to rationalise after the fact:
correspondence analysis is the classical ordination of vegetation
science, PMI methods come with a clean information-theoretic story, and
any of them can be made to look best by quoting the metric it happens to
win. A shipped protocol removes that freedom. It fits every method on
the same training plots, scores them on the same held-out pairs, and
applies one verdict rule decided in advance, so the answer to “which
method” is a number rather than a preference.

specvec registers five recipes, and they do genuinely different things
to the same plot table. Correspondence analysis factorizes a
chi-square-residual contingency, the classical unimodal-gradient
ordination. Plain PMI and abundance-weighted PMI build a
pointwise-mutual-information operator and take its leading eigenvectors;
the abundance variant replaces presence with the square root of cover,
so a co-occurring plot contributes the geometric mean of the two covers.
GloVe fits a weighted least-squares objective over the raw co-occurrence
counts, with bias terms and a count ceiling. The compositional `clr`
route works on a centred-log-ratio covariance, which treats each plot as
a composition and removes the constant-sum constraint before measuring
how species vary together. Each is defensible on paper, and each has
regimes where it wins. The way to find out which one fits a particular
dataset is to fit them all under one protocol and read the held-out
numbers, which is what this function does.

``` r

specvec_methods()
#> [1] "abund_pmi" "ca"        "clr"       "glove"     "pmi"
```

The default method, `abund_pmi`, was chosen this way. The bake-off ran
on continental vegetation data, and abundance-weighted PMI came out
ahead on the neutral co-occurrence tasks. That decision is not binding
on your data. A checklist survey with no reliable cover, a single
habitat type with a narrow species pool, or a small regional dataset can
all shift the ordering, and the function exists so you can check rather
than assume.

## A worked dataset

The simulation below draws species presence from smooth latent niches,
so species with nearby niches co-occur and carry cover that tracks how
central the plot is to each niche. A known-truth simulation is worth the
few extra lines: because we control the niche map that generates the
data, we know what the embedding ought to recover, and we can attach the
niche coordinates as traits and the environment quadrant as a habitat
label without inventing anything. The co-occurrence metrics need only
the species table; the habitat and trait metrics need the labels and the
trait matrix, and the simulation produces both as a by-product of the
same niche structure.

``` r

sim_plots <- function(M = 600, S = 150, seed = 1, spread = 1.5) {
  set.seed(seed)
  mu  <- matrix(rnorm(S * 2) * spread, S, 2)
  env <- matrix(rnorm(M * 2) * spread, M, 2)
  rows <- lapply(seq_len(M), function(p) {
    d2   <- rowSums((mu - matrix(env[p, ], S, 2, byrow = TRUE))^2)
    prob <- exp(-d2 / 2)
    present <- which(runif(S) < prob)
    if (length(present) < 2) present <- order(prob, decreasing = TRUE)[1:2]
    quad <- paste0(if (env[p, 1] >= 0) "E" else "W",
                   if (env[p, 2] >= 0) "N" else "S")
    data.frame(plot = paste0("p", p), species = paste0("sp", present),
               cover = round(100 * prob[present] / max(prob[present]), 1),
               habitat = quad)
  })
  df <- do.call(rbind, rows)
  attr(df, "mu") <- mu
  df
}
df <- sim_plots(M = 600, S = 150)
head(df, 4)
#>   plot species cover habitat
#> 1   p1     sp8  38.2      ES
#> 2   p1     sp9  88.4      ES
#> 3   p1    sp18 100.0      ES
#> 4   p1    sp25  18.2      ES
```

Each plot draws its species from the niches whose centres sit near the
plot’s environment, and the cover of a present species scales with how
close its niche is to that environment. The `spread` factor widens the
niche centres and the plot environments together, which thins out
co-occurrence: with 150 species spread across a wider plane, many
species pairs never share a plot, and that is the regime the
link-prediction metric needs to have genuine never-seen pairs to score
against. The `habitat` column labels each plot by the quadrant of the
environment plane it fell in (`EN`, `ES`, `WN`, `WS`), which gives four
classes with enough plots each for the habitat metric. The latent niche
centres `mu` ride along on the data frame as an attribute, ready to
become a species trait table.

We pass the column roles to
[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md) and
carry the habitat column as a plot-level label so it is available to the
EUNIS metric later. Cover reads on the percent scale by default.

``` r

x <- specvec(df, plot = "plot", species = "species",
             abundance = "cover", labels = "habitat")
x
#> <specvec_data> plots=600  species=150
#>   presence: nnz=16270  density=18.0778%
#>   abundance: yes (cover_scale=percent)  duplicates=max
#>   labels: habitat
```

The print line reports 600 plots and 150 species, cover present, and the
label column carried along. The label is stored per plot (one habitat
value for each plot id), so the EUNIS metric reads a clean plot-level
class vector when it runs.

A species-by-trait matrix completes the worked dataset. Each species’
niche centre is its position on the two environment axes, so the two
columns of `mu` are exactly the kind of continuous trait the trait
metric is built to recover. We name the rows by species id so the matrix
lines up with the species the embedding learns.

``` r

mu <- attr(df, "mu")
trait_matrix <- mu
rownames(trait_matrix) <- paste0("sp", seq_len(nrow(mu)))
colnames(trait_matrix) <- c("niche_x", "niche_y")
head(trait_matrix, 3)
#>        niche_x     niche_y
#> sp1 -0.9396807  0.67528065
#> sp2  0.2754650 -0.02783975
#> sp3 -1.2534429 -0.47710256
```

Each row of `trait_matrix` is one species’ position in the
two-dimensional niche space that generated the data. If the embedding
geometry carries the niche gradient, species that sit near each other in
the embedding will share trait values, and a nearest-neighbour
regression over the embedding will predict `niche_x` and `niche_y` from
a species’ neighbours. That is the test the trait metric runs.

## Running the co-occurrence comparison

The three co-occurrence metrics need only the data object. We compare
correspondence analysis (`ca`), presence PMI (`pmi`), and the
abundance-weighted default (`abund_pmi`) across three seeds, at a small
dimension so the example fits quickly.

``` r

b <- compare_embeddings(
  x,
  methods = c("ca", "pmi", "abund_pmi"),
  metrics = c("cooc_ppmi", "cooc_raw", "link_auc"),
  dim = 16, seeds = 1:3, min_occurrence = 3
)
b
#> <specvec_benchmark> plots=600  species=150  dim=16  seeds=1,2,3
#>   methods: ca, pmi, abund_pmi
#> 
#>   method            cooc_ppmi         cooc_raw         link_auc
#>   ca             0.771+-0.003     0.690+-0.011     0.935+-0.004
#>   abund_pmi      0.818+-0.006     0.757+-0.006     0.843+-0.002
#>   pmi            0.801+-0.004     0.700+-0.006     0.816+-0.002
#> 
#>   verdict vs 'ca' (beats on cooc_raw & link_auc by >2 pooled SDs):
#>     pmi        cooc_raw+0.010[n]  link_auc-0.119[n]  -> no advantage
#>     abund_pmi  cooc_raw+0.068[Y]  link_auc-0.092[n]  -> partial
```

The print method opens with the run configuration: the number of plots
and kept species, the embedding dimension, and the seeds. Then it gives
one row per method, sorted by `link_auc`, with each cell showing the
across-seed mean and standard deviation in the form `mean+-sd`. The
closing block is the verdict, which we unpack in its own section.

Read the rows against what each metric measures. A method that scores
high on `cooc_raw` recovers the held-out raw co-occurrence pattern; a
method with `link_auc` well above 0.5 separates co-occurring pairs from
never-seen pairs. The `+-sd` part is the across-seed spread, and a wide
spread relative to the gap between two methods means the ordering is not
yet stable at this sample size.

The latent-niche simulation is the unimodal-gradient regime
correspondence analysis is built for, so `ca` tends to score well here.
That is the regime, not the verdict: the abundance-weighted advantage
that sets `abund_pmi` as the package default shows up on continental
vegetation data where co-occurrence is sparse and the species pool is
large, a setting this 600-plot toy does not reproduce.

## Reading every metric

The function reports three co-occurrence metrics, two of them neutral
and one home-field. They share a scoring spine: for every species pair,
the embedding’s gram score (the dot product of the two species vectors)
is compared against a held-out measure of how the pair actually behaves
in the test plots.

### cooc_raw: the neutral co-occurrence test

`cooc_raw` is the Spearman correlation between the embedding’s pair
scores and the raw co-occurrence counts in the held-out plots. For each
scored pair, the held-out count is the number of test plots in which
both species appear, and the embedding score is the dot product of their
vectors. A high Spearman means the geometry ranks pairs by how often
they co-occur in plots the model never saw during fitting.

This is the fair test for a PMI method, and the reason is worth stating
plainly. A PMI method optimises an association objective, the log-odds
of co-occurrence against chance, not the raw count. Scoring it against
held-out PMI would be home-field: it would reward the method for
matching the very transform it was trained to produce. Scoring against
raw counts asks something the method was not directly handed: can the
learned geometry reconstruct the plain count of shared plots. A method
that scores well on raw-count recovery has encoded co-occurrence
structure that generalises beyond its own objective, which is the
property worth selecting on.

``` r

b$summary[b$summary$metric == "cooc_raw", c("method", "mean", "sd", "n")]
#>      method      mean          sd n
#> 2        ca 0.6898363 0.010882868 3
#> 5       pmi 0.6995917 0.005676294 3
#> 8 abund_pmi 0.7574371 0.005656675 3
```

The `n` column is the number of seeds that produced a finite score for
that method and metric. With three seeds it should read 3 across the
board; a smaller number flags a method that failed on some split,
usually because the kept-species set or the dimension left the
eigensolver with nothing to return.

### link_auc: ranking real pairs above never-seen pairs

`link_auc` is a rank-based Mann-Whitney AUC, the probability that a
randomly chosen positive pair scores above a randomly chosen negative
pair. The positive and negative sets are constructed carefully so the
test cannot be gamed. A positive is a pair that co-occurs in the test
split: two species that share at least one held-out plot. A negative is
a genuinely never-seen pair: two species that co-occur in neither the
training split nor the test split. The negatives are not merely absent
from the test plots; they are absent from everything the model was
shown, so a high AUC means the embedding pushes pairs that will appear
together above pairs it has no co-occurrence evidence for at all.

``` r

b$summary[b$summary$metric == "link_auc", c("method", "mean", "sd")]
#>      method      mean          sd
#> 3        ca 0.9348572 0.003780011
#> 6       pmi 0.8162387 0.002486162
#> 9 abund_pmi 0.8430082 0.001817510
```

An AUC of 0.5 is chance, the score a coin flip would earn. Values toward
1.0 mean the gram score reliably orders pairs that co-occur above pairs
that do not. Two guards keep the metric well behaved: the positive and
negative sets each need at least 20 members for the AUC to be computed,
and each is capped (50,000 by default) so the rank computation stays
cheap on large species pools. Below the floor the metric returns `NA`
for that seed and drops out of the mean.

### cooc_ppmi: the home-field reference

`cooc_ppmi` is the same Spearman correlation as `cooc_raw`, scored
against held-out positive PMI rather than raw counts. The held-out PPMI
is computed on the test plots with the same clipping the operator uses,
so it downweights common species that dominate raw counts by sheer
ubiquity. PMI methods optimise a PPMI-style objective, so this metric is
home-field for them, and it is reported for reference rather than used
in the verdict.

``` r

b$summary[b$summary$metric == "cooc_ppmi", c("method", "mean", "sd")]
#>      method      mean          sd
#> 1        ca 0.7705946 0.002894009
#> 4       pmi 0.8014725 0.003620886
#> 7 abund_pmi 0.8175119 0.006493026
```

It is useful to see, because a large gap between a PMI method’s
`cooc_ppmi` and its `cooc_raw` tells you the method is fitting its own
objective well but generalising less to plain counts. Correspondence
analysis, which does not optimise PMI, is not advantaged by this metric,
so a method that beats `ca` on `cooc_ppmi` alone has not shown a neutral
advantage. That is exactly why the verdict ignores it.

## The protocol in depth

The scores are only as trustworthy as the protocol that produces them,
and the function fixes several things on purpose so the comparison is
fair across methods.

The species set is filtered once, globally, before any plot split. The
`min_occurrence` filter drops species occurring in fewer plots than the
threshold, and it runs on the full dataset, so every method is handed
the identical species index. No method silently drops rarer taxa on its
own and scores against an easier pair set. The kept count is recorded on
the run config.

``` r

b$config$n_species
#> [1] 150
b$config$n_plots
#> [1] 600
```

`n_species` is the count after the global `min_occurrence` filter, and
`n_plots` is the full plot count before any split. Both are fixed for
the run, so when the table lists three methods they were all fit and
scored on the same species index and the same 600 plots. A quick check
that the species set really is shared is that the summary carries one
row per method per metric with the same `n`, drawn from one filtered
index.

The plots are split 80/20 per seed: a random 20% become test plots, the
rest train. Each method is then fit on the training plots through the
same engine path
[`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)
uses, the weighting-then-factorization pipeline, so the only thing that
varies between methods is the recipe, not the plumbing. The held-out 20%
supply the co-occurrence counts the metrics score against.

Each seed is an independent split, which is where the uncertainty
estimate comes from. There is no posterior and no resampling machinery
here; the spread reported as `sd` is simply the standard deviation of a
method’s score across the seeds, each seed being a fresh train/test
partition of the same plots. Three seeds give a coarse estimate of that
spread, and more seeds tighten it. The mean across seeds is the score,
and the sd next to it is how much that score wobbles when the split
changes.

``` r

b$raw[b$raw$method == "abund_pmi", c("seed", "cooc_raw", "link_auc")]
#>   seed  cooc_raw  link_auc
#> 3    1 0.7528298 0.8440719
#> 6    2 0.7637505 0.8409096
#> 9    3 0.7557311 0.8440431
```

The raw per-seed rows are kept on the object, so you can see the
individual splits behind each mean. A method whose per-seed scores
cluster tightly has a small sd and a stable ordering; one whose scores
swing from seed to seed has a wide sd, and any ranking that depends on
it should be read with caution.

## The verdict rule

The verdict applies a rule fixed before the numbers are seen: a method
beats the reference (`ca` by default) only if it exceeds it on both
neutral metrics (`cooc_raw` and `link_auc`) by more than two pooled
standard deviations. The pooled standard deviation combines the two
methods’ across-seed spreads, so the gap has to clear roughly two sigma
of the noise on each metric, and it has to do so on both metrics at
once.

``` r

v <- b$verdict
v$reference
#> [1] "ca"
v$rows
#> $pmi
#> $pmi$method
#> [1] "pmi"
#> 
#> $pmi$beats
#> cooc_raw link_auc 
#>    FALSE    FALSE 
#> 
#> $pmi$delta
#>     cooc_raw     link_auc 
#>  0.009755387 -0.118618448 
#> 
#> $pmi$all
#> [1] FALSE
#> 
#> $pmi$any
#> [1] FALSE
#> 
#> 
#> $abund_pmi
#> $abund_pmi$method
#> [1] "abund_pmi"
#> 
#> $abund_pmi$beats
#> cooc_raw link_auc 
#>     TRUE    FALSE 
#> 
#> $abund_pmi$delta
#>    cooc_raw    link_auc 
#>  0.06760078 -0.09184900 
#> 
#> $abund_pmi$all
#> [1] FALSE
#> 
#> $abund_pmi$any
#> [1] TRUE
```

Fixing the criterion in advance is what keeps the comparison honest.
There is no scanning across metrics for the one on which a favoured
method happens to look best, and no relaxing the threshold after seeing
the result. The same rule is reported whatever the outcome, so a method
that ties, or wins one neutral metric while tying the other, is labelled
`partial` or `no advantage` rather than quietly promoted. The print tags
each non-reference method `BEATS`, `partial`, or `no advantage`, and
shows the per-metric delta with a `[Y]` or `[n]` for whether that single
metric cleared the two-sigma bar.

On this dataset `abund_pmi` clears the bar on `cooc_raw` but not on
`link_auc`, where correspondence analysis happens to lead, so it earns a
`partial` tag rather than a `BEATS`. That is the criterion working as
intended: a method that wins one neutral task and loses the other has
not shown the across-the-board advantage the rule asks for. On a few
hundred plots a clean two-sigma gap on both metrics at once is uncommon,
and `partial` or `no advantage` is the usual verdict. That is the honest
reading of a small dataset: the ordering exists but the evidence is not
strong enough to call a winner. The published continental result is
where the criterion bites. On vegetation data at 10,000 and 200,000
plots, `abund_pmi` wins both neutral metrics over correspondence
analysis, and the lead widens with the number of plots: more plots
shrink the across-seed spread and sharpen the gap, so the two-sigma bar
is cleared comfortably at continental scale even where a 600-plot sample
shows only a hint.

## Habitat recovery

Two further metrics activate when you supply labels. The first is
habitat recovery, named EUNIS after the European habitat classification
it was built for, though it works with any plot-level class label. It is
a community-embedding test: each plot is pooled into a community vector
by presence-mean pooling, a k-nearest-neighbour vote over the training
plots predicts each test plot’s class, and the prediction is scored by
macro-F1 and accuracy against the majority-class baseline.

The label is the habitat column we carried on the `specvec` object, so
we add `"eunis"` to `metrics` and point `labels` at that column.

``` r

be <- compare_embeddings(
  x,
  methods = c("ca", "pmi", "abund_pmi"),
  metrics = c("cooc_raw", "link_auc", "eunis"),
  labels  = "habitat",
  dim = 16, seeds = 1:3, min_occurrence = 3
)
be$summary[be$summary$metric %in% c("eunis_f1", "eunis_acc", "eunis_base"),
           c("method", "metric", "mean", "sd")]
#>       method     metric      mean          sd
#> 3         ca   eunis_f1 0.8982345 0.015474872
#> 4         ca  eunis_acc 0.9000000 0.014433757
#> 5         ca eunis_base 0.2388889 0.026787919
#> 8        pmi   eunis_f1 0.8826122 0.017979302
#> 9        pmi  eunis_acc 0.8861111 0.019245009
#> 10       pmi eunis_base 0.2388889 0.026787919
#> 13 abund_pmi   eunis_f1 0.8976419 0.009270240
#> 14 abund_pmi  eunis_acc 0.9000000 0.008333333
#> 15 abund_pmi eunis_base 0.2388889 0.026787919
```

The metric returns three numbers per method. `eunis_f1` is the
macro-averaged F1 across the habitat classes, which weights every class
equally regardless of how many plots it holds. `eunis_acc` is plain
accuracy, the fraction of test plots classified correctly. `eunis_base`
is the majority-class baseline, the accuracy a classifier earns by
always guessing the most common habitat. The number that matters is
`eunis_f1` against `eunis_base`: a macro-F1 well above the majority
baseline means the community geometry carries habitat structure that a
nearest-neighbour vote can read, and it does so across all four
quadrants rather than only the largest one. In the worked run every
method scores a macro-F1 near 0.90 against a 0.24 baseline, so all three
embeddings place plots from the same habitat near each other in
community space; that the methods are close on this task says the
habitat signal is strong enough that the choice of recipe barely matters
for it.

Presence-mean pooling is the deliberate choice for this metric. The
habitat test asks whether the species list of a plot, mapped into the
embedding, predicts the plot’s habitat, and presence pooling gives every
recorded species an equal vote so a dominant species cannot drag a
plot’s vector toward its own corner. Using cover pooling here would mix
two effects, the species composition and the abundance profile, and the
metric is meant to isolate composition. The pooling is the same across
methods, so the habitat scores differ only through the species geometry
each method produces, which is the comparison the metric is built to
make.

The habitat metric needs enough plots per class for the kNN vote to have
neighbours of the right label, which is why the simulation uses 600
plots over four roughly balanced quadrants. The default neighbourhood is
15, set through `control`, and a class with fewer than that many
training plots contributes little. When the labels are missing or
unusable, the metric prints a note and drops out cleanly rather than
failing the whole run.

## Trait recovery

The second label-driven metric is trait recovery. It is a
species-embedding test: a k-nearest-neighbour regression over the
species vectors predicts each held-out species’ trait value from its
neighbours, and the prediction is scored by R^2. We pass the
niche-coordinate matrix through `traits` and add `"trait"` to `metrics`.

``` r

bt <- compare_embeddings(
  x,
  methods = c("ca", "pmi", "abund_pmi"),
  metrics = c("cooc_raw", "link_auc", "trait"),
  traits  = trait_matrix,
  dim = 16, seeds = 1:3, min_occurrence = 3
)
bt$summary[bt$summary$metric == "trait_r2", c("method", "mean", "sd", "n")]
#>      method      mean          sd n
#> 3        ca 0.6907666 0.028598393 3
#> 6       pmi 0.7425845 0.027460105 3
#> 9 abund_pmi 0.7507589 0.008946875 3
```

`trait_r2` is the coefficient of determination of the kNN regression,
averaged over the trait columns. A positive R^2 means the embedding
geometry predicts a species’ niche coordinates better than the training
mean does, so species that sit near each other in the embedding really
do share trait values. That is direct evidence the species geometry
carries the niche gradient that generated the data, not just the
co-occurrence counts. An R^2 near zero means the neighbours are no more
informative than a constant, and a negative R^2 means the geometry
actively misleads the regression.

The trait metric splits species, not plots, into train and test, and it
needs roughly twenty training species with finite trait values for the
regression to have enough neighbours; below that floor it returns `NA`
for the seed. Both niche columns are continuous and defined for every
species, so the regression has full coverage here. On real data a trait
table is usually sparse, with many species lacking a measured value, and
the metric simply scores over the species that have one.

## Sensitivity to dimension and seeds

The scores depend on the run settings, and it is worth seeing how. The
embedding dimension sets how many axes the geometry has to encode
co-occurrence; too few and the operator is compressed below its real
rank, too many and later axes fit noise. We run the same comparison at
`dim = 8` and `dim = 32` and read the neutral metric.

``` r

b8  <- compare_embeddings(x, methods = c("ca", "abund_pmi"),
                          metrics = "cooc_raw", dim = 8,  seeds = 1:3,
                          min_occurrence = 3)
b32 <- compare_embeddings(x, methods = c("ca", "abund_pmi"),
                          metrics = "cooc_raw", dim = 32, seeds = 1:3,
                          min_occurrence = 3)
rbind(d8  = b8$summary$mean[b8$summary$metric == "cooc_raw"],
      d32 = b32$summary$mean[b32$summary$metric == "cooc_raw"])
#>          [,1]      [,2]
#> d8  0.6295247 0.7458313
#> d32 0.6910582 0.7574371
```

Raising the dimension usually lifts raw-count recovery up to a point,
after which extra axes add little because the leading eigenvalues
already carry the structure. The simulation has a two-dimensional niche
cause, so a modest dimension captures most of the signal and the gain
from 8 to 32 is small. On data with a richer gradient structure the
curve keeps climbing further before it flattens.

The number of seeds changes the uncertainty estimate, not the expected
score. A single seed gives a point with no spread; more seeds estimate
the across-seed sd more tightly and make the verdict’s two-sigma test
meaningful.

``` r

b1 <- compare_embeddings(x, methods = "abund_pmi", metrics = "cooc_raw",
                         dim = 16, seeds = 1, min_occurrence = 3)
b5 <- compare_embeddings(x, methods = "abund_pmi", metrics = "cooc_raw",
                         dim = 16, seeds = 1:5, min_occurrence = 3)
c(seeds1_sd = b1$summary$sd, seeds5_sd = b5$summary$sd)
#>   seeds1_sd   seeds5_sd 
#> 0.000000000 0.006195587
```

With one seed the sd is reported as 0 because there is nothing to vary
over; with five seeds it takes a real value that reflects how much the
score moves between splits. Raising `min_occurrence` is the third knob:
a higher threshold keeps only the better-sampled species, which usually
steadies the scores, while a lower threshold admits rarer species whose
unstable vectors tend to lower the raw Spearman. The filter also sets
how many species survive, and the run errors if the kept count drops to
the dimension or below.

## The GloVe method

The fifth recipe, `glove`, factorizes the same co-occurrence counts
through a weighted least-squares objective with bias terms. It is
included by listing it in `methods`, and it activates only when the
`text2vec` package is installed, since that package supplies the GloVe
fitter. When `text2vec` is absent, the method prints a note and is
skipped, and the rest of the comparison runs as usual.

``` r

bg <- compare_embeddings(
  x,
  methods = c("ca", "abund_pmi", "glove"),
  metrics = c("cooc_raw", "link_auc"),
  dim = 16, seeds = 1:3, min_occurrence = 3
)
bg$config$methods
#> [1] "ca"        "abund_pmi" "glove"
```

`bg$config$methods` lists the methods that actually ran. If `text2vec`
is installed it includes `glove`; if not, `glove` is absent and you will
have seen the skip note. The call is written the same way either way, so
a script that lists `glove` knits on a machine without `text2vec` rather
than failing, and picks up the GloVe row automatically where the package
is present.

## Practical guidance

A few rules of thumb cover most uses of the bake-off.

Use at least three seeds, and prefer five for a verdict you mean to act
on. One seed reports a score with no uncertainty and cannot trigger the
two-sigma test; three is the minimum that gives a spread, and five
sharpens it without much extra cost. The verdict’s pooled standard
deviation is only as good as the seed count behind it, so a `BEATS` tag
earned on two seeds is weaker than the same tag on five.

Set the dimension in the range 16 to 64. Below 16 the geometry is often
compressed below the real rank of the co-occurrence structure, and the
raw Spearman suffers; above 64 the later axes mostly fit noise and the
gain flattens. For a quick comparison 16 keeps the fits fast; for a
result that feeds a downstream analysis 32 or 64 is the usual choice,
matching the package default of 64.

Match `min_occurrence` to your sampling, and keep it from starving the
run. A threshold of 3 to 5 drops the rarest species whose vectors are
unstable, which usually steadies every method’s score; the package
default is 5. The filter must leave more species than the dimension, or
the run errors, so a small dataset needs either a lower dimension or a
lower threshold.

Expect the verdict to stabilise only once you have enough plots. On a
few hundred plots the across-seed spread is wide and the verdict is
usually `partial` or `no advantage`, which is the honest reading rather
than a defect. The two-sigma gap becomes reliable in the low thousands
of plots and is decisive at the tens of thousands the continental
comparison used.

When the verdict reads `partial` or `no advantage`, the conclusion is
that the methods are not separable on your data at this sample size, not
that the reference is best. A `partial` tag means one neutral metric
cleared the bar and the other did not, which often points to a real but
narrow advantage worth confirming with more seeds or more plots. A
`no advantage` tag with overlapping spreads simply means the choice of
method does not matter much here, and the cheapest or most interpretable
recipe is a reasonable default.

Do not read too much into a single-seed win. Because each seed is one
random split, a method can lead on one split and trail on another, and
the per-seed rows on the object will show this directly. A lead that
survives across seeds, with the gap exceeding the spread, is the only
kind worth acting on, which is exactly what the two-sigma verdict
encodes.

## Using your own data

[`compare_embeddings()`](https://gcol33.github.io/specvec/reference/compare_embeddings.md)
takes any `specvec` object, so the bake-off runs on a vegetation
database export. A European Vegetation Archive (EVA) species table has a
plot identifier, a taxon name, a cover value, and usually a habitat
classification column, which map onto the
[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md)
arguments directly. The plot id becomes `plot`, the taxon becomes
`species`, the cover becomes `abundance`, and the EUNIS column becomes a
label carried for the habitat metric.

``` r

library(data.table)
sp <- fread("species_export.csv")
sp <- sp[rank == "species"]                       # taxonomic scope is the caller's

x <- specvec(sp, plot = "PlotObservationID", species = "taxon",
             abundance = "cover", labels = "Eunis_lvl2")

compare_embeddings(
  x,
  methods = c("ca", "pmi", "abund_pmi", "glove"),
  metrics = c("cooc_raw", "link_auc", "eunis"),
  labels  = "Eunis_lvl2",
  dim = 64, seeds = 1:5, min_occurrence = 5
)
```

A few mapping points carry over from the worked example. The taxonomic
scope is the caller’s to set; the filter to species rank above is one
common choice, and aggregating subspecies or removing non-vascular taxa
happens before the table reaches
[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md).
The cover column should be on a consistent scale, set through
`cover_scale` if it is a proportion or a Braun-Blanquet code rather than
percent. The habitat label must be a plot-level column, one value per
plot, which an EUNIS classification already is. To add the trait metric,
build a species-by-trait matrix with row names matching the taxon ids in
the table and pass it through `traits`, exactly as the niche matrix was
passed above.

EVA is an access-controlled database, available on request rather than
as an open download, and specvec carries none of it; the comparison is
meant to run on an export you already have access to. The simulated
example here exists to make every metric runnable and to show what a
clean recovery looks like, so that when you read the numbers off your
own export you know what each column is telling you.
