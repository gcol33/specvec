# Extending specvec: weightings, factorizations, methods

``` r

library(specvec)
```

The methods that ship with SpecVec are points on a grid, and the grid is
open. A new weighting, a new factorization, or a new ordinal cover scale
enters by registration, and from that moment it behaves like a built-in.
This vignette adds four things that run: a van der Maarel cover scale, a
Jaccard association weighting with its method, a diffusion-map
factorization that reads vectors out of the symmetric operator a
different way, and a deliberate capability mismatch that the engine
rejects. The first two add operands; the third adds a readout over
operands the package already builds, which is the other axis of the
grid. The companion
[`vignette("specvec-methods")`](https://gcol33.github.io/specvec/articles/specvec-methods.md)
covers the built-in operators those examples build on.

## The registry model

A method in SpecVec is a registered name that points at a
`(weighting, factorization)` pair together with a capability descriptor:
its input shape, its native output, whether it reads cover, and a
default dimension.
[`register_method()`](https://gcol33.github.io/specvec/reference/register_method.md)
writes that record; nothing is copied. Adding the hundredth method costs
the same as adding the second, because dispatch reads the record rather
than branching on a method name. The five shipped methods are exactly
five such records. `pmi` and `abund_pmi` reuse the `eigen` factorization
with different weightings, `ca` and `glove` pair their own operands with
`svd` and the GloVe objective, and a sixth method that reused an
existing weighting would be one more
[`register_method()`](https://gcol33.github.io/specvec/reference/register_method.md)
line rather than one more code path.

``` r

specvec_methods()
#> [1] "abund_pmi" "ca"        "clr"       "glove"     "pmi"
```

Three registries sit behind those names.
[`register_weighting()`](https://gcol33.github.io/specvec/reference/register_weighting.md)
records the function that turns a `specvec_data` object into the
operator a factorization consumes.
[`register_factorization()`](https://gcol33.github.io/specvec/reference/register_factorization.md)
records the function that turns that operator into a species by
dimension matrix.
[`register_method()`](https://gcol33.github.io/specvec/reference/register_method.md)
records the pairing plus its capability flags. A weighting takes
`(data, ks, n_plots, min_cooccurrence)`, where `data$P` is the sparse
plot by species presence matrix, `data$COV` is cover or `NULL`,
`data$species` holds the character ids, and `ks` indexes the species
columns that survived the rare-species filter. It returns a list with a
`kind` field.

That `kind` field is the contract between the two halves. A weighting
declares one of three kinds, and that declaration alone decides which
factorizations can consume it. `"sym"` is a symmetric species by species
matrix `$M`, the form an eigendecomposition wants. `"counts"` is a
species by species co-occurrence matrix `$M`, the raw operand the GloVe
objective brings its own loss to. `"implicit"` carries matrix-vector
closures rather than a stored matrix, so a large contingency table is
never densified, which is how correspondence analysis runs. Alongside
`$M` (or the closures), a weighting also returns `species`, the
character ids of the kept species in the row and column order of the
operator, so the factorization can attach names to the rows of the
matrix it produces. A factorization declares, through its `accepts`
field, exactly one kind it will consume: `eigen` accepts `"sym"`, `svd`
accepts `"implicit"`, `glove` accepts `"counts"`. Pairing a weighting
with a factorization that accepts a different kind is the one error the
fitter raises before doing any numerical work, and the worked mismatch
later in this vignette shows the message it produces.

The capability descriptor on a method carries a few more flags. `input`
is `"species_species"` when the operand is a square species matrix and
`"plot_species"` when it is the plot by species contingency that
correspondence analysis factorizes implicitly. `native_output` is
`"species"` for a method that produces species vectors directly, or
`"species_and_community"` for one that yields plot vectors in the same
pass. `supports_abundance` records whether the method reads cover, which
lets a caller know an embedding will respond to the `COV` matrix rather
than only presence. `default_dim` is the dimension used when none is
given. These flags are descriptive: they let the print methods report
what a method does without running it, and they are how a method
declares its shape rather than discovering it mid-fit.

The dataset for the rest of the vignette is a small latent-niche
simulation. Each of 60 plots draws a niche center on the unit interval;
a species enters a plot with probability falling off as the squared
distance from its own niche optimum, and cover is a random percent.
Species with similar optima end up co-occurring, which is the structure
an embedding should recover. The rare-species filter stays low
(`min_occurrence` of 2) because the simulation is small. Twelve species
and sixty plots keep every operator small enough to print a corner of,
and the embedding dimension stays at six or fewer, which is more than
the latent structure needs and small enough to fit instantly.
`set.seed(1)` fixes the simulation so the numbers below reproduce.

``` r

set.seed(1)
n_plots <- 60L; n_sp <- 12L
niche <- runif(n_sp)
rows <- lapply(seq_len(n_plots), function(p) {
  prob    <- exp(-((niche - runif(1))^2) / 0.02)
  present <- which(runif(n_sp) < prob)
  if (length(present) < 2L) present <- order(prob, decreasing = TRUE)[1:2]
  data.frame(plot = sprintf("p%02d", p),
             species = sprintf("sp%02d", present),
             cover = round(runif(length(present), 1, 100)))
})
df <- do.call(rbind, rows)
x  <- specvec(df, plot = "plot", species = "species", abundance = "cover")
x
#> <specvec_data> plots=60  species=12
#>   presence: nnz=172  density=23.8889%
#>   abundance: yes (cover_scale=percent)  duplicates=max
```

## Inspecting operators with cooc_matrix()

Before writing a weighting it helps to see what one produces.
[`cooc_matrix()`](https://gcol33.github.io/specvec/reference/cooc_matrix.md)
runs a built-in weighting and hands back the operator, so the species by
species form is visible without fitting anything. The counts operator is
the plain co-occurrence matrix: entry `[a, b]` is the number of plots
holding both species, and the diagonal is each species’ occurrence
count.

``` r

Cc <- cooc_matrix(x, "counts", min_occurrence = 2)
class(Cc)
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"
dim(Cc)
#> [1] 12 12
```

The result is a `dgCMatrix`, sparse and square over the kept species.
Most species pairs never share a plot in 60 records, so most entries are
structural zeros that the sparse format never stores. The PPMI operator
runs the same co-occurrence count through a positive pointwise mutual
information transform: the count is compared against the product of the
two marginals, the logarithm of the ratio is taken, and anything below
zero is clipped to zero.

``` r

Cp <- cooc_matrix(x, "ppmi", min_occurrence = 2)
dim(Cp)
#> [1] 12 12
as.matrix(Cp)[1:5, 1:5]
#>          sp01      sp02      sp03 sp04     sp05
#> sp01 0.000000 0.5978370 0.0000000    0 1.121085
#> sp02 0.597837 0.0000000 0.6649763    0 0.000000
#> sp03 0.000000 0.6649763 0.0000000    0 0.000000
#> sp04 0.000000 0.0000000 0.0000000    0 0.000000
#> sp05 1.121085 0.0000000 0.0000000    0 0.000000
```

The dense corner shows the symmetry the eigensolver relies on:
`[sp01, sp02]` equals `[sp02, sp01]`, the diagonal is zero because a
species carries no mutual information with itself, and a pair that never
co-occurs reads as zero. That symmetric, non-negative, sparse matrix is
precisely the `kind = "sym"` operator a new weighting needs to produce
if it wants to feed `eigen`. The counts matrix above is the same idea
one transform earlier, the `kind = "counts"` operand. With those two
shapes in hand, writing a weighting becomes filling in the transform
between them.

[`cooc_matrix()`](https://gcol33.github.io/specvec/reference/cooc_matrix.md)
is also the place to confirm that `min_occurrence` and
`min_cooccurrence` do what an extension expects. `min_occurrence` is
applied before the operator is built: species seen in fewer plots than
the cutoff are dropped from `ks`, so they never reach the weighting at
all, and both operators above are square over the same surviving
species. `min_cooccurrence` is the weighting’s own business; it arrives
as the fourth argument, and a species by species weighting decides
whether to prune pairs that share too few plots. A new weighting
receives the same two controls every built-in does, which keeps the
rare-species behaviour consistent across methods rather than
reimplemented in each.

## Worked example A: a new cover scale

Vegetation surveys record abundance on ordinal scales, and the one a
survey used is rarely the three SpecVec ships with. The van der Maarel
1-9 ordinal scale maps nine class codes to percent-cover midpoints.
Registering it makes those codes convertible the same way Braun-Blanquet
codes already are. A cover-scale function takes the raw column `x` and
an optional `mapping`, and returns a numeric proportion in `[0, 1]`,
with `NA` for any value it does not recognize. The `NA` is a signal the
imputation step reads: unrecognized or missing entries are filled
downstream with the column median.

``` r

register_cover_scale("vdm", function(x, mapping = NULL) {
  m <- c(`1` = 0.1, `2` = 0.5, `3` = 1.5, `4` = 3, `5` = 8,
         `6` = 18, `7` = 37.5, `8` = 62.5, `9` = 87.5)
  key <- as.character(x)
  out <- unname(m[key])
  pmin(pmax(out, 0), 100) / 100
})
```

The body looks up each code in the named midpoint vector `m`, where
class 1 is 0.1 percent cover and class 9 is 87.5 percent.
`as.character(x)` coerces the column to the lookup keys, `m[key]`
returns the midpoints (with `NA` for any code absent from `m`),
`pmin(pmax(out, 0), 100)` clamps to the valid percent range, and
dividing by 100 returns the proportion. The scale is now usable through
[`cover_from_scale()`](https://gcol33.github.io/specvec/reference/cover_from_scale.md),
the front door every scale shares.

``` r

cover_from_scale(c("1", "5", "9"), "vdm")
#> [1] 0.001 0.080 0.875
```

Class 1 returns 0.001, class 5 returns 0.08, class 9 returns 0.875,
matching the midpoints divided by 100. To feed a survey recorded on this
scale into
[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md),
convert the column first and pass the result as proportions. The
`cover_scale` argument of
[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md) is
restricted to the three built-in names, so a freshly registered scale is
applied at the
[`cover_from_scale()`](https://gcol33.github.io/specvec/reference/cover_from_scale.md)
step rather than passed by name into
[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md).

``` r

veg <- data.frame(plot = df$plot, species = df$species,
                  vdm = sample(as.character(1:9), nrow(df), replace = TRUE))
veg$cover_prop <- cover_from_scale(veg$vdm, "vdm")
xv <- specvec(veg, plot = "plot", species = "species",
              abundance = "cover_prop", cover_scale = "proportion")
xv$meta$has_abundance
#> [1] TRUE
```

The converted column carries proportions in `[0, 1]`, so
`cover_scale = "proportion"` is the right reader for it, and the
resulting object holds cover the same as any other. A custom scale
registered once is then a one-line conversion wherever that survey’s
codes appear.

The `mapping` argument is the reason a scale is a function rather than a
fixed table. Braun-Blanquet conventions vary by author, and the split
classes (2a, 2b, 2m) are not universal, so the built-in scale accepts a
`mapping` that overrides its default midpoints. The van der Maarel scale
above ignores `mapping` because its nine classes are fixed, but a scale
whose midpoints a study revises can read the override and fall back to a
default when it is `NULL`. Returning `NA` for an unrecognized code keeps
the conversion honest: a stray code becomes a missing value the
imputation step handles, so a typo in one record stays contained to that
record.

## Worked example B: a new weighting and method

A Jaccard association is a natural weighting that SpecVec does not ship:
the fraction of plots holding either of two species that hold both. It
is symmetric and bounded in `[0, 1]`, which makes it a `kind = "sym"`
operator that `eigen` can take directly. The function below builds it.

``` r

jaccard_weighting <- function(data, ks, n_plots, min_cooccurrence) {
  P <- data$P[, ks, drop = FALSE]
  C <- Matrix::crossprod(P)                 # species x species co-occurrence counts; diag = occurrence
  f <- Matrix::diag(C)
  C <- methods::as(C, "TsparseMatrix")
  i <- C@i + 1L; j <- C@j + 1L; x <- C@x
  keep <- i != j & x > 0
  if (min_cooccurrence > 1L) keep <- keep & x >= min_cooccurrence
  i <- i[keep]; j <- j[keep]; x <- x[keep]
  jac <- x / (f[i] + f[j] - x)
  M <- Matrix::sparseMatrix(i = i, j = j, x = jac, dims = dim(C),
                            dimnames = dimnames(C))
  list(kind = "sym", M = M, species = data$species[ks])
}
```

`P` is the presence matrix restricted to the kept species.
`C <- Matrix::crossprod(P)` is the species by species co-occurrence
count, where `C[a, b]` counts the plots holding both `a` and `b` and the
diagonal `C[a, a]` is the occurrence count of `a`. The diagonal is
pulled out as `f`, so `f[a]` is the number of plots holding species `a`.
Coercing `C` to a `TsparseMatrix` exposes its stored entries as
triplets: `C@i` and `C@j` are the zero-based row and column indices,
raised to one-based with `+ 1L`, and `C@x` is the count at each. The
`keep` mask drops the diagonal (`i != j`) and any zero count, and
applies the `min_cooccurrence` cutoff only on the surviving off-diagonal
entries. The Jaccard value `jac <- x / (f[i] + f[j] - x)` divides each
co-occurrence count by the size of the union of the two species’ plot
sets, since `f[a] + f[b] - C[a, b]` counts plots holding either. Those
values are reassembled into a sparse matrix with the same dimensions and
species names as `C`, and the function returns it with `kind = "sym"`.

Two registrations wire it in.
[`register_weighting()`](https://gcol33.github.io/specvec/reference/register_weighting.md)
records the function under a name, declaring its input shape and that it
does not read cover.
[`register_method()`](https://gcol33.github.io/specvec/reference/register_method.md)
pairs that weighting with the `eigen` factorization under a friendly
method name.

``` r

register_weighting("jaccard", jaccard_weighting,
                   input = "species_species", supports_abundance = FALSE)
register_method("jaccard", weighting = "jaccard", factorization = "eigen",
                input = "species_species", native_output = "species")
"jaccard" %in% specvec_methods()
#> [1] TRUE
```

The method now lists alongside the built-ins and fits through the
ordinary front door.
[`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)
reads the method record, runs `jaccard_weighting` to build the symmetric
operator, checks that `eigen` accepts the `"sym"` kind it produced, and
takes the top eigenpairs.

``` r

emb <- species_embedding(x, method = "jaccard", dim = 6, min_occurrence = 2)
dim(emb$V)
#> [1] 12  6
nearest_species(emb, rownames(emb$V)[1], n = 3)
#>   species similarity
#> 1    sp11  0.9969326
#> 2    sp05  0.9933682
#> 3    sp12  0.9921126
```

The embedding `$V` is a 12 by 6 matrix, one row per kept species, six
dimensions as requested. `eigen` accepts this operator because its
`accepts` field is `"sym"` and the weighting declared `kind = "sym"`: a
symmetric matrix has a real eigendecomposition, and the embedding is its
leading eigenvectors scaled by the square roots of the positive
eigenvalues.
[`nearest_species()`](https://gcol33.github.io/specvec/reference/nearest_species.md)
reads the same matrix and returns the species whose vectors point most
nearly the same way, which for this simulation are the species sharing a
niche optimum.

The fitted object carries more than the matrix. `emb$method` records
`"jaccard"`, `emb$weighting` and `emb$factorization` record the pair
that produced it, and `emb$preprocessing` records the rare-species
cutoffs and the species counts, so a Jaccard embedding is
self-describing in the same way an `abund_pmi` one is. The new method
gained every downstream consumer at once: anything that takes a
`specvec_embedding` (nearest neighbours, similarities, community
pooling) works on the Jaccard vectors without a line written for the new
operand, because those functions read the matrix and its metadata rather
than the method that built it.

## The capability check

The kind contract is enforced at fit time. Pairing the Jaccard
weighting, which produces `"sym"`, with the `svd` factorization, which
accepts `"implicit"`, is a combination the engine refuses. Calling
[`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)
with `method = NULL` and an explicit `weighting`/`factorization` pair
forces exactly that mismatch.

``` r

species_embedding(x, weighting = "jaccard", factorization = "svd",
                  method = NULL, dim = 6, min_occurrence = 2)
#> Error:
#> ! factorization 'svd' consumes operator kind 'implicit', but weighting 'jaccard' produced 'sym'.
```

The fitter runs the weighting, reads the `kind` of the operator it
returned, compares it against the factorization’s `accepts`, and stops
before any decomposition when they differ. The message names both sides:
`svd` consumes `"implicit"`, but `jaccard` produced `"sym"`. The check
is cheap and total. A stored symmetric matrix has no matrix-vector
closures for `svd` to call, and an implicit operator has no `$M` for
`eigen` to decompose, so a silent type error later is converted into a
clear refusal up front. The same gate is why every built-in method pairs
a weighting with a factorization that names the matching kind.

Forcing the pair through the `weighting` and `factorization` arguments
is the explicit route a method record normally takes implicitly. When
`method` is supplied, the record names a weighting and a factorization
that already agree on kind, so the check passes by construction. Setting
`method = NULL` and naming the two halves directly is how an extension
author tries a combination the registry does not bless, and the check is
what tells them, immediately, whether the operand a weighting produces
is the operand a factorization reads. A new factorization that consumes
`"sym"` would compose with `jaccard`, `pmi`, and `abund_pmi` the moment
it registered, with no change to any of them, because the kind is the
only thing the pairing depends on. The next section makes that claim
concrete by writing one.

## Worked example C: a new factorization

The first two examples added new operands. A factorization is the other
half of the grid: the operand stays the same, and what changes is how
vectors are read out of it. The built-in `eigen` takes the symmetric
PPMI matrix and returns its leading eigenvectors scaled by the square
roots of the positive eigenvalues. That choice gives most weight to the
directions of largest association. A different and equally standard
reading normalizes the operator by the species degrees first, which
down-weights the high-occurrence species that otherwise dominate the top
eigenvectors, then takes the eigenvectors of that normalized matrix.
This is the symmetric-normalized graph Laplacian eigenmap, and it
consumes the same `kind = "sym"` operator every PMI-family weighting
already produces. We write it as a factorization and register it once.

A factorization function takes `(op, dim, ...)`, where `op` is the
operator a weighting returned and `op$M` is the symmetric matrix for a
`"sym"` operand. It returns a species by dimension matrix whose row
names are `op$species`. The engine handles the kind check before calling
it and the sign orientation after, so the function only needs to be
deterministic. The body below uses base R, `Matrix`, and `RSpectra` and
nothing else. The function `diffusion_map` here is the author’s own
code, not part of the SpecVec API.

``` r

diffusion_map <- function(op, dim, ...) {
  M <- op$M
  d <- as.numeric(Matrix::rowSums(M))            # species degree
  dinv <- ifelse(d > 0, 1 / sqrt(d), 0)
  Dn <- Matrix::Diagonal(x = dinv)
  Mn <- Dn %*% M %*% Dn                           # D^-1/2 M D^-1/2
  Mn <- methods::as((Mn + Matrix::t(Mn)) / 2, "CsparseMatrix")
  k  <- min(dim, nrow(Mn) - 1L)
  e  <- RSpectra::eigs_sym(Mn, k = k, which = "LA")
  V  <- e$vectors[, seq_len(k), drop = FALSE]
  if (k < dim) V <- cbind(V, matrix(0, nrow(V), dim - k))
  rownames(V) <- op$species
  V
}
```

`d` is the row sum of the symmetric matrix, the total association mass
each species carries, and `dinv` is its inverse square root with zero
where a species has no association at all. `Dn` is the diagonal scaling,
and `Mn = Dn %*% M %*% Dn` is the symmetric-normalized matrix
`D^{-1/2} M D^{-1/2}`. The symmetrize line removes the floating-point
asymmetry the two products introduce, so the result reaching the solver
is exactly symmetric and the fit stays deterministic. `k` is the
requested dimension capped one below the species count, the bound a
symmetric eigensolver respects.
`RSpectra::eigs_sym(Mn, k = k, which = "LA")` returns the `k` largest
algebraic eigenpairs of the normalized operator, the same solver `eigen`
uses on the raw operator. The eigenvectors become the embedding columns,
zero-padded when fewer are requested than the data supports, and the row
names carry the species ids through. Returning the eigenvectors without
an eigenvalue scaling is the deliberate difference from `eigen`: the
Laplacian eigenmap reads geometry from the normalized spectrum, while
`eigen` reads amplitude from the raw one. One registration adds it to
the factorization registry.

``` r

register_factorization("diffusion", diffusion_map,
                       kind = "matrix", accepts = "sym")
"diffusion" %in% specvec_methods()
#> [1] FALSE
```

The last line returns `FALSE`, which is correct: a factorization is
registered, but no method preset names it yet.
[`specvec_methods()`](https://gcol33.github.io/specvec/reference/specvec_methods.md)
lists method records, and the diffusion factorization is reachable
through the `factorization` argument until a
[`register_method()`](https://gcol33.github.io/specvec/reference/register_method.md)
call gives it a friendly name. `accepts = "sym"` is the entire contract.
The factorization names the one operator kind it reads, and from that
moment it pairs with every weighting that produces that kind. We fit it
against the PPMI operator by naming the two halves directly, with
`method = NULL`.

``` r

emb_d <- species_embedding(x, weighting = "ppmi", factorization = "diffusion",
                           method = NULL, dim = 6, min_occurrence = 2)
dim(emb_d$V)
#> [1] 12  6
nearest_species(emb_d, rownames(emb_d$V)[1], n = 3)
#>   species similarity
#> 1    sp11  0.5870648
#> 2    sp12  0.4854369
#> 3    sp02  0.3658322
```

The embedding is a 12 by 6 matrix, one row per kept species, and
[`nearest_species()`](https://gcol33.github.io/specvec/reference/nearest_species.md)
reads it the same way it reads any other embedding. The neighbours are
the species sharing a niche optimum, recovered through the normalized
spectrum rather than the raw one. The factorization did not name a
weighting anywhere. It declared `accepts = "sym"`, and the PPMI
weighting declared `kind = "sym"`, so the engine paired them, ran the
kind check, and applied the sign orientation, exactly as it does for the
built-in `pmi` method.

The O(1) point for the factorization half is that the same line works
against every `"sym"` weighting at once. The abundance-weighted operator
and the Jaccard operator registered earlier are both `"sym"`, so the new
factorization reads them with no edit to either.

``` r

emb_a <- species_embedding(x, weighting = "abundance_pmi",
                           factorization = "diffusion",
                           method = NULL, dim = 6, min_occurrence = 2)
emb_j <- species_embedding(x, weighting = "jaccard", factorization = "diffusion",
                           method = NULL, dim = 6, min_occurrence = 2)
c(abund = nrow(emb_a$V), jaccard = nrow(emb_j$V))
#>   abund jaccard 
#>      12      12
```

Both fits return a species by dimension matrix, and neither the
abundance weighting nor the Jaccard weighting changed by a line. Three
weightings now reach the diffusion readout through one registration, and
the four shipped `"sym"` weightings (`ppmi`, `abundance_pmi`, `clr`, and
the Jaccard added above) all reach it the same way. A method record
makes one such pairing the default for a friendly name;
[`register_method()`](https://gcol33.github.io/specvec/reference/register_method.md)
would record `("ppmi", "diffusion")` under a name of its own and give it
a capability descriptor, the subject of the next section.

## Capability descriptor flags

[`register_method()`](https://gcol33.github.io/specvec/reference/register_method.md)
records four descriptor flags alongside the weighting and factorization
names, and each one drives a specific behaviour. `input` is the operand
shape: `"species_species"` for a square species matrix, the form the
PMI-family methods build, and `"plot_species"` for the plot by species
contingency. Correspondence analysis is the one shipped method that
declares `input = "plot_species"`, because its weighting returns
matrix-vector closures over the contingency table and never forms a
square species matrix. The flag records that the operand is the
rectangular table, a fact a consumer reads without building anything.

`native_output` is `"species"` for a method that produces species
vectors and pools community vectors from them afterward, or
`"species_and_community"` for one that yields plot vectors in the same
factorization. Every shipped method is `"species"`, with the uniform
community readout pooling species vectors per plot, which keeps plot
vectors comparable across methods. `supports_abundance` records whether
the weighting reads the `COV` matrix; `abund_pmi`, `clr`, and `ca`
differ here from `pmi` and `glove`, and a caller reads the flag to know
whether an embedding will respond to cover or only to presence.
`default_dim` is the dimension a fit uses when none is passed.

These flags are read, not computed at fit time. The print method for an
embedding reports them from the method record, so a method declares its
shape once at registration rather than rediscovering it on every fit.
Registering the diffusion factorization as a named method would set
`input = "species_species"`, `native_output = "species"`, and
`supports_abundance` to match its paired weighting, and those four
values would then describe the method wherever it is listed or printed.

## Practical guidance

A few rules keep an extension clean. Add a weighting when the operand is
new: a different way to turn plots into a species by species (or plot by
species) matrix, as Jaccard turns co-occurrence into an overlap
fraction. Add a factorization when the operand is the same but the
readout differs: a different way to read vectors out of a matrix the
existing weightings already produce, as the diffusion map reads the
normalized spectrum where `eigen` reads the raw one. Add an argument to
an existing weighting or factorization when the change is a threshold or
a switch rather than a new operand or readout. Most ecological
extensions are weightings, because the modelling choice usually lives in
how association is measured, and a single factorization then serves
every weighting of its kind.

Match the kind to the factorization. A `"sym"` operator needs a
symmetric matrix and pairs with `eigen`; an asymmetric matrix breaks the
real eigendecomposition `eigen` assumes, so symmetrize it inside the
weighting before returning. A `"counts"` operator is the raw
co-occurrence matrix that an objective factorization such as `glove`
consumes with its own loss. An `"implicit"` operator carries
matrix-vector closures for a factorization that never densifies a large
operand, the route correspondence analysis takes. Declaring a kind the
paired factorization does not accept fails at the capability check, by
design.

Keep the result deterministic. The `eigen` and `svd` factorizations
return vectors that are unique only up to sign, and the engine fixes
each column’s sign so its largest-magnitude entry is positive. A
weighting that returns a stable operator therefore yields a stable
embedding: two fits of the same data return identical matrices.
Randomness inside a weighting, an unordered set used as a key, or a tie
broken by floating-point noise would leak into the vectors and break
that guarantee, so a weighting that needs a random component should seed
it. A factorization carries the same obligation. The diffusion map
symmetrizes its normalized operator before the solve precisely so the
matrix handed to
[`RSpectra::eigs_sym`](https://rdrr.io/pkg/RSpectra/man/eigs.html) is
exactly symmetric, since an asymmetric input would give a complex
spectrum and an unstable readout. The sign fix is applied after the
factorization runs, so an extension does not handle it; returning a
deterministic matrix from a deterministic operator is enough.

Registrations belong in the `.onLoad` of an extending package, where
they run once when the package loads, or in a script sourced before
fitting; the registry is shared package state, not per-call
configuration. Registering inside the function that fits would
re-register on every call, and a name registered twice silently replaces
the earlier entry, so a typo that shadows a built-in is easy to make and
quiet to miss. Keeping registrations in one place, run once, keeps the
set of available methods inspectable through
[`specvec_methods()`](https://gcol33.github.io/specvec/reference/specvec_methods.md)
at any point.

Do not register a new method for what an argument on an existing
weighting already covers. A weighting that differs from a built-in only
by a threshold or a switch belongs as a parameter to that weighting. Two
near-identical weightings split the maintenance and tempt a caller to
compare results that differ by one line of hidden code; one weighting
with an argument keeps the choice visible at the call. The Jaccard
weighting earned its own registration because the overlap fraction is a
genuinely different operand from PPMI, with its own bounded definition.

The built-in operators and the reasoning behind each are described in
[`vignette("specvec-methods")`](https://gcol33.github.io/specvec/articles/specvec-methods.md),
which walks through the cover-weighted PMI default, the compositional
`clr` route, correspondence analysis, and the GloVe objective. Reading
what those operands measure is the fastest way to decide whether a new
modelling idea is a fresh weighting, an argument on an existing one, or
a new factorization over an operand SpecVec already builds.
