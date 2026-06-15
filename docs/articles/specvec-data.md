# Preparing data and cover scales

``` r

library(specvec)
```

Everything in SpecVec starts from one table and one call to
[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md).
This vignette covers the table’s shape, the object the call returns, the
cover scales that turn an abundance column into a proportion, the rules
for duplicated records, the optional time and label columns, and what
happens to invalid or missing values. The companion vignette
[`vignette("specvec-methods")`](https://gcol33.github.io/specvec/articles/specvec-methods.md)
picks up from the object and defines the embedding maths.

The one call does the formatting work that vegetation analyses usually
spread across several preprocessing steps: it sorts the species and
plots into a stable index, builds the sparse presence and cover matrices
the engine reads, converts whatever cover scale the source used, and
collapses any repeated records. The sections that follow each take one
of those steps and show its inputs, its outputs, and the argument that
controls it. The aim is that by the end you can look at a vegetation
table from any source and know which arguments to set.

## Expected input

[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md)
reads a long table: one row per plot-species record. A plot column
identifies the site, a species column names the taxon, and an optional
abundance column carries cover. The same plot id repeats across its
species, and the same species id repeats across the plots that hold it.
This is the shape most vegetation databases export, where each row is
one observation of one taxon in one relevé, so in practice the export
usually needs no reshaping before it reaches
[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md).

The long form is what the co-occurrence operator wants. A wide
site-by-species table would need pivoting back to pairs of co-occurring
species, and the long table already lists them. Three plots and four
species fit on the page; the European examples later in this vignette
run to hundreds of thousands of rows with the same three columns.

``` r

df <- data.frame(
  plot    = c("p1", "p1", "p2", "p2", "p3", "p3"),
  species = c("Festuca", "Trifolium", "Festuca", "Bromus", "Festuca", "Carex"),
  cover   = c(40, 10, 60, 5, 30, 80)
)
df
#>   plot   species cover
#> 1   p1   Festuca    40
#> 2   p1 Trifolium    10
#> 3   p2   Festuca    60
#> 4   p2    Bromus     5
#> 5   p3   Festuca    30
#> 6   p3     Carex    80
```

Column roles are passed by name, so the columns can be called anything.
Here `plot` is the site id, `species` is the taxon, and `cover` is the
abundance, read on the default percent scale.

``` r

x <- specvec(df, plot = "plot", species = "species", abundance = "cover")
x
#> <specvec_data> plots=3  species=4
#>   presence: nnz=6  density=50.0000%
#>   abundance: yes (cover_scale=percent)  duplicates=max
```

The first two arguments after `data` are the plot and species column
names, and both are required. Passing them by name keeps the call
self-documenting and lets the columns carry whatever names the source
gave them, `PlotObservationID` and `taxon` as readily as `plot` and
`species`. The abundance column is optional: without it, SpecVec builds
presence only and the methods that need cover fall back to presence. The
`time` and `labels` arguments are optional too, covered further down. A
`specvec_data` object prints a short summary of its plots, species, and
settings, and the sections below open it up.

## What specvec() builds

The returned object is a list with named components. `P` and `COV` are
the two sparse matrices the engine reads; the rest are index maps and
metadata.

``` r

names(x)
#> [1] "P"       "COV"     "species" "plots"   "time"    "labels"  "meta"
x$species
#> [1] "Bromus"    "Carex"     "Festuca"   "Trifolium"
dim(x$P)
#> [1] 3 4
x$meta
#> $duplicates
#> [1] "max"
#> 
#> $cover_scale
#> [1] "percent"
#> 
#> $n_obs
#> [1] 6
#> 
#> $n_duplicates
#> [1] 0
#> 
#> $has_abundance
#> [1] TRUE
```

`x$species` is the sorted vector of unique species ids, and `x$plots` is
the sorted vector of unique plot ids. These define the row and column
order of every matrix on the object, so a species’ column position is
the same in `P`, in `COV`, and in the embedding the engine returns.
Sorting fixes that order up front, which is what makes two runs on the
same data return identical matrices. `dim(x$P)` is plots by species: 3
plots and 4 species in this example.

`P` is the presence matrix, plot by species, with a 1 where a species
occurs in a plot and an implicit 0 elsewhere. `COV` is the cover matrix,
plot by species, with each entry a cover proportion in `[0, 1]`. `COV`
is `NULL` when no abundance column is supplied. Both are stored sparse,
so a plot that holds 30 of 50,000 species stores 30 entries, not 50,000,
and an absent species costs nothing in memory. The presence methods read
`P`; the cover-weighted methods read `COV`; the same object serves both.

``` r

as.matrix(x$P)
#>    Bromus Carex Festuca Trifolium
#> p1      0     0       1         1
#> p2      1     0       1         0
#> p3      0     1       1         0
round(as.matrix(x$COV), 2)
#>    Bromus Carex Festuca Trifolium
#> p1   0.00   0.0     0.4       0.1
#> p2   0.05   0.0     0.6       0.0
#> p3   0.00   0.8     0.3       0.0
```

The percent covers from `df` divide by 100: a cover of 40 becomes 0.40,
80 becomes 0.80. `x$meta` records the settings used, including the
duplicate rule, the cover scale, the number of records read (`n_obs`),
the number of duplicated plot-species pairs found (`n_duplicates`), and
whether an abundance column was present (`has_abundance`). These fields
are worth a glance after every build: they confirm how many records
survived validation and how the cover was read, which is exactly the
kind of input check that catches a wrong-scale or wrong-column mistake
before it reaches the embedding.

## Cover scales

The `cover_scale` argument tells SpecVec how to read the abundance
column. Three scales ship: `"percent"`, `"proportion"`, and
`"braun_blanquet"`. Each maps the raw values to a cover proportion in
`[0, 1]`, which is the common currency the cover-weighted methods
expect, so the choice of scale is a reading instruction rather than a
transformation you apply yourself.
[`cover_from_scale()`](https://gcol33.github.io/specvec/reference/cover_from_scale.md)
exposes that conversion on its own, which is handy for checking how a
column will be read before committing it to a build.

`"percent"` is the default and divides by 100, with values clamped to
the `[0, 100]` range first, so a stray 120 lands at 1.0 rather than
overshooting. `"proportion"` treats the column as already on `[0, 1]`
and clamps to that range. The difference between the two is one factor
of 100, and picking the wrong one is the most common cover mistake,
which is why the standalone converter is the quick way to see the values
before they go in.

``` r

cover_from_scale(c(0, 50, 100), scale = "percent")
#> [1] 0.0 0.5 1.0
cover_from_scale(c(0, 0.5, 1), scale = "proportion")
#> [1] 0.0 0.5 1.0
```

`"braun_blanquet"` reads ordinal cover-abundance codes and replaces each
with the percent midpoint of its cover class, then divides by 100 like
the percent scale. The codes are matched case-insensitively after
trimming whitespace, so `"R"` and `" r "` both read as `r`. The shipped
codes are `r`, `+`, `1`, `2`, `2m`, `2a`, `2b`, `3`, `4`, and `5`, with
`2m`/`2a`/`2b` the split forms of class 2. The midpoints run from 0.1
percent for `r` up to 87.5 percent for `5`, which is the centre of the
75 to 100 percent class.

``` r

cover_from_scale(c("r", "+", "2", "5"), scale = "braun_blanquet")
#> [1] 0.001 0.005 0.150 0.875
```

The Braun-Blanquet midpoints are a documented convention rather than a
fixed standard: authors differ on the exact class boundaries, and the
split classes `2m`, `2a`, and `2b` are not universal. SpecVec ships one
widely used set of midpoints and makes them overridable so your analysis
can match whichever convention your data follow. The `cover_mapping`
argument replaces the lookup with your own named numeric vector of
percent midpoints. The names are the codes, the values are the percents,
and the same vector works through
[`cover_from_scale()`](https://gcol33.github.io/specvec/reference/cover_from_scale.md)
or through
[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md).
Here is a van der Maarel 1-9 ordinal scale:

``` r

vdm <- c("1" = 0.1, "2" = 0.5, "3" = 3, "4" = 10, "5" = 18,
         "6" = 37.5, "7" = 62.5, "8" = 87.5, "9" = 97.5)
cover_from_scale(c("1", "5", "9"), scale = "braun_blanquet", mapping = vdm)
#> [1] 0.001 0.180 0.975
```

Passed to
[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md),
the same mapping drives the cover column:

``` r

dv <- data.frame(
  plot = c("q1", "q1", "q2"),
  taxon = c("Festuca", "Carex", "Festuca"),
  vdm = c("9", "5", "1")
)
xv <- specvec(dv, "plot", "taxon", abundance = "vdm",
              cover_scale = "braun_blanquet", cover_mapping = vdm)
round(as.matrix(xv$COV), 3)
#>    Carex Festuca
#> q1  0.18   0.975
#> q2  0.00   0.001
```

Any percent-midpoint scale, Domin or Londo or another, fits the same
way: a named numeric vector handed to `cover_mapping`. The mechanism is
one lookup, so a new scale is a new vector, not new code. A code absent
from the mapping reads as `NA` and is filled downstream, covered in the
validation section below, so a mapping that omits a rarely used class
still builds rather than failing.

## Duplicate plot-species rows

A long table can carry the same plot-species pair more than once: two
cover estimates for one taxon in one plot, a taxonomy harmonization that
merged two names into one, or a join that left duplicates. The presence
matrix `P` is unaffected, since a species either occurs in a plot or it
does not, but the cover matrix needs a single value per pair. The
`duplicates` argument decides how those rows collapse into one cover.
The default is `"max"`, which suits the merged-name case where the
largest cover is the safest single representative.

``` r

dd <- data.frame(
  plot    = c("p1", "p1", "p1", "p2"),
  species = c("Festuca", "Festuca", "Trifolium", "Festuca"),
  cover   = c(30, 50, 10, 20)
)
```

`Festuca` in `p1` appears at 30 and 50. `"max"` keeps 0.50, `"sum"` adds
the covers to 0.80, and `"first"` keeps the row order’s first value,
0.30.

``` r

c(
  max   = as.matrix(specvec(dd, "plot", "species", abundance = "cover",
                            duplicates = "max")$COV)["p1", "Festuca"],
  sum   = as.matrix(specvec(dd, "plot", "species", abundance = "cover",
                            duplicates = "sum")$COV)["p1", "Festuca"],
  first = as.matrix(specvec(dd, "plot", "species", abundance = "cover",
                            duplicates = "first")$COV)["p1", "Festuca"]
)
#>   max   sum first 
#>   0.5   0.8   0.3
```

Which rule fits depends on what the duplicate means. `"sum"` reads two
rows as two parts of one plant’s cover and adds them; `"max"` and
`"first"` read them as competing estimates and keep one. The count of
duplicated pairs is recorded on the object regardless of the rule, so
you can see whether duplicates were present at all:

``` r

specvec(dd, "plot", "species", abundance = "cover")$meta$n_duplicates
#> [1] 1
```

`duplicates = "error"` stops with a message naming how many duplicated
pairs it found, instead of aggregating. That is the strict choice for a
dataset where a duplicate signals an upstream problem, not a real repeat
measurement: the build fails loudly so you fix the source, with no quiet
averaging over a bug.

## Time and labels

Two more columns ride along when present. `time` is a plot-level time
column, decade or year, stored as a named vector keyed by plot id.
`labels` is a character vector of plot-level column names carried as a
data frame, one row per plot. Both are read once per plot, taken from
the first record of each plot, so they belong to the plot rather than to
any single species in it.

``` r

dt <- data.frame(
  plot    = c("p1", "p1", "p2", "p2"),
  species = c("Festuca", "Bromus", "Festuca", "Carex"),
  cover   = c(40, 20, 60, 30),
  decade  = c(1980, 1980, 2010, 2010),
  habitat = c("grassland", "grassland", "fen", "fen")
)
xt <- specvec(dt, "plot", "species", abundance = "cover",
              time = "decade", labels = c("habitat"))
xt$time
#>   p1   p2 
#> 1980 2010
xt$labels
#>      habitat
#> p1 grassland
#> p2       fen
```

A stored `time` column powers the windowed embeddings and the species
and community trajectories, where the data are split into time slices
and a species’ vector is tracked across them;
[`vignette("specvec-temporal")`](https://gcol33.github.io/specvec/articles/specvec-temporal.md)
runs those. Labels are carried for the benchmark, where the EUNIS
habitat metric scores how well the embedding separates known habitat
classes in the species space;
[`vignette("specvec-benchmark")`](https://gcol33.github.io/specvec/articles/specvec-benchmark.md)
uses them. The core embedding reads neither column, so a dataset without
time or habitat information still builds and embeds, and you add these
columns only for the analyses that use them.

## Validation and missing data

Before building anything,
[`specvec()`](https://gcol33.github.io/specvec/reference/specvec.md)
drops rows with a missing or empty plot id or species id, treating both
`NA` and the empty string `""` as missing. A record with no taxon name
carries no co-occurrence and a record with no plot id cannot be placed
in a site, so neither contributes and both leave the count. If every row
turns out invalid, the build stops rather than returning an empty
object.

``` r

dn <- data.frame(
  plot    = c("p1", "p1", "p2", NA),
  species = c("Festuca", "Trifolium", "Festuca", "Carex"),
  cover   = c(40, 10, 60, 80)
)
specvec(dn, "plot", "species", abundance = "cover")$meta$n_obs
#> [1] 3
```

The table has 4 rows; one has an `NA` plot, so `n_obs` reports the 3
valid records that survived. Reading `n_obs` against the row count of
the source table is the input check this vignette keeps returning to: a
number lower than expected means rows were dropped, and the gap tells
you how many. An unrecognized Braun-Blanquet code becomes `NA` with a
warning that names the offending codes, so a typo surfaces instead of
silently mapping to a cover.

``` r

cover_from_scale(c("r", "2", "xx"), scale = "braun_blanquet")
#> Warning: braun_blanquet: 1 unrecognized code(s) set to NA (e.g. xx)
#> [1] 0.001 0.150    NA
```

Missing cover values, whether absent in the source or produced by an
unrecognized code, are filled with the median of the observed covers and
then clamped to `[0, 1]`, so `COV` has no holes for the engine to trip
on. The median is a neutral stand-in that does not pull the cover
distribution toward an extreme, and the warning still fires, so an
imputed value is visible rather than silent.

## On your own data

A continental table such as the European Vegetation Archive (EVA) maps
straight onto the same call. EVA is an access-controlled database,
available on request rather than as an open download, and specvec
bundles none of it. The archive exports one row per plot-species
observation already, so the work is naming columns: point each role at
its column and read the cover on the scale the archive uses.
[`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)
handles a multi-million-row export at a speed base R reading cannot
match, and the sparse matrices keep the built object small even when the
species count runs into the tens of thousands.

``` r

library(data.table)
eva <- fread("eva_export.csv")
x <- specvec(
  eva,
  plot      = "PlotObservationID",
  species   = "taxon",
  abundance = "cover",
  time      = "Decade",
  labels    = c("EUNIS"),
  cover_scale = "percent"
)
emb <- species_embedding(x, method = "abund_pmi", dim = 64, min_occurrence = 5)
```

From here the object feeds
[`species_embedding()`](https://gcol33.github.io/specvec/reference/species_embedding.md)
and the rest of the pipeline.
[`vignette("specvec-methods")`](https://gcol33.github.io/specvec/articles/specvec-methods.md)
defines what the embedding computes.
