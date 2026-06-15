# Nearest species

Rank the species closest to a focal species in the embedding. The sticky
demo function: ask which species an embedding places next to Robinia.

## Usage

``` r
nearest_species(emb, species, n = 10L, metric = c("cosine", "euclidean"))
```

## Arguments

- emb:

  A `specvec_embedding`.

- species:

  Focal species id (a row name of the embedding).

- n:

  Number of neighbours to return.

- metric:

  `"cosine"` (default) or `"euclidean"`.

## Value

A tidy data frame ranked nearest-first, with a `species` column and a
`similarity` column (cosine) or `distance` column (euclidean).

## Examples

``` r
df <- data.frame(plot = rep(paste0("p", 1:8), each = 2),
  species = c("A","B","A","B","A","C","B","C","A","B","B","C","A","C","A","B"))
emb <- species_embedding(specvec(df, "plot", "species"),
                         method = "pmi", dim = 3, min_occurrence = 1)
nearest_species(emb, "A", n = 2)
```
