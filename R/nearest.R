#' Nearest species
#'
#' Rank the species closest to a focal species in the embedding. The sticky
#' demo function: ask which species an embedding places next to Robinia.
#'
#' @param emb A `specvec_embedding`.
#' @param species Focal species id (a row name of the embedding).
#' @param n Number of neighbours to return.
#' @param metric `"cosine"` (default) or `"euclidean"`.
#' @return A tidy data frame ranked nearest-first, with a `species` column and a
#'   `similarity` column (cosine) or `distance` column (euclidean).
#' @export
#' @examples
#' df <- data.frame(plot = rep(paste0("p", 1:8), each = 2),
#'   species = c("A","B","A","B","A","C","B","C","A","B","B","C","A","C","A","B"))
#' emb <- species_embedding(specvec(df, "plot", "species"),
#'                          method = "pmi", dim = 3, min_occurrence = 1)
#' nearest_species(emb, "A", n = 2)
nearest_species <- function(emb, species, n = 10L, metric = c("cosine", "euclidean")) {
  if (!inherits(emb, "specvec_embedding")) stop("`emb` must be a specvec_embedding.", call. = FALSE)
  metric <- match.arg(metric)
  V <- emb$V
  i <- match(species, rownames(V))
  if (is.na(i)) stop(sprintf("species '%s' not in embedding.", species), call. = FALSE)
  n <- min(as.integer(n), nrow(V) - 1L)

  if (metric == "cosine") {
    Vn <- .l2_normalize_rows(V)
    score <- as.numeric(Vn %*% Vn[i, ])
    score[i] <- -Inf
    ord <- order(score, decreasing = TRUE)[seq_len(n)]
    data.frame(species = rownames(V)[ord], similarity = score[ord],
               stringsAsFactors = FALSE, row.names = NULL)
  } else {
    d <- sqrt(rowSums((V - matrix(V[i, ], nrow(V), ncol(V), byrow = TRUE))^2))
    d[i] <- Inf
    ord <- order(d, decreasing = FALSE)[seq_len(n)]
    data.frame(species = rownames(V)[ord], distance = d[ord],
               stringsAsFactors = FALSE, row.names = NULL)
  }
}

#' Species similarity
#'
#' Similarity of one species to another, or to all species.
#'
#' @param emb A `specvec_embedding`.
#' @param a Focal species id.
#' @param b Optional second species id; if `NULL`, returns the named similarity
#'   (cosine) or distance (euclidean) of `a` to every species.
#' @param metric `"cosine"` (default) or `"euclidean"`.
#' @return A scalar when `b` is supplied, otherwise a named numeric vector.
#' @export
species_similarity <- function(emb, a, b = NULL, metric = c("cosine", "euclidean")) {
  if (!inherits(emb, "specvec_embedding")) stop("`emb` must be a specvec_embedding.", call. = FALSE)
  metric <- match.arg(metric)
  V <- emb$V
  i <- match(a, rownames(V))
  if (is.na(i)) stop(sprintf("species '%s' not in embedding.", a), call. = FALSE)

  if (is.null(b)) {
    if (metric == "cosine") {
      Vn <- .l2_normalize_rows(V)
      out <- as.numeric(Vn %*% Vn[i, ])
    } else {
      out <- sqrt(rowSums((V - matrix(V[i, ], nrow(V), ncol(V), byrow = TRUE))^2))
    }
    names(out) <- rownames(V)
    return(out)
  }
  j <- match(b, rownames(V))
  if (is.na(j)) stop(sprintf("species '%s' not in embedding.", b), call. = FALSE)
  if (metric == "cosine") {
    sum(V[i, ] * V[j, ]) / sqrt(sum(V[i, ]^2) * sum(V[j, ]^2))
  } else {
    sqrt(sum((V[i, ] - V[j, ])^2))
  }
}
