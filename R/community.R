#' Community embedding
#'
#' Place each plot in the species embedding space by pooling the vectors of the
#' species it contains. The readout is uniform across methods (pooled species
#' vectors), so plot embeddings stay comparable. `weights = "cover"` pools by
#' cover (the abundance moat at plot level); `weights = "presence"` is the plain
#' mean of present species vectors.
#'
#' @param x A `specvec_data` object.
#' @param method Embedding method (see [specvec_methods()]); ignored if
#'   `embedding` is supplied.
#' @param dim Embedding dimension; ignored if `embedding` is supplied.
#' @param weights `"cover"` (default when cover is present) or `"presence"`.
#' @param embedding Optional pre-fitted `specvec_embedding` to pool, instead of
#'   fitting one here.
#' @param normalize If `TRUE`, L2-normalize each plot vector after pooling.
#' @param time Optional time window (see [species_embedding()]): pool (and, when
#'   fitting, train on) only plots in the window. `NULL` (default) uses all
#'   plots. For per-window communities in one shared frame, use
#'   [community_trajectory()].
#' @param min_occurrence,min_cooccurrence Passed to [species_embedding()] when
#'   fitting.
#' @param ... Passed to [species_embedding()].
#' @return A `specvec_community`: plot x dim matrix `U` (plot row names) plus
#'   pooling and provenance metadata.
#' @export
#' @examples
#' df <- data.frame(plot = rep(paste0("p", 1:8), each = 2),
#'   species = c("A","B","A","B","A","C","B","C","A","B","B","C","A","C","A","B"),
#'   cover = c(80,20,50,50,60,40,50,50,70,30,40,60,90,10,55,45))
#' x <- specvec(df, "plot", "species", abundance = "cover")
#' community_embedding(x, method = "abund_pmi", dim = 3, min_occurrence = 1)
community_embedding <- function(x, method = "abund_pmi", dim = 64L,
                                weights = c("cover", "presence"), embedding = NULL,
                                normalize = FALSE, time = NULL, min_occurrence = 5L,
                                min_cooccurrence = 1L, ...) {
  if (!inherits(x, "specvec_data")) stop("`x` must be a specvec_data object.", call. = FALSE)
  weights <- match.arg(weights)
  if (weights == "cover" && is.null(x$COV)) {
    message("specvec: no cover in data; community pooling uses presence.")
    weights <- "presence"
  }
  if (is.null(embedding)) {
    embedding <- species_embedding(x, method = method, dim = dim, time = time,
                                   min_occurrence = min_occurrence,
                                   min_cooccurrence = min_cooccurrence, ...)
  } else if (!inherits(embedding, "specvec_embedding")) {
    stop("`embedding` must be a specvec_embedding.", call. = FALSE)
  }

  rows <- if (is.null(time)) seq_len(nrow(x$P)) else .time_rows(x, time)
  ks <- match(embedding$species, x$species)
  V <- embedding$V
  W <- if (weights == "cover") x$COV[rows, ks, drop = FALSE] else x$P[rows, ks, drop = FALSE]
  U <- .pool_rows(W, V)
  if (normalize) U <- .l2_normalize_rows(U)
  rownames(U) <- x$plots[rows]

  structure(list(
    U = U, pooling = weights, normalized = normalize,
    from = list(method = embedding$method, dim = embedding$dim,
                weighting = embedding$weighting, factorization = embedding$factorization,
                time_window = time),
    plots = x$plots[rows]
  ), class = "specvec_community")
}

#' Community similarity
#'
#' Pairwise similarity among plot embeddings, or of each plot to a reference set.
#'
#' @param x A `specvec_community`.
#' @param reference Optional `specvec_community` or matrix to compare against;
#'   if `NULL`, compares the object's plots to themselves.
#' @param metric `"cosine"` (default) returns a similarity matrix; `"euclidean"`
#'   returns a distance matrix.
#' @return A matrix of `x` rows by `reference` rows.
#' @export
community_similarity <- function(x, reference = NULL,
                                 metric = c("cosine", "euclidean")) {
  if (!inherits(x, "specvec_community")) stop("`x` must be a specvec_community.", call. = FALSE)
  metric <- match.arg(metric)
  U <- x$U
  R <- if (is.null(reference)) U else .community_matrix(reference, "reference")
  S <- if (metric == "cosine") .cosine_sim(U, R) else .pairwise_euclid(U, R)
  rownames(S) <- rownames(U); colnames(S) <- rownames(R)
  S
}

#' Community novelty
#'
#' Per-plot novelty: the mean distance to the `k` nearest reference communities.
#' A plot far from everything in the reference set scores high.
#'
#' @param x A `specvec_community` (the plots to score).
#' @param reference A `specvec_community` or matrix of reference communities.
#' @param k Number of nearest reference communities to average over.
#' @return A named numeric vector of novelty per plot.
#' @export
community_novelty <- function(x, reference, k = 5L) {
  Q <- .community_matrix(x, "x")
  R <- .community_matrix(reference, "reference")
  k <- min(as.integer(k), nrow(R))
  if (k < 1L) stop("reference has no communities.", call. = FALSE)

  if (requireNamespace("FNN", quietly = TRUE) && nrow(R) > 200L) {
    d <- FNN::knnx.dist(R, Q, k = k)
    nov <- rowMeans(d)
  } else {
    D <- .pairwise_euclid(Q, R)
    nov <- apply(D, 1L, function(row) mean(sort(row, partial = seq_len(k))[seq_len(k)]))
  }
  names(nov) <- rownames(Q)
  nov
}

#' @export
print.specvec_community <- function(x, ...) {
  cat(sprintf("<specvec_community> plots=%d  dim=%d  pooling=%s%s\n",
              nrow(x$U), ncol(x$U), x$pooling,
              if (isTRUE(x$normalized)) " (L2-normalized)" else ""))
  cat(sprintf("  from: method=%s  weighting=%s  factorization=%s\n",
              x$from$method %||% "custom", x$from$weighting, x$from$factorization))
  invisible(x)
}
