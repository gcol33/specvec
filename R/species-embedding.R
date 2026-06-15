#' Species embedding
#'
#' Learn species vectors from co-occurrence. `method` is sugar for a registered
#' `(weighting, factorization)` pair; pass `weighting`/`factorization` to
#' override for power use. The default `"abund_pmi"` is abundance-weighted PMI.
#'
#' @param x A `specvec_data` object.
#' @param method One of [specvec_methods()]; default `"abund_pmi"`. Pass `NULL`
#'   to drive the engine purely from `weighting` + `factorization`.
#' @param dim Embedding dimension.
#' @param weighting Optional weighting name overriding the method's.
#' @param factorization Optional factorization name overriding the method's.
#' @param time Optional time window: fit on plots whose stored `time` falls in
#'   the window. Numeric time selects the closed interval `[min(time),
#'   max(time)]`; otherwise set membership. `NULL` (default) uses all plots. For
#'   trajectories comparable across windows, use [species_trajectory()].
#' @param min_occurrence Drop species occurring in fewer than this many plots.
#' @param min_cooccurrence Keep a species pair only if it co-occurs in at least
#'   this many plots.
#' @param glove_iter Iterations for the GloVe factorizer.
#' @return A `specvec_embedding`: species x dim matrix `V` (species row names)
#'   plus method, capability, and preprocessing metadata.
#' @export
#' @examples
#' df <- data.frame(plot = rep(paste0("p", 1:6), each = 2),
#'                  species = c("A","B","A","B","A","C","B","C","A","B","B","C"))
#' x <- specvec(df, "plot", "species")
#' species_embedding(x, method = "pmi", dim = 2, min_occurrence = 1)
species_embedding <- function(x, method = "abund_pmi", dim = 64L,
                              weighting = NULL, factorization = NULL, time = NULL,
                              min_occurrence = 5L, min_cooccurrence = 1L,
                              glove_iter = 20L) {
  if (!inherits(x, "specvec_data")) stop("`x` must be a specvec_data object.", call. = FALSE)
  dim <- as.integer(dim)
  min_occurrence <- as.integer(min_occurrence)
  min_cooccurrence <- as.integer(min_cooccurrence)

  if (!is.null(method)) {
    m <- .get_method(method)
    w_name <- weighting %||% m$weighting
    f_name <- factorization %||% m$factorization
    capability <- m
  } else {
    if (is.null(weighting) || is.null(factorization))
      stop("supply either `method` or both `weighting` and `factorization`.", call. = FALSE)
    w_name <- weighting; f_name <- factorization
    capability <- NULL
  }
  xfit <- if (is.null(time)) x else .subset_plots(x, .time_rows(x, time))
  ks <- .kept_species(xfit, min_occurrence)
  n_plots <- nrow(xfit$P)
  V <- .fit_embedding(xfit, ks, w_name, f_name, dim, n_plots,
                      min_cooccurrence = min_cooccurrence, glove_iter = glove_iter)

  structure(list(
    V = V,
    species = rownames(V),
    method = method,
    weighting = w_name,
    factorization = f_name,
    dim = dim,
    capability = capability,
    preprocessing = list(min_occurrence = min_occurrence,
                         min_cooccurrence = min_cooccurrence,
                         n_species_kept = length(ks),
                         n_species_total = ncol(x$P),
                         n_plots = n_plots,
                         time_window = time),
    spectrum = NULL,
    call = match.call()
  ), class = "specvec_embedding")
}

## Build operator -> factorize -> sign-orient, for a fixed kept-species index
## `ks` over `data` (a specvec_data or any list with $P/$COV/$species). Shared by
## species_embedding (ks from the global rare-species filter) and the benchmark
## (ks fixed globally, data restricted to training plots), so the two never copy
## the weighting/factorization plumbing. Returns a species x dim matrix in `ks`
## order with species row names.
.fit_embedding <- function(data, ks, w_name, f_name, dim, n_plots,
                           min_cooccurrence = 1L, glove_iter = 20L) {
  w <- .get_weighting(w_name)
  f <- .get_factorization(f_name)
  op <- w$fn(data, ks, n_plots, min_cooccurrence)
  if (!identical(op$kind, f$accepts))
    stop(sprintf("factorization '%s' consumes operator kind '%s', but weighting '%s' produced '%s'.",
                 f_name, f$accepts, w_name, op$kind), call. = FALSE)
  V <- f$fn(op, dim, glove_iter = glove_iter)
  V <- .sign_orient(V)
  rownames(V) <- op$species
  V
}
