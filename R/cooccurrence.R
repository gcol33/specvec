#' Co-occurrence operator (engine primitive)
#'
#' Build the weighted species x species operator a factorization would consume
#' (for `chi_square`, return the plot x species contingency, whose residual is
#' applied implicitly at factorization). Exported for power users; most callers
#' use [species_embedding()].
#'
#' @param x A `specvec_data` object.
#' @param weighting One of `"counts"`, `"ppmi"`, `"abundance_pmi"`,
#'   `"chi_square"`.
#' @param min_occurrence Drop species occurring in fewer than this many plots.
#' @param min_cooccurrence Keep a species pair only if it co-occurs in at least
#'   this many plots (species x species weightings only).
#' @param time Optional time window: restrict to plots whose stored `time` falls
#'   in the window before building the operator. Numeric time selects the closed
#'   interval `[min(time), max(time)]`; otherwise set membership. `NULL` uses all
#'   plots.
#' @return A sparse matrix: species x species for `counts`/`ppmi`/
#'   `abundance_pmi`, or plot x species for `chi_square`.
#' @export
#' @examples
#' df <- data.frame(plot = c("p1","p1","p2","p2","p3","p3"),
#'                  species = c("A","B","A","B","A","C"))
#' x <- specvec(df, "plot", "species")
#' cooc_matrix(x, "ppmi", min_occurrence = 1)
cooc_matrix <- function(x, weighting = c("counts", "ppmi", "abundance_pmi", "chi_square"),
                        min_occurrence = 5L, min_cooccurrence = 1L, time = NULL) {
  if (!inherits(x, "specvec_data")) stop("`x` must be a specvec_data object.", call. = FALSE)
  weighting <- match.arg(weighting)
  if (!is.null(time)) x <- .subset_plots(x, .time_rows(x, time))
  ks <- .kept_species(x, min_occurrence)
  w <- .get_weighting(weighting)
  op <- w$fn(x, ks, nrow(x$P), as.integer(min_cooccurrence))
  if (op$kind == "implicit") {
    out <- x$P[, ks, drop = FALSE]
    return(out)
  }
  M <- op$M
  dimnames(M) <- list(op$species, op$species)
  M
}

.kept_species <- function(x, min_occurrence) {
  occ <- as.numeric(Matrix::colSums(x$P))
  ks <- which(occ >= min_occurrence)
  if (length(ks) < 2L)
    stop(sprintf("fewer than 2 species pass min_occurrence=%d; lower the threshold.",
                 min_occurrence), call. = FALSE)
  ks
}
