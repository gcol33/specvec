#' @export
print.specvec_data <- function(x, ...) {
  cat(sprintf("<specvec_data> plots=%d  species=%d\n",
              length(x$plots), length(x$species)))
  nnz <- length(x$P@x)
  dens <- nnz / (as.numeric(nrow(x$P)) * ncol(x$P))
  cat(sprintf("  presence: nnz=%d  density=%.4f%%\n", nnz, 100 * dens))
  cat(sprintf("  abundance: %s (cover_scale=%s)  duplicates=%s\n",
              if (x$meta$has_abundance) "yes" else "no",
              x$meta$cover_scale, x$meta$duplicates))
  if (x$meta$n_duplicates > 0L)
    cat(sprintf("  aggregated %d duplicated plot x species pairs\n", x$meta$n_duplicates))
  if (!is.null(x$time)) {
    tv <- x$time[!is.na(x$time)]
    nd <- length(unique(tv))
    if (is.numeric(tv) && length(tv)) {
      cat(sprintf("  time: %d distinct value(s), range %s-%s\n",
                  nd, format(min(tv)), format(max(tv))))
    } else {
      cat(sprintf("  time: %d distinct value(s)\n", nd))
    }
  }
  if (!is.null(x$labels)) cat(sprintf("  labels: %s\n", paste(names(x$labels), collapse = ", ")))
  invisible(x)
}

#' @export
print.specvec_embedding <- function(x, ...) {
  tag <- x$method %||% paste(x$weighting, x$factorization, sep = "+")
  cat(sprintf("<specvec_embedding> method=%s  dim=%d  species=%d\n",
              tag, x$dim, nrow(x$V)))
  cat(sprintf("  weighting=%s  factorization=%s\n", x$weighting, x$factorization))
  pp <- x$preprocessing
  cat(sprintf("  kept %d/%d species (min_occurrence=%d, min_cooccurrence=%d)  plots=%d\n",
              pp$n_species_kept, pp$n_species_total, pp$min_occurrence,
              pp$min_cooccurrence, pp$n_plots))
  if (!is.null(pp$time_window))
    cat(sprintf("  time window: %s\n", paste(pp$time_window, collapse = "-")))
  invisible(x)
}
