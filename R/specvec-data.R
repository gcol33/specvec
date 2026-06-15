#' Build a specvec data object
#'
#' Turn a long `plot, species (, abundance)` table into the sparse matrices the
#' embedding engine needs. Column roles are passed by name. Abundance, time,
#' and labels are optional. A stored `time` column powers the windowed
#' embeddings and trajectories ([species_trajectory()], [community_trajectory()]).
#'
#' @param data A data frame in long form (one row per plot-species record).
#' @param plot Name of the plot/site id column.
#' @param species Name of the species column.
#' @param abundance Optional name of the cover/abundance column.
#' @param time Optional name of a plot-level time column (e.g. decade or year),
#'   stored for windowed embeddings and trajectories.
#' @param labels Optional character vector of plot-level label column names
#'   (carried for the benchmark; the core embedding never needs them).
#' @param cover_scale How to read `abundance`: `"percent"` (default),
#'   `"proportion"`, or `"braun_blanquet"`. See [cover_from_scale()].
#' @param duplicates How to aggregate duplicated `plot x species` rows:
#'   `"max"` (default), `"sum"`, `"first"`, or `"error"`.
#' @param cover_mapping Optional named numeric vector overriding an ordinal
#'   cover-scale lookup.
#' @return A `specvec_data` object: sparse `P` (presence), `COV` (cover or
#'   `NULL`), sorted `species`/`plots` id maps, optional `time`/`labels`, `meta`.
#' @export
#' @examples
#' df <- data.frame(
#'   plot = c("p1","p1","p2","p2","p3","p3"),
#'   species = c("A","B","A","B","A","C"),
#'   cover = c(40, 10, 60, 5, 30, 80)
#' )
#' specvec(df, plot = "plot", species = "species", abundance = "cover")
specvec <- function(data, plot, species, abundance = NULL, time = NULL,
                    labels = NULL,
                    cover_scale = c("percent", "proportion", "braun_blanquet"),
                    duplicates = c("max", "sum", "first", "error"),
                    cover_mapping = NULL) {
  cover_scale <- match.arg(cover_scale)
  duplicates  <- match.arg(duplicates)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  cols <- names(data)
  need <- function(nm, role) {
    if (is.null(nm)) return(invisible())
    if (!nm %in% cols) stop(sprintf("column '%s' (%s) not found in data.", nm, role),
                            call. = FALSE)
  }
  if (missing(plot) || missing(species))
    stop("both `plot` and `species` column names are required.", call. = FALSE)
  need(plot, "plot"); need(species, "species")
  need(abundance, "abundance"); need(time, "time")
  if (!is.null(labels)) for (l in labels) need(l, "label")

  dt <- data.table::as.data.table(data)
  pv <- as.character(dt[[plot]]); sv <- as.character(dt[[species]])
  ok <- !is.na(pv) & !is.na(sv) & sv != "" & pv != ""
  dt <- dt[ok]; pv <- pv[ok]; sv <- sv[ok]
  if (nrow(dt) == 0L) stop("no valid plot-species records after dropping NA/empty.", call. = FALSE)

  plots   <- sort(unique(pv))
  species_ids <- sort(unique(sv))
  pi_ <- match(pv, plots); si_ <- match(sv, species_ids)
  M <- length(plots); S <- length(species_ids)

  pa <- unique(data.table::data.table(p = pi_, s = si_))
  P <- Matrix::sparseMatrix(i = pa$p, j = pa$s, x = 1, dims = c(M, S),
                            dimnames = list(plots, species_ids))

  COV <- NULL; n_dups <- 0L; has_ab <- !is.null(abundance)
  if (has_ab) {
    cov_raw <- cover_from_scale(dt[[abundance]], scale = cover_scale, mapping = cover_mapping)
    agg <- data.table::data.table(p = pi_, s = si_, cov = as.numeric(cov_raw))
    n_dups <- agg[, .N, by = c("p", "s")][, sum(N > 1L)]
    if (duplicates == "error" && n_dups > 0L)
      stop(sprintf("%d duplicated plot x species pairs; set `duplicates=` to aggregate.",
                   n_dups), call. = FALSE)
    cagg <- switch(
      duplicates,
      max   = agg[, list(cov = suppressWarnings(max(cov, na.rm = TRUE))), by = c("p", "s")],
      sum   = agg[, list(cov = sum(cov, na.rm = TRUE)), by = c("p", "s")],
      first = agg[, list(cov = cov[1L]), by = c("p", "s")],
      error = agg[, list(cov = cov[1L]), by = c("p", "s")]
    )
    cagg[!is.finite(cov), cov := NA_real_]
    mc <- stats::median(cagg$cov, na.rm = TRUE)
    if (!is.finite(mc)) mc <- 0.01
    cagg[is.na(cov), cov := mc]
    cagg[, cov := pmin(pmax(cov, 0), 1)]
    COV <- Matrix::sparseMatrix(i = cagg$p, j = cagg$s, x = cagg$cov, dims = c(M, S),
                                dimnames = list(plots, species_ids))
  }

  time_vec <- NULL
  if (!is.null(time)) {
    idx <- match(plots, pv)
    time_vec <- dt[[time]][idx]
    names(time_vec) <- plots
  }
  label_df <- NULL
  if (!is.null(labels)) {
    idx <- match(plots, pv)
    label_df <- as.data.frame(dt[idx, labels, with = FALSE], stringsAsFactors = FALSE)
    rownames(label_df) <- plots
  }

  structure(list(
    P = P, COV = COV, species = species_ids, plots = plots,
    time = time_vec, labels = label_df,
    meta = list(duplicates = duplicates, cover_scale = cover_scale,
                n_obs = length(pv), n_duplicates = n_dups, has_abundance = has_ab)
  ), class = "specvec_data")
}

#' @rdname specvec
#' @export
as_specvec <- function(data, plot, species, abundance = NULL, time = NULL,
                       labels = NULL,
                       cover_scale = c("percent", "proportion", "braun_blanquet"),
                       duplicates = c("max", "sum", "first", "error"),
                       cover_mapping = NULL) {
  specvec(data, plot = plot, species = species, abundance = abundance, time = time,
          labels = labels, cover_scale = cover_scale, duplicates = duplicates,
          cover_mapping = cover_mapping)
}
