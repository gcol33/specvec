## Temporal layer (v0.2). Two windowing concepts share these internals:
##   .time_rows()    : one window (a slice filter) -> plot row indices. Used by
##                     cooc_matrix(), species_embedding(), community_embedding().
##   .time_windows() : many windows (split the time axis) -> per-window row sets.
##                     Used by species_trajectory(), community_trajectory().
## Both read the plot-level `time` vector stored by specvec(). No new maths: the
## frame embedding reuses species_embedding(), pooling reuses .pool_rows(), and
## over-time novelty reuses community_novelty().

## Plot row indices for a single time window. Numeric time -> closed interval
## [min(time), max(time)] (the "windowed" reading; a single value selects that
## value exactly). Non-numeric time -> set membership. NA times are excluded.
.time_rows <- function(x, time) {
  if (is.null(x$time))
    stop("this specvec_data has no `time`; rebuild with specvec(..., time = ).",
         call. = FALSE)
  t <- x$time
  ok <- !is.na(t)
  if (is.numeric(t) && is.numeric(time)) {
    lo <- min(time); hi <- max(time)
    keep <- ok & t >= lo & t <= hi
  } else {
    keep <- ok & (as.character(t) %in% as.character(time))
  }
  rows <- which(keep)
  if (length(rows) < 2L)
    stop(sprintf("time window selects %d plot(s); need at least 2.", length(rows)),
         call. = FALSE)
  rows
}

## Split the time axis into ordered windows. `by = NULL`: one window per distinct
## time value (natural for decade-coded data). `by =` a numeric break vector:
## bin with cut(include.lowest). Returns aligned lists of plot-row vectors,
## labels, and numeric centers (for ordering and plotting); empty windows drop.
.time_windows <- function(x, by = NULL) {
  if (is.null(x$time))
    stop("this specvec_data has no `time`; rebuild with specvec(..., time = ).",
         call. = FALSE)
  t <- x$time
  ok <- !is.na(t)
  if (is.null(by)) {
    vals <- sort(unique(t[ok]))
    rows <- lapply(vals, function(v) which(ok & t == v))
    labels  <- as.character(vals)
    centers <- if (is.numeric(vals)) as.numeric(vals) else seq_along(vals)
  } else {
    by <- sort(unique(as.numeric(by)))
    if (length(by) < 2L) stop("`by` needs at least two breaks.", call. = FALSE)
    g <- cut(as.numeric(t), breaks = by, include.lowest = TRUE)
    labels  <- levels(g)
    rows <- lapply(labels, function(l) which(ok & as.character(g) == l))
    centers <- (utils::head(by, -1L) + utils::tail(by, -1L)) / 2
  }
  keep <- vapply(rows, length, integer(1)) > 0L
  if (!any(keep)) stop("no non-empty time windows.", call. = FALSE)
  list(rows = rows[keep], labels = labels[keep], centers = centers[keep])
}

## A specvec_data restricted to a set of species columns, keeping every plot row
## so the time/labels stay aligned. The column analogue of .subset_plots(); used
## to fit a frame embedding over a chosen species subset through species_embedding().
.subset_species <- function(x, cols) {
  structure(list(
    P = x$P[, cols, drop = FALSE],
    COV = if (is.null(x$COV)) NULL else x$COV[, cols, drop = FALSE],
    species = x$species[cols], plots = x$plots,
    time = x$time, labels = x$labels, meta = x$meta
  ), class = "specvec_data")
}

#' Species trajectory through a fixed embedding frame
#'
#' Track focal species through time without the cross-time alignment problem.
#' A single *frame* embedding is fitted once (the stable coordinate system), and
#' each focal species is then placed, per time window, at the cover-weighted
#' centroid of the frame species it co-occurs with in that window. Because every
#' window is read out in one fixed frame, the trajectory points are directly
#' comparable across time -- no per-window rotation, no Procrustes alignment, no
#' circularity. This is the engine behind alien-integration trajectories (v0.3).
#'
#' The frame defaults to every species except the focal ones, so the focal
#' species move through the background community rather than helping define it.
#' Pass `frame` to fix the frame to a reference set (e.g. native species only).
#'
#' @param x A `specvec_data` object built with a `time` column.
#' @param species Character vector of focal species ids to trace.
#' @param frame Optional character vector of species defining the fixed frame;
#'   default is every species except `species`.
#' @param by Window definition passed to the time splitter: `NULL` (default) for
#'   one window per distinct time value, or a numeric break vector to bin time.
#' @param method Frame embedding method (see [specvec_methods()]).
#' @param dim Frame embedding dimension.
#' @param weights `"cover"` (default when cover is present) weights co-occurrence
#'   by the geometric mean of covers, matching AbundPMI; `"presence"` counts
#'   co-occurring plots.
#' @param min_occurrence,min_cooccurrence Frame species/pair filters.
#' @param frame_embedding Optional pre-fitted `specvec_embedding` to reuse as the
#'   frame instead of fitting one (its species become the frame).
#' @param ... Passed to [species_embedding()] when fitting the frame.
#' @return A `specvec_trajectory`: `U` (focal x window x dim array, `NA` where a
#'   focal species has no co-occurrence in a window), `support` (focal x window
#'   count of window plots containing the focal species), the fixed `frame`
#'   embedding, the `windows` table, and pooling provenance.
#' @seealso [community_trajectory()], [species_embedding()]
#' @export
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   plot = rep(paste0("p", 1:60), each = 3),
#'   species = sample(c("focal", paste0("s", 1:12)), 180, replace = TRUE),
#'   decade = rep(c(1990, 2000, 2010), each = 60)
#' )
#' x <- specvec(df, "plot", "species", time = "decade")
#' tr <- species_trajectory(x, species = "focal", dim = 4, min_occurrence = 1)
#' as.data.frame(tr)
species_trajectory <- function(x, species, frame = NULL, by = NULL,
                               method = "abund_pmi", dim = 64L,
                               weights = c("cover", "presence"),
                               min_occurrence = 5L, min_cooccurrence = 1L,
                               frame_embedding = NULL, ...) {
  if (!inherits(x, "specvec_data")) stop("`x` must be a specvec_data object.", call. = FALSE)
  weights <- match.arg(weights)
  species <- as.character(species)
  if (length(species) < 1L) stop("`species` must name at least one focal species.", call. = FALSE)
  miss <- setdiff(species, x$species)
  if (length(miss))
    stop(sprintf("focal species not in data: %s", paste(utils::head(miss, 5L), collapse = ", ")),
         call. = FALSE)
  if (weights == "cover" && is.null(x$COV)) {
    message("specvec: no cover in data; trajectory uses presence co-occurrence.")
    weights <- "presence"
  }

  if (is.null(frame_embedding)) {
    frame_species <- if (is.null(frame)) setdiff(x$species, species) else as.character(frame)
    cols <- which(x$species %in% frame_species)
    if (length(cols) < 2L) stop("the frame has fewer than 2 species.", call. = FALSE)
    frame_embedding <- species_embedding(.subset_species(x, cols), method = method, dim = dim,
                                         min_occurrence = min_occurrence,
                                         min_cooccurrence = min_cooccurrence, ...)
  } else if (!inherits(frame_embedding, "specvec_embedding")) {
    stop("`frame_embedding` must be a specvec_embedding.", call. = FALSE)
  }

  Vf <- frame_embedding$V
  fr_idx <- match(rownames(Vf), x$species)
  focal_idx <- match(species, x$species)
  use_cover <- weights == "cover" && !is.null(x$COV)
  Mfull <- if (use_cover) x$COV else x$P

  win <- .time_windows(x, by)
  nF <- length(species); nW <- length(win$rows); D <- ncol(Vf)
  U <- array(NA_real_, dim = c(nF, nW, D), dimnames = list(species, win$labels, NULL))
  support <- matrix(0L, nF, nW, dimnames = list(species, win$labels))
  self_in_frame <- match(species, rownames(Vf))   # NA when focal not in frame

  for (w in seq_len(nW)) {
    rows <- win$rows[[w]]
    Wfoc <- Mfull[rows, focal_idx, drop = FALSE]    # window plots x focal
    Wfr  <- Mfull[rows, fr_idx,    drop = FALSE]    # window plots x frame
    if (use_cover) { Wfoc <- sqrt(Wfoc); Wfr <- sqrt(Wfr) }
    A <- as.matrix(Matrix::crossprod(Wfr, Wfoc))    # frame x focal co-occurrence
    present <- as.numeric(Matrix::colSums(Mfull[rows, focal_idx, drop = FALSE] > 0))
    for (f in seq_len(nF)) {
      support[f, w] <- as.integer(present[f])
      a <- A[, f]
      if (!is.na(self_in_frame[f])) a[self_in_frame[f]] <- 0   # never self-place
      tot <- sum(a)
      if (tot > 0) U[f, w, ] <- as.numeric(a %*% Vf) / tot
    }
  }

  structure(list(
    U = U, support = support, frame = frame_embedding,
    windows = data.frame(label = win$labels, center = win$centers,
                         n_plots = vapply(win$rows, length, integer(1)),
                         stringsAsFactors = FALSE),
    species = species, weights = weights,
    from = list(method = frame_embedding$method %||% method, dim = D,
                frame_species = nrow(Vf))
  ), class = "specvec_trajectory")
}

#' @export
print.specvec_trajectory <- function(x, ...) {
  cat(sprintf("<specvec_trajectory> focal=%d  windows=%d  dim=%d  weights=%s\n",
              length(x$species), nrow(x$windows), dim(x$U)[3L], x$weights))
  cat(sprintf("  frame: method=%s  species=%d\n",
              x$from$method %||% "custom", x$from$frame_species))
  sr <- range(x$support)
  cat(sprintf("  windows: %s\n", paste(x$windows$label, collapse = ", ")))
  cat(sprintf("  support per cell: %d-%d plots\n", sr[1L], sr[2L]))
  invisible(x)
}

#' Tidy a species trajectory
#'
#' Long-to-wide data frame with one row per focal-species-by-window cell: the
#' window label and center, the support (plots backing the cell), and the `dim`
#' coordinate columns `d1..dD`.
#'
#' @param x A `specvec_trajectory`.
#' @param ... Unused.
#' @param na.rm Drop cells with no co-occurrence (all-`NA` coordinates).
#' @return A data frame.
#' @export
as.data.frame.specvec_trajectory <- function(x, ..., na.rm = FALSE) {
  sp <- dimnames(x$U)[[1L]]; wl <- dimnames(x$U)[[2L]]
  nF <- length(sp); nW <- length(wl); D <- dim(x$U)[3L]
  grid <- expand.grid(species = sp, window = wl, stringsAsFactors = FALSE,
                      KEEP.OUT.ATTRS = FALSE)
  si <- match(grid$species, sp); wi <- match(grid$window, wl)
  out <- data.frame(
    species = grid$species, window = grid$window,
    center  = x$windows$center[match(grid$window, x$windows$label)],
    support = x$support[cbind(si, wi)],
    stringsAsFactors = FALSE)
  for (d in seq_len(D)) {
    Md <- matrix(x$U[, , d], nF, nW)
    out[[paste0("d", d)]] <- Md[cbind(si, wi)]
  }
  if (isTRUE(na.rm)) {
    dcols <- grep("^d[0-9]+$", names(out))
    out <- out[stats::complete.cases(out[, dcols, drop = FALSE]), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

#' Community trajectory and novelty over time
#'
#' Embed each time window's communities in one fixed frame and score how novel
#' they are relative to a reference window. The species frame is fitted once on
#' all plots, so the per-window community embeddings live in the same space and
#' novelty is comparable across time (the temporal counterpart of
#' [community_novelty()]).
#'
#' @param x A `specvec_data` object built with a `time` column.
#' @param by Window definition (see [species_trajectory()]).
#' @param reference Reference window the novelty is measured against: `NULL`
#'   (default) uses the first window, a window label or index picks one, or a
#'   `specvec_community`/matrix supplies an external baseline.
#' @param method,dim,min_occurrence,min_cooccurrence Frame embedding controls.
#' @param weights `"cover"` (default when present) or `"presence"` pooling.
#' @param normalize If `TRUE`, L2-normalize each community vector before scoring.
#' @param k Neighbours averaged in the novelty distance.
#' @param frame_embedding Optional pre-fitted frame embedding to reuse.
#' @param ... Passed to [species_embedding()] when fitting the frame.
#' @return A `specvec_community_trajectory`: `communities` (per-window plot x dim
#'   matrices in the shared frame), a `novelty` table (window, center, n_plots,
#'   mean and median per-plot novelty vs the reference), the fixed `frame`, and
#'   the resolved `reference` label.
#' @seealso [species_trajectory()], [community_novelty()]
#' @export
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   plot = rep(paste0("p", 1:90), each = 3),
#'   species = sample(paste0("s", 1:15), 270, replace = TRUE),
#'   decade = rep(c(1990, 2000, 2010), each = 90)
#' )
#' x <- specvec(df, "plot", "species", time = "decade")
#' community_trajectory(x, dim = 4, min_occurrence = 1, k = 3)
community_trajectory <- function(x, by = NULL, reference = NULL,
                                 method = "abund_pmi", dim = 64L,
                                 weights = c("cover", "presence"), normalize = FALSE,
                                 k = 5L, min_occurrence = 5L, min_cooccurrence = 1L,
                                 frame_embedding = NULL, ...) {
  if (!inherits(x, "specvec_data")) stop("`x` must be a specvec_data object.", call. = FALSE)
  weights <- match.arg(weights)
  if (weights == "cover" && is.null(x$COV)) {
    message("specvec: no cover in data; community pooling uses presence.")
    weights <- "presence"
  }
  if (is.null(frame_embedding)) {
    frame_embedding <- species_embedding(x, method = method, dim = dim,
                                         min_occurrence = min_occurrence,
                                         min_cooccurrence = min_cooccurrence, ...)
  } else if (!inherits(frame_embedding, "specvec_embedding")) {
    stop("`frame_embedding` must be a specvec_embedding.", call. = FALSE)
  }

  Vf <- frame_embedding$V
  ks <- match(rownames(Vf), x$species)
  use_cover <- weights == "cover" && !is.null(x$COV)
  Mfull <- if (use_cover) x$COV else x$P

  win <- .time_windows(x, by)
  nW <- length(win$rows)
  communities <- vector("list", nW)
  names(communities) <- win$labels
  for (w in seq_len(nW)) {
    rows <- win$rows[[w]]
    U <- .pool_rows(Mfull[rows, ks, drop = FALSE], Vf)
    if (normalize) U <- .l2_normalize_rows(U)
    rownames(U) <- x$plots[rows]
    communities[[w]] <- U
  }

  ref_label <- .resolve_reference(reference, win, communities)
  Uref <- ref_label$matrix
  k <- min(as.integer(k), nrow(Uref))

  nov_mean <- nov_med <- numeric(nW)
  for (w in seq_len(nW)) {
    obj <- structure(list(U = communities[[w]]), class = "specvec_community")
    nv <- community_novelty(obj, Uref, k = k)
    nov_mean[w] <- mean(nv); nov_med[w] <- stats::median(nv)
  }
  novelty <- data.frame(
    window = win$labels, center = win$centers,
    n_plots = vapply(win$rows, length, integer(1)),
    mean_novelty = nov_mean, median_novelty = nov_med,
    stringsAsFactors = FALSE)

  structure(list(
    communities = communities, novelty = novelty, frame = frame_embedding,
    reference = ref_label$label, weights = weights,
    from = list(method = frame_embedding$method %||% method, dim = ncol(Vf))
  ), class = "specvec_community_trajectory")
}

## Resolve the novelty reference into a plot x dim matrix and a label. NULL ->
## first window; a window label/index -> that window; a community/matrix ->
## external baseline.
.resolve_reference <- function(reference, win, communities) {
  if (is.null(reference))
    return(list(matrix = communities[[1L]], label = win$labels[1L]))
  if (inherits(reference, "specvec_community") || is.matrix(reference))
    return(list(matrix = .community_matrix(reference, "reference"), label = "external"))
  if (is.numeric(reference) && length(reference) == 1L) {
    i <- as.integer(reference)
    if (i < 1L || i > length(communities)) stop("`reference` window index out of range.", call. = FALSE)
    return(list(matrix = communities[[i]], label = win$labels[i]))
  }
  i <- match(as.character(reference), win$labels)
  if (is.na(i)) stop(sprintf("`reference` window '%s' not found.", reference), call. = FALSE)
  list(matrix = communities[[i]], label = win$labels[i])
}

#' @export
print.specvec_community_trajectory <- function(x, ...) {
  cat(sprintf("<specvec_community_trajectory> windows=%d  dim=%d  weights=%s\n",
              nrow(x$novelty), x$from$dim, x$weights))
  cat(sprintf("  frame: method=%s  species=%d\n",
              x$from$method %||% "custom", nrow(x$frame$V)))
  cat(sprintf("  reference window: %s\n", x$reference))
  nv <- x$novelty
  cat(sprintf("  %-12s %8s %14s\n", "window", "n_plots", "mean_novelty"))
  for (i in seq_len(nrow(nv)))
    cat(sprintf("  %-12s %8d %14.4f\n", nv$window[i], nv$n_plots[i], nv$mean_novelty[i]))
  invisible(x)
}
