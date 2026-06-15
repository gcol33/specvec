## Alien integration trajectory (v0.3, the flagship application). A thin
## specialization of the v0.2 fixed-frame primitive: the focal neophyte is placed
## per window by species_trajectory(), the native-community centroid is pooled per
## window with .pool_rows() in the same frame, and the readout is the distance
## between them over time. No placement or windowing logic is re-implemented here.

#' Alien integration trajectory
#'
#' Track how a neophyte moves toward (or away from) the native community over
#' time, in one fixed embedding frame. The focal species is placed per time
#' window at the cover-weighted centroid of the species it co-occurs with (the
#' [species_trajectory()] projection), the native community is placed per window
#' at the centroid of its pooled plot embeddings ([community_embedding()]'s
#' readout), and the trajectory is the focal-to-native distance per window. A
#' falling distance is integration: the neophyte's associates shift from where it
#' arrived toward the resident native community.
#'
#' @details
#' Two species sets shape the measurement.
#'
#' * `frame` is the fixed coordinate system, fitted once. It defaults to every
#'   species except the focal ones, so the neophyte moves through the resident
#'   pool rather than helping define it. Keeping the whole resident pool in the
#'   frame (not the native subset alone) means the pole the neophyte starts from
#'   -- the disturbed or ruderal flora it arrives with -- is represented, so an
#'   early position is well defined and the motion toward the natives is real.
#' * `native` is the subset of the frame whose centroid is the integration
#'   target. It defaults to the whole frame (the resident community at large); for
#'   a meaningful measurement supply the species that are native at the study
#'   region (e.g. the rows flagged native in an EVA `STATUS` column).
#'
#' **ReSurvey anchoring.** Raw multi-decade plot data confounds integration with
#' where and what was sampled each decade. Restrict `x` to resampled plots (the
#' `ReSurvey plot (Y/N)` flag in EVA header exports) before building the
#' `specvec_data`, so the same locations are tracked through time. The function
#' takes whatever plots it is given and does not hard-code a survey-design column;
#' the restriction is a data-preparation step, demonstrated in the science
#' vignette.
#'
#' @param x A `specvec_data` object built with a `time` column.
#' @param species Character vector of focal neophyte ids to trace.
#' @param native Character vector of native species defining the integration
#'   target; default is every frame species (the resident community at large).
#' @param frame Optional character vector fixing the coordinate system; default
#'   is every species except `species`.
#' @param by Window definition passed to the time splitter: `NULL` (default) for
#'   one window per distinct time value, or a numeric break vector to bin time.
#' @param method Frame embedding method (see [specvec_methods()]).
#' @param dim Frame embedding dimension.
#' @param weights `"cover"` (default when cover is present) weights co-occurrence
#'   and pooling by cover; `"presence"` counts plots.
#' @param metric Distance from the focal species to the native centroid:
#'   `"euclidean"` (default) or `"cosine"` (returned as `1 - cosine similarity`,
#'   so a smaller value still means closer for both).
#' @param min_occurrence,min_cooccurrence Frame species/pair filters.
#' @param frame_embedding Optional pre-fitted `specvec_embedding` to reuse as the
#'   frame instead of fitting one.
#' @param ... Passed to [species_embedding()] when fitting the frame.
#' @return A `specvec_integration`: `distance` (focal x window distance to the
#'   native centroid, `NA` where the focal has no co-occurrence or a window holds
#'   no native community), `support` (focal x window plot counts), the per-window
#'   `native_centroid` and `native_support`, the `windows` table, the underlying
#'   `specvec_trajectory`, and the fixed `frame`. `as.data.frame()` tidies it.
#' @seealso [species_trajectory()], [community_trajectory()]
#' @export
#' @examples
#' set.seed(1)
#' native  <- paste0("nat", 1:6)
#' ruderal <- paste0("rud", 1:6)
#' rows <- list(); pid <- 0L
#' for (dec in c(1990, 2000, 2010)) {
#'   share_nat <- (dec - 1990) / 20             # the neophyte shifts ruderal -> native
#'   for (i in 1:50) {
#'     pid <- pid + 1L
#'     pool <- if (stats::runif(1) < share_nat) native else ruderal
#'     rows[[pid]] <- data.frame(plot = paste0("p", pid),
#'                               species = c(sample(pool, 3), "alien"), decade = dec)
#'   }
#'   for (i in 1:30) {                          # resident native backbone, every decade
#'     pid <- pid + 1L
#'     rows[[pid]] <- data.frame(plot = paste0("p", pid),
#'                               species = sample(native, 3), decade = dec)
#'   }
#' }
#' df <- do.call(rbind, rows)
#' x <- specvec(df, "plot", "species", time = "decade")
#' it <- integration_trajectory(x, species = "alien", native = native,
#'                              dim = 4, weights = "presence", min_occurrence = 1)
#' as.data.frame(it)
integration_trajectory <- function(x, species, native = NULL, frame = NULL,
                                   by = NULL, method = "abund_pmi", dim = 64L,
                                   weights = c("cover", "presence"),
                                   metric = c("euclidean", "cosine"),
                                   min_occurrence = 5L, min_cooccurrence = 1L,
                                   frame_embedding = NULL, ...) {
  if (!inherits(x, "specvec_data")) stop("`x` must be a specvec_data object.", call. = FALSE)
  weights <- match.arg(weights)
  metric  <- match.arg(metric)
  species <- as.character(species)

  ## Focal placement + the fixed frame come straight from species_trajectory: the
  ## neophyte is placed per window at the cover-weighted centroid of the frame
  ## species it co-occurs with. The placement maths is not re-implemented here.
  tr <- species_trajectory(x, species = species, frame = frame, by = by,
                           method = method, dim = dim, weights = weights,
                           min_occurrence = min_occurrence,
                           min_cooccurrence = min_cooccurrence,
                           frame_embedding = frame_embedding, ...)

  Vf <- tr$frame$V
  frame_species <- rownames(Vf)
  native <- if (is.null(native)) frame_species else as.character(native)
  native_in_frame <- intersect(native, frame_species)
  if (!length(native_in_frame))
    stop("no `native` species are in the fitted frame; check `native` and `min_occurrence`.",
         call. = FALSE)

  ## Native-community centroid per window: pool each window plot's native species
  ## in the same fixed frame (community_embedding's readout via .pool_rows), then
  ## average the plots that hold at least one native species. Same windows and
  ## same frame as the focal, so the focal-to-native distance lives in one space.
  use_cover <- tr$weights == "cover" && !is.null(x$COV)
  Mfull <- if (use_cover) x$COV else x$P
  nat_cols <- match(native_in_frame, x$species)
  Vnat <- Vf[native_in_frame, , drop = FALSE]

  win <- .time_windows(x, by)
  nW <- length(win$rows); D <- ncol(Vf)
  Cnat <- matrix(NA_real_, nW, D, dimnames = list(win$labels, NULL))
  nat_support <- stats::setNames(integer(nW), win$labels)
  for (w in seq_len(nW)) {
    Wn <- Mfull[win$rows[[w]], nat_cols, drop = FALSE]   # window plots x native
    has <- as.numeric(Matrix::rowSums(Wn)) > 0
    nat_support[w] <- sum(has)
    if (any(has)) Cnat[w, ] <- colMeans(.pool_rows(Wn[has, , drop = FALSE], Vnat))
  }

  ## focal x window distance to the native centroid.
  sp <- dimnames(tr$U)[[1L]]; nF <- length(sp)
  Dist <- matrix(NA_real_, nF, nW, dimnames = list(sp, win$labels))
  for (w in seq_len(nW)) {
    cen <- Cnat[w, ]
    if (anyNA(cen)) next
    for (f in seq_len(nF)) {
      u <- tr$U[f, w, ]
      if (!anyNA(u)) Dist[f, w] <- .vec_dist(u, cen, metric)
    }
  }

  structure(list(
    species = species, native = native_in_frame,
    distance = Dist, support = tr$support, metric = metric,
    native_centroid = Cnat, native_support = nat_support,
    windows = tr$windows, trajectory = tr, frame = tr$frame, weights = tr$weights,
    from = list(method = tr$from$method %||% method, dim = D,
                frame_species = length(frame_species),
                native_species = length(native_in_frame))
  ), class = "specvec_integration")
}

#' @export
print.specvec_integration <- function(x, ...) {
  cat(sprintf("<specvec_integration> focal=%d  windows=%d  dim=%d  metric=%s  weights=%s\n",
              length(x$species), nrow(x$windows), x$from$dim, x$metric, x$weights))
  cat(sprintf("  frame: method=%s  species=%d  native=%d\n",
              x$from$method %||% "custom", x$from$frame_species, x$from$native_species))
  for (f in x$species) {
    cat(sprintf("  %s -> native centroid:\n", f))
    cat(sprintf("    %-12s %8s %8s %12s\n", "window", "n_plots", "support", x$metric))
    for (i in seq_len(nrow(x$windows))) {
      d <- x$distance[f, i]
      cat(sprintf("    %-12s %8d %8d %12s\n",
                  x$windows$label[i], x$windows$n_plots[i], x$support[f, i],
                  if (is.na(d)) "NA" else formatC(d, format = "f", digits = 4)))
    }
  }
  invisible(x)
}

#' Tidy an integration trajectory
#'
#' Long data frame with one row per focal-species-by-window cell: the window
#' label and center, the plots in the window, the focal support, the count of
#' native-bearing plots, and the focal-to-native `distance`.
#'
#' @param x A `specvec_integration`.
#' @param ... Unused.
#' @param na.rm Drop cells with no measured distance.
#' @return A data frame.
#' @export
as.data.frame.specvec_integration <- function(x, ..., na.rm = FALSE) {
  sp <- x$species; wl <- x$windows$label
  grid <- expand.grid(species = sp, window = wl, stringsAsFactors = FALSE,
                      KEEP.OUT.ATTRS = FALSE)
  si <- match(grid$species, sp); wi <- match(grid$window, wl)
  out <- data.frame(
    species = grid$species, window = grid$window,
    center  = x$windows$center[wi],
    n_plots = x$windows$n_plots[wi],
    support = x$support[cbind(si, wi)],
    native_support = x$native_support[wi],
    distance = x$distance[cbind(si, wi)],
    stringsAsFactors = FALSE)
  if (isTRUE(na.rm)) out <- out[!is.na(out$distance), , drop = FALSE]
  out <- out[order(match(out$species, sp), out$center), , drop = FALSE]
  rownames(out) <- NULL
  out
}
