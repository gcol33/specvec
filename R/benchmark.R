## compare_embeddings(): the shipped bake-off. Scores registered methods on
## neutral ecological tasks under one fair protocol (section 6.7): a plot-level
## train/test split, one scored pair set per seed, never-seen negatives, and a
## single presence-mean community definition shared across methods. The species
## set is filtered once (globally) before splitting; each method is then fit on
## the training plots through the same engine path as species_embedding().

## Scale guards and kNN neighbourhood. Defaults reproduce bakeoff.R.
.benchmark_control <- function(pair_all_max = 2.5e6, pair_sample = 2.0e6,
                               auc_cap = 50000L, eunis_ref_cap = 30000L,
                               eunis_q_cap = 10000L, knn_k = 15L) {
  list(pair_all_max = pair_all_max, pair_sample = pair_sample,
       auc_cap = as.integer(auc_cap), eunis_ref_cap = as.integer(eunis_ref_cap),
       eunis_q_cap = as.integer(eunis_q_cap), knn_k = as.integer(knn_k))
}

## A specvec_data restricted to a set of plot rows, keeping every species column
## so the global kept-species index stays valid. Fed to .fit_embedding().
.subset_plots <- function(x, rows) {
  list(P = x$P[rows, , drop = FALSE],
       COV = if (is.null(x$COV)) NULL else x$COV[rows, , drop = FALSE],
       species = x$species, plots = x$plots[rows])
}

## Plot-level label vector for EUNIS, aligned to plot row order. Accepts a column
## name in x$labels, a full-length vector, or the first label column by default.
.resolve_plot_labels <- function(x, labels) {
  col <- NULL
  if (is.null(labels)) {
    if (!is.null(x$labels) && ncol(x$labels) >= 1L) col <- x$labels[[1L]]
  } else if (length(labels) == 1L && is.character(labels)) {
    if (!is.null(x$labels) && labels %in% names(x$labels)) col <- x$labels[[labels]]
  } else if (length(labels) == nrow(x$P)) {
    col <- labels
  }
  if (is.null(col)) return(NULL)
  v <- as.character(col)
  v[is.na(v) | v %in% c("", "~", "NA")] <- NA
  if (all(is.na(v))) return(NULL)
  v
}

## Species x trait numeric matrix aligned to the kept species, or NULL.
.resolve_traits <- function(traits, species) {
  if (is.null(traits)) return(NULL)
  if (is.data.frame(traits)) traits <- as.matrix(traits)
  if (!is.matrix(traits) || is.null(rownames(traits))) return(NULL)
  common <- intersect(species, rownames(traits))
  if (length(common) == 0L) return(NULL)
  out <- matrix(NA_real_, length(species), ncol(traits),
                dimnames = list(species, colnames(traits)))
  tt <- traits[common, , drop = FALSE]
  out[common, ] <- matrix(suppressWarnings(as.numeric(tt)), nrow(tt), ncol(tt))
  out
}

.metric_map <- list(
  cooc_ppmi = "cooc_ppmi", cooc_raw = "cooc_raw", link_auc = "link_auc",
  eunis = c("eunis_f1", "eunis_acc", "eunis_base"), trait = "trait_r2"
)

#' Benchmark embedding methods
#'
#' Score registered embedding methods on neutral ecological tasks under one fair
#' protocol: a plot-level train/test split, one scored species-pair set per seed,
#' never-seen negatives for link prediction, and a single presence-mean community
#' definition shared across methods. Co-occurrence metrics need only `x`; the
#' EUNIS and trait metrics activate when labels or a trait table are supplied.
#' The species set is filtered once (`min_occurrence`) before splitting, so every
#' method sees the same species. Runs on any user dataset, not just the reference
#' data that tuned the package default.
#'
#' @param x A `specvec_data` object.
#' @param methods Method names to compare (see [specvec_methods()]). `"glove"`
#'   is skipped with a note if `text2vec` is not installed.
#' @param metrics Which metrics to report: any of `"cooc_ppmi"`, `"cooc_raw"`,
#'   `"link_auc"`, `"eunis"`, `"trait"`. EUNIS/trait skip cleanly when their
#'   inputs are absent.
#' @param dim Embedding dimension.
#' @param seeds Integer seeds; each gives an independent split. Reported as
#'   across-seed mean and SD.
#' @param test_frac Held-out plot fraction (and species fraction for traits).
#' @param min_occurrence,min_cooccurrence Species/pair filters, applied globally
#'   before splitting.
#' @param labels Plot-level label for the EUNIS metric: a column name in
#'   `x$labels`, a plot-ordered vector, or `NULL` to use the first label column.
#' @param traits Species-level numeric attributes for the trait metric: a matrix
#'   or data frame with species row names, or `NULL` to skip.
#' @param reference Method the verdict compares against (default `"ca"`).
#' @param glove_iter Iterations for the GloVe factorizer.
#' @param control Scale guards and kNN neighbourhood; see the defaults in the
#'   source. Override only for very large data.
#' @return A `specvec_benchmark`: a tidy `method x metric` summary (mean, sd, n),
#'   the per-seed raw rows, the run config, and the fixed-criterion verdict. A
#'   method "beats" the reference only if it exceeds it on both neutral metrics
#'   (`cooc_raw` and `link_auc`) by more than two pooled SDs.
#' @export
#' @examples
#' df <- data.frame(
#'   plot = rep(paste0("p", 1:40), each = 3),
#'   species = sample(paste0("s", 1:15), 120, replace = TRUE),
#'   cover = round(runif(120, 1, 100), 1)
#' )
#' x <- specvec(df, "plot", "species", abundance = "cover")
#' \donttest{
#' compare_embeddings(x, methods = c("ca", "pmi", "abund_pmi"),
#'                    dim = 4, seeds = 1:2, min_occurrence = 1)
#' }
compare_embeddings <- function(x,
                               methods = c("ca", "pmi", "abund_pmi", "glove"),
                               metrics = c("cooc_ppmi", "cooc_raw", "link_auc",
                                           "eunis", "trait"),
                               dim = 64L, seeds = 1:3, test_frac = 0.2,
                               min_occurrence = 5L, min_cooccurrence = 1L,
                               labels = NULL, traits = NULL,
                               reference = "ca", glove_iter = 20L,
                               control = .benchmark_control()) {
  if (!inherits(x, "specvec_data")) stop("`x` must be a specvec_data object.", call. = FALSE)
  metrics <- match.arg(metrics, several.ok = TRUE)
  dim <- as.integer(dim); seeds <- as.integer(seeds)
  min_occurrence <- as.integer(min_occurrence); min_cooccurrence <- as.integer(min_cooccurrence)
  methods <- unique(methods)

  specs <- list()
  for (m in methods) {
    mm <- .get_method(m)
    if (mm$factorization == "glove" && !requireNamespace("text2vec", quietly = TRUE)) {
      message(sprintf("specvec: method '%s' needs 'text2vec' (Suggests); skipping.", m))
      next
    }
    specs[[m]] <- mm
  }
  if (length(specs) == 0L) stop("no runnable methods.", call. = FALSE)

  ks <- .kept_species(x, min_occurrence)
  S <- length(ks)
  species <- x$species[ks]
  BIN <- x$P[, ks, drop = FALSE]
  M <- nrow(BIN)
  if (S <= dim + 1L)
    stop(sprintf("species kept (%d) must exceed dim+1 (%d); lower dim or min_occurrence.",
                 S, dim + 1L), call. = FALSE)

  want_eunis <- "eunis" %in% metrics
  want_trait <- "trait" %in% metrics
  eun_map <- if (want_eunis) .resolve_plot_labels(x, labels) else NULL
  if (want_eunis && is.null(eun_map)) {
    message("specvec: no usable plot labels for 'eunis'; skipping that metric.")
    want_eunis <- FALSE
  }
  TRAITM <- if (want_trait) .resolve_traits(traits, species) else NULL
  if (want_trait && is.null(TRAITM)) {
    message("specvec: no species traits for 'trait'; skipping that metric.")
    want_trait <- FALSE
  }

  rows <- list()
  for (seed in seeds) {
    set.seed(seed)
    te_plot <- sort(sample.int(M, floor(M * test_frac)))
    tr_plot <- setdiff(seq_len(M), te_plot)
    sp_te <- sort(sample.int(S, floor(S * test_frac)))
    sp_tr <- setdiff(seq_len(S), sp_te)
    if (length(tr_plot) < 2L || length(te_plot) < 1L)
      stop("test_frac leaves too few train/test plots.", call. = FALSE)

    ev <- .make_eval(BIN, tr_plot, te_plot, control)
    data_tr <- .subset_plots(x, tr_plot)

    for (m in names(specs)) {
      mm <- specs[[m]]
      emb <- tryCatch(
        .fit_embedding(data_tr, ks, mm$weighting, mm$factorization, dim,
                       n_plots = length(tr_plot),
                       min_cooccurrence = min_cooccurrence, glove_iter = glove_iter),
        error = function(e) {
          message(sprintf("specvec: method '%s' (seed %d) failed: %s",
                          m, seed, conditionMessage(e)))
          NULL
        })
      if (is.null(emb)) next
      co  <- .co_metrics(emb, ev)
      eu  <- if (want_eunis) .eunis_metric(emb, BIN, eun_map, tr_plot, te_plot, control)
             else c(eunis_f1 = NA_real_, eunis_acc = NA_real_, eunis_base = NA_real_)
      trr <- if (want_trait) .trait_metric(emb, TRAITM, sp_tr, sp_te, control) else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        seed = seed, method = m,
        cooc_ppmi = unname(co["cooc_ppmi"]), cooc_raw = unname(co["cooc_raw"]),
        link_auc = unname(co["link_auc"]),
        eunis_f1 = unname(eu["eunis_f1"]), eunis_acc = unname(eu["eunis_acc"]),
        eunis_base = unname(eu["eunis_base"]), trait_r2 = trr,
        stringsAsFactors = FALSE)
    }
  }
  if (length(rows) == 0L) stop("no method produced results.", call. = FALSE)
  raw <- do.call(rbind, rows)

  display_cols <- unique(unlist(.metric_map[metrics]))
  display_cols <- display_cols[vapply(display_cols, function(cc) any(is.finite(raw[[cc]])),
                                      logical(1))]
  agg <- list()
  for (m in unique(raw$method)) for (cc in display_cols) {
    v <- raw[raw$method == m, cc]; v <- v[is.finite(v)]
    agg[[length(agg) + 1L]] <- data.frame(
      method = m, metric = cc, mean = if (length(v)) mean(v) else NA_real_,
      sd = if (length(v) > 1L) stats::sd(v) else 0, n = length(v),
      stringsAsFactors = FALSE)
  }
  summary_tbl <- do.call(rbind, agg)
  verdict <- .benchmark_verdict(summary_tbl, reference, neutral = c("cooc_raw", "link_auc"))

  structure(list(
    summary = summary_tbl, raw = raw, verdict = verdict,
    config = list(methods = names(specs), metrics = metrics, dim = dim, seeds = seeds,
                  test_frac = test_frac, min_occurrence = min_occurrence,
                  min_cooccurrence = min_cooccurrence, reference = reference,
                  n_plots = M, n_species = S, control = control)
  ), class = "specvec_benchmark")
}

#' @rdname compare_embeddings
#' @export
specvec_benchmark <- compare_embeddings

## Fixed win criterion: a method beats the reference iff it exceeds it on every
## neutral metric by more than two pooled SDs. Reported regardless of the result.
.benchmark_verdict <- function(summary_tbl, reference, neutral) {
  methods <- unique(summary_tbl$method)
  get1 <- function(m, cc, fld) {
    r <- summary_tbl[summary_tbl$method == m & summary_tbl$metric == cc, ]
    if (nrow(r) == 0L) NA_real_ else r[[fld]]
  }
  beats <- function(m1, s1, m0, s0) {
    if (any(!is.finite(c(m1, s1, m0, s0)))) return(NA)
    isTRUE((m1 - m0) > 2 * sqrt((s1^2 + s0^2) / 2))
  }
  has_ref <- reference %in% methods &&
    all(neutral %in% summary_tbl$metric[summary_tbl$method == reference])
  out <- list(reference = reference, has_reference = has_ref, neutral = neutral, rows = list())
  if (!has_ref) return(out)
  for (m in setdiff(methods, reference)) {
    per <- vapply(neutral, function(cc)
      beats(get1(m, cc, "mean"), get1(m, cc, "sd"),
            get1(reference, cc, "mean"), get1(reference, cc, "sd")), logical(1))
    delta <- vapply(neutral, function(cc)
      get1(m, cc, "mean") - get1(reference, cc, "mean"), numeric(1))
    out$rows[[m]] <- list(method = m, beats = per, delta = delta,
                          all = isTRUE(all(per)), any = isTRUE(any(per)))
  }
  out
}

#' @export
print.specvec_benchmark <- function(x, ...) {
  cfg <- x$config
  cat(sprintf("<specvec_benchmark> plots=%d  species=%d  dim=%d  seeds=%s\n",
              cfg$n_plots, cfg$n_species, cfg$dim, paste(cfg$seeds, collapse = ",")))
  cat(sprintf("  methods: %s\n", paste(cfg$methods, collapse = ", ")))

  s <- x$summary
  metrics <- unique(s$metric)
  ord_metric <- if ("link_auc" %in% metrics) "link_auc" else metrics[1]
  ms <- unique(s$method)
  ord <- vapply(ms, function(m) {
    r <- s[s$method == m & s$metric == ord_metric, ]
    if (nrow(r)) r$mean else -Inf
  }, numeric(1))
  ms <- ms[order(-ord)]

  cat("\n", sprintf("  %-10s", "method"), sep = "")
  for (cc in metrics) cat(sprintf(" %16s", cc))
  cat("\n")
  for (m in ms) {
    cat(sprintf("  %-10s", m))
    for (cc in metrics) {
      r <- s[s$method == m & s$metric == cc, ]
      cell <- if (nrow(r) && is.finite(r$mean)) sprintf("%.3f+-%.3f", r$mean, r$sd) else "-"
      cat(sprintf(" %16s", cell))
    }
    cat("\n")
  }

  v <- x$verdict
  cat(sprintf("\n  verdict vs '%s' (beats on %s by >2 pooled SDs):\n",
              v$reference, paste(v$neutral, collapse = " & ")))
  if (!isTRUE(v$has_reference)) {
    cat(sprintf("    reference '%s' absent; ordering only.\n", v$reference))
  } else if (length(v$rows) == 0L) {
    cat("    no comparison methods.\n")
  } else {
    for (m in names(v$rows)) {
      rr <- v$rows[[m]]
      tag <- if (rr$all) "BEATS" else if (rr$any) "partial" else "no advantage"
      d <- paste(sprintf("%s%+.3f[%s]", v$neutral, rr$delta,
                         ifelse(is.na(rr$beats), "?", ifelse(rr$beats, "Y", "n"))),
                 collapse = "  ")
      cat(sprintf("    %-10s %s  -> %s\n", m, d, tag))
    }
  }
  invisible(x)
}
