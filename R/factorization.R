## Factorizations consume the operator a weighting produced and return a
## species x dim matrix with species row names.

## Symmetric eigendecomposition of a PPMI-style operator; embedding = top
## eigenvectors scaled by sqrt of the positive eigenvalues, zero-padded to dim.
.f_eigen <- function(op, dim, ...) {
  PP <- op$M
  S <- nrow(PP)
  nnz <- if (methods::is(PP, "sparseMatrix")) length(PP@x) else sum(PP != 0)
  if (nnz == 0L) return(.zeros(S, dim, op$species))
  k <- min(dim, S - 1L)
  e <- NULL
  if (k >= 1L && S > 2L) {
    e <- tryCatch(RSpectra::eigs_sym(PP, k = k, which = "LA"),
                  error = function(err) NULL)
  }
  if (is.null(e) || all(!is.finite(e$values)) || max(e$values, na.rm = TRUE) <= 0) {
    PPd <- as.matrix(PP); PPd <- (PPd + t(PPd)) / 2
    e <- eigen(PPd, symmetric = TRUE)
  }
  pos <- which(e$values > 1e-8)
  k2 <- min(dim, length(pos))
  if (k2 < 1L) return(.zeros(S, dim, op$species))
  idx <- pos[seq_len(k2)]
  emb <- e$vectors[, idx, drop = FALSE] %*% diag(sqrt(e$values[idx]), k2, k2)
  if (k2 < dim) emb <- cbind(emb, matrix(0, S, dim - k2))
  rownames(emb) <- op$species
  emb
}

## Implicit-matvec SVD (CA): no densify of the plot x species contingency.
.f_svd <- function(op, dim, ...) {
  k <- min(dim, op$ncol - 1L, op$nrow - 1L)
  if (k < 1L) return(.zeros(op$ncol, dim, op$species))
  sv <- tryCatch(
    RSpectra::svds(op$Af, k = k, nu = 0, nv = k, dim = c(op$nrow, op$ncol),
                   Atrans = op$Atf, args = NULL),
    error = function(err) NULL)
  if (is.null(sv)) return(.zeros(op$ncol, dim, op$species))
  emb <- op$readout(sv$v, sv$d)
  if (ncol(emb) < dim) emb <- cbind(emb, matrix(0, nrow(emb), dim - ncol(emb)))
  rownames(emb) <- op$species
  emb
}

## GloVe objective-factorizer: own weighted least-squares loss, bias terms,
## x_max. Brings its own weighting, so it consumes the raw count operator.
.f_glove <- function(op, dim, glove_iter = 20L, x_max = 10, ...) {
  if (!requireNamespace("text2vec", quietly = TRUE)) {
    stop("method 'glove' requires the 'text2vec' package; install it or choose ",
         "another method (e.g. 'abund_pmi').", call. = FALSE)
  }
  ## GloVe logs per-epoch loss at INFO through the rsparse logger; it propagates
  ## to lgr's root console appender, which neither suppressMessages nor a logger
  ## threshold silences. Raise that appender's threshold for the fit and restore.
  if (requireNamespace("lgr", quietly = TRUE)) {
    ap <- lgr::get_logger("root")$appenders$console
    if (!is.null(ap)) {
      old_ap <- ap$threshold
      ap$set_threshold("warn")
      on.exit(try(ap$set_threshold(old_ap), silent = TRUE), add = TRUE)
    }
  }
  tcm <- methods::as(op$M, "TsparseMatrix")
  Matrix::diag(tcm) <- 0
  dimnames(tcm) <- list(op$species, op$species)
  gv <- text2vec::GlobalVectors$new(rank = dim, x_max = x_max)
  wv <- suppressMessages(gv$fit_transform(tcm, n_iter = glove_iter,
                                          convergence_tol = 0.01, progressbar = FALSE))
  emb <- (wv + t(gv$components))[op$species, , drop = FALSE]
  emb <- as.matrix(emb)
  rownames(emb) <- op$species
  emb
}
