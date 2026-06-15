## Weightings. Each builds the operator its paired factorization consumes:
##   kind "sym"      : $M sparse symmetric species x species         (eigen)
##   kind "counts"   : $M sparse species x species co-occurrence     (glove)
##   kind "implicit" : $Af/$Atf matvec + $readout, no densify        (svd / CA)

.cooc_counts <- function(P) Matrix::crossprod(P)   # species x species; diag = occurrence

## Cleaned sparse PPMI of a species x species association `C`.
## f = diag(C) is the marginal (occurrence count, or cover mass for AbundPMI);
## PPMI[a,b] = max(log(C[a,b] * n_plots / (f[a] f[b])), 0), a != b.
## `count_for_threshold` (binary co-occurrence counts) gates min_cooccurrence
## by number of plots, independent of the cover-weighted magnitude.
.ppmi_sparse <- function(C, n_plots, min_cooccurrence = 1L, count_for_threshold = NULL) {
  ## crossprod() yields a symmetric matrix that stores a single triangle; force
  ## to general so both (a,b) and (b,a) entries are present and the result is a
  ## genuinely symmetric operator for the eigensolver.
  C <- methods::as(methods::as(C, "generalMatrix"), "TsparseMatrix")
  i <- C@i + 1L; j <- C@j + 1L; x <- C@x
  f <- Matrix::diag(methods::as(C, "CsparseMatrix"))
  off <- i != j & x > 0 & f[i] > 0 & f[j] > 0
  if (min_cooccurrence > 1L) {
    if (!is.null(count_for_threshold)) {
      ct <- as.numeric(count_for_threshold[cbind(i, j)])
      off <- off & is.finite(ct) & ct >= min_cooccurrence
    } else {
      off <- off & x >= min_cooccurrence
    }
  }
  i <- i[off]; j <- j[off]; xx <- x[off]
  pp <- pmax(log(xx * n_plots / (f[i] * f[j])), 0)
  keep <- pp > 0
  Matrix::sparseMatrix(i = i[keep], j = j[keep], x = pp[keep], dims = dim(C),
                       dimnames = dimnames(C))
}

.w_counts <- function(data, ks, n_plots, min_cooccurrence) {
  P <- data$P[, ks, drop = FALSE]
  C <- methods::as(.cooc_counts(P), "generalMatrix")   # both triangles stored
  if (min_cooccurrence > 1L) {
    C <- methods::as(C, "TsparseMatrix")
    drop_edge <- C@i != C@j & C@x < min_cooccurrence
    C@x[drop_edge] <- 0
    C <- Matrix::drop0(C)
  }
  list(kind = "counts", M = methods::as(C, "CsparseMatrix"), species = data$species[ks])
}

.w_ppmi <- function(data, ks, n_plots, min_cooccurrence) {
  P <- data$P[, ks, drop = FALSE]
  C <- .cooc_counts(P)
  PP <- .ppmi_sparse(C, n_plots, min_cooccurrence,
                     count_for_threshold = if (min_cooccurrence > 1L) C else NULL)
  list(kind = "sym", M = PP, species = data$species[ks])
}

.w_abundance_pmi <- function(data, ks, n_plots, min_cooccurrence) {
  if (is.null(data$COV)) {
    message("specvec: no abundance/cover in data; 'abund_pmi' falls back to presence PMI.")
    return(.w_ppmi(data, ks, n_plots, min_cooccurrence))
  }
  W <- sqrt(data$COV[, ks, drop = FALSE])
  A <- Matrix::crossprod(W)             # A[a,b] = sum_p sqrt(cov[p,a] * cov[p,b])
  ct <- if (min_cooccurrence > 1L) .cooc_counts(data$P[, ks, drop = FALSE]) else NULL
  PP <- .ppmi_sparse(A, n_plots, min_cooccurrence, count_for_threshold = ct)
  list(kind = "sym", M = PP, species = data$species[ks])
}

## CA: chi-square-residual SVD applied implicitly to the plot x species
## contingency (no densify). Readout maps the species singular vectors to
## species coordinates. min_cooccurrence does not apply (operand is plot x species).
.w_chi_square <- function(data, ks, n_plots, min_cooccurrence) {
  N <- data$P[, ks, drop = FALSE]
  tot <- sum(N)
  r <- as.numeric(Matrix::rowSums(N)) / tot
  c <- as.numeric(Matrix::colSums(N)) / tot
  rs <- ifelse(r > 0, r, NA_real_)
  cs <- ifelse(c > 0, c, NA_real_)
  Af <- function(x, args) {
    w <- x / sqrt(cs); w[!is.finite(w)] <- 0
    pw <- as.numeric(N %*% w) / tot
    o <- (pw - r * sum(c * w)) / sqrt(rs); o[!is.finite(o)] <- 0; o
  }
  Atf <- function(x, args) {
    z <- x / sqrt(rs); z[!is.finite(z)] <- 0
    ptz <- as.numeric(Matrix::crossprod(N, z)) / tot
    o <- (ptz - c * sum(r * z)) / sqrt(cs); o[!is.finite(o)] <- 0; o
  }
  readout <- function(v, d) {
    emb <- sweep(v, 1, sqrt(cs), "/")
    emb <- sweep(emb, 2, d, "*")
    emb[!is.finite(emb)] <- 0
    emb
  }
  list(kind = "implicit", Af = Af, Atf = Atf, nrow = nrow(N), ncol = ncol(N),
       readout = readout, species = data$species[ks])
}
