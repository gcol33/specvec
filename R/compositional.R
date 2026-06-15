## Compositional weighting (philosophy B, the research-track fork of section 10).
## Philosophy A (abund_pmi) weights co-occurrence by cover and builds weighted
## PMI. Philosophy B treats each plot as a composition and works in Aitchison
## geometry: the species x species operator is the centered-log-ratio (CLR)
## covariance, estimated from the variation matrix of pairwise log-ratios. This
## is a genuinely distinct operand (it can be negative, it is dense), which is
## why it registers through the same capability machinery rather than reusing the
## PMI path. The benchmark (compare_embeddings) decides whether A or B is the
## default; the comparison is itself the methods-paper result.

## Variation matrix T over the kept species: T[i,j] = Var_p( log(cov[p,i] /
## cov[p,j]) ) over plots where both i and j are present. The CLR covariance is
## recovered by double-centering, G = -1/2 J T J, since clr_i - clr_j =
## log(cov_i / cov_j) makes T[i,j] = Sigma_ii + Sigma_jj - 2 Sigma_ij (the
## classical distance-to-Gram identity). Pairs co-present in fewer than two plots
## have no defined log-ratio variance and are filled with the largest observed
## variation (compositionally most dissimilar). The pairwise co-present estimate
## of T (each pair on its own plot set) is the same approximation SparCC makes;
## it is what lets the operator be built from sparse data without zero-imputing
## the full plot x species table.
.clr_covariance <- function(COV, species, min_copresent = 2L) {
  S <- length(species)
  COVt <- methods::as(methods::as(COV, "generalMatrix"), "TsparseMatrix")
  ii <- COVt@i + 1L; jj <- COVt@j + 1L
  ## presence and log-cover on the present pattern. log(cover)=0 at cover=1 is a
  ## harmless zero: it contributes 0 to every sum it appears in, and presence is
  ## tracked separately through B.
  B <- Matrix::sparseMatrix(i = ii, j = jj, x = rep(1, length(ii)),
                            dims = dim(COV))
  L <- Matrix::sparseMatrix(i = ii, j = jj, x = log(COVt@x), dims = dim(COV))
  Nn <- as.matrix(Matrix::crossprod(B))            # co-present plot counts
  M1 <- as.matrix(Matrix::crossprod(L, B))         # sum over co-present of log cov_i
  P2 <- as.matrix(Matrix::crossprod(L))            # sum over co-present of log cov_i log cov_j
  Q  <- as.matrix(Matrix::crossprod(L * L, B))     # sum over co-present of (log cov_i)^2

  sum_diff <- M1 - t(M1)                            # sum over co-present of (log cov_i - log cov_j)
  sum_sq   <- Q + t(Q) - 2 * P2                     # sum over co-present of (log cov_i - log cov_j)^2
  Tvar <- sum_sq / Nn - (sum_diff / Nn)^2           # population variance over co-present plots
  defined <- Nn >= min_copresent
  Tvar[!defined] <- NA_real_
  diag(Tvar) <- 0
  fill <- if (any(is.finite(Tvar))) max(Tvar[is.finite(Tvar)]) else 0
  Tvar[!is.finite(Tvar)] <- fill
  Tvar <- (Tvar + t(Tvar)) / 2                      # symmetrize numerical drift

  ## Double-center -> CLR covariance: G = -1/2 J T J with J = I - 11'/S is, for
  ## symmetric T, the O(S^2) row/column-mean centering (no S^3 matrix product).
  rm <- rowMeans(Tvar); gm <- mean(rm)
  G <- -0.5 * (Tvar - outer(rm, rm, "+") + gm)
  G <- (G + t(G)) / 2
  dimnames(G) <- list(species, species)
  G
}

.w_clr <- function(data, ks, n_plots, min_cooccurrence) {
  if (is.null(data$COV)) {
    message("specvec: no cover in data; 'clr' has no composition and falls back to presence PMI.")
    return(.w_ppmi(data, ks, n_plots, min_cooccurrence))
  }
  species <- data$species[ks]
  min_co <- max(2L, as.integer(min_cooccurrence))   # a pair's log-ratio variance needs >= 2 co-present
  G <- .clr_covariance(data$COV[, ks, drop = FALSE], species, min_copresent = min_co)
  list(kind = "sym", M = G, species = species)
}
