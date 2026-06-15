`%||%` <- function(a, b) if (is.null(a)) b else a

## Fix the sign ambiguity of eigen/SVD columns: flip each column so its
## largest-magnitude entry is positive. Deterministic, no RNG, so two fits of
## the same operator give byte-identical matrices. Magnitude-invariant for
## dot products and cosine, so it never changes similarity results.
.sign_orient <- function(M) {
  if (!is.matrix(M) || nrow(M) == 0L || ncol(M) == 0L) return(M)
  for (j in seq_len(ncol(M))) {
    col <- M[, j]
    i <- which.max(abs(col))
    if (length(i) == 1L && is.finite(col[i]) && col[i] < 0) M[, j] <- -col
  }
  M
}

.zeros <- function(nr, dim, species) {
  m <- matrix(0, nr, dim)
  rownames(m) <- species
  m
}

.l2_normalize_rows <- function(M) {
  n <- sqrt(rowSums(M * M))
  n[!is.finite(n) | n == 0] <- 1
  M / n
}

## Cosine similarity between rows of A (and B, or A with itself).
.cosine_sim <- function(A, B = NULL) {
  An <- .l2_normalize_rows(A)
  if (is.null(B)) return(tcrossprod(An))
  tcrossprod(An, .l2_normalize_rows(B))
}

## Pool species vectors into one row per plot: weighted mean of the species
## vectors present in each plot. W is plot x species (cover or presence), V is
## species x dim. The single readout shared by community_embedding and the
## benchmark's presence-mean community definition.
.pool_rows <- function(W, V) {
  as.matrix(W %*% V) / pmax(as.numeric(Matrix::rowSums(W)), 1)
}

## Distance between two embedding points. Euclidean by default; "cosine" returns
## 1 - cosine similarity, so a smaller value means closer under either metric.
.vec_dist <- function(u, v, metric = c("euclidean", "cosine")) {
  metric <- match.arg(metric)
  if (metric == "cosine") {
    nu <- sqrt(sum(u * u)); nv <- sqrt(sum(v * v))
    if (nu == 0 || nv == 0) return(NA_real_)
    return(1 - sum(u * v) / (nu * nv))
  }
  sqrt(sum((u - v)^2))
}

## Pairwise Euclidean distance between rows of A and rows of B.
.pairwise_euclid <- function(A, B) {
  a2 <- rowSums(A * A); b2 <- rowSums(B * B)
  d2 <- outer(a2, b2, "+") - 2 * tcrossprod(A, B)
  d2[d2 < 0] <- 0
  sqrt(d2)
}

## Coerce a specvec_community or a plain matrix to its plot x dim matrix.
.community_matrix <- function(z, arg = "object") {
  if (inherits(z, "specvec_community")) return(z$U)
  if (is.matrix(z)) return(z)
  stop(sprintf("`%s` must be a specvec_community or a matrix.", arg), call. = FALSE)
}
