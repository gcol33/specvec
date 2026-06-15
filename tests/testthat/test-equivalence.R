## Dense CA reference: same chi-square-residual SVD, computed by densifying.
ca_dense <- function(N, dim) {
  P <- as.matrix(N); tot <- sum(P); P <- P / tot
  r <- rowSums(P); c <- colSums(P); cs <- ifelse(c > 0, c, NA)
  E <- outer(r, c); Sr <- (P - E) / sqrt(E); Sr[!is.finite(Sr)] <- 0
  sv <- svd(Sr, nu = 0, nv = dim)
  d <- sv$d[seq_len(dim)]; v <- sv$v[, seq_len(dim), drop = FALSE]
  emb <- sweep(v, 1, sqrt(cs), "/"); emb <- sweep(emb, 2, d, "*")
  emb[!is.finite(emb)] <- 0
  emb
}

test_that("implicit-matvec CA equals dense CA (gram, sign/rotation invariant)", {
  df <- sim_communities(M = 400L, S = 36L, seed = 3L)
  x <- specvec(df, "plot", "species")
  dim <- 10L

  emb_i <- species_embedding(x, method = "ca", dim = dim, min_occurrence = 1L)$V
  ks <- which(as.numeric(Matrix::colSums(x$P)) >= 1L)
  Nfull <- x$P[, ks, drop = FALSE]
  emb_d <- ca_dense(Nfull, dim)

  ## Gram matrix V V^T is invariant to per-column sign flips and to rotation
  ## within equal-singular-value subspaces, so it is the fair equivalence check.
  G_i <- emb_i %*% t(emb_i)
  G_d <- emb_d %*% t(emb_d)
  expect_gt(cor(as.numeric(G_i), as.numeric(G_d)), 0.999)
  expect_lt(max(abs(G_i - G_d)), 1e-3)
})

test_that("the same fit is deterministic (byte-identical V)", {
  df <- sim_communities(M = 300L, S = 30L, seed = 5L)
  x <- specvec(df, "plot", "species")
  a <- species_embedding(x, method = "abund_pmi", dim = 8L, min_occurrence = 1L)$V
  b <- species_embedding(x, method = "abund_pmi", dim = 8L, min_occurrence = 1L)$V
  expect_identical(a, b)
})
