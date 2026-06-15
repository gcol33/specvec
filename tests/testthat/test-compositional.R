## Compositional method (philosophy B, the section-10 research fork). The clr
## operator is the CLR covariance estimated from the variation matrix of pairwise
## log-ratios, double-centered. Two bars: it recovers the latent niche geometry
## from simulated compositional data (recovery-first), and its operator matches an
## independent brute-force variation matrix (correctness).

test_that("clr recovers latent niche geometry from simulated cover", {
  recov <- function(seed, dim = 12L) {
    df <- sim_communities(M = 800L, S = 60L, seed = seed, abundance = TRUE)
    x <- specvec(df, "plot", "species", abundance = "cover")
    emb <- species_embedding(x, method = "clr", dim = dim, min_occurrence = 3L)
    mu <- mu_for(emb$species, attr(df, "mu"))
    suppressWarnings(cor(as.numeric(dist(mu)), as.numeric(dist(emb$V)), method = "spearman"))
  }
  vals <- vapply(1:3, recov, numeric(1))
  expect_true(all(vals > 0.6),
              info = paste("clr recovery spearman:", paste(round(vals, 3), collapse = ", ")))
})

test_that(".clr_covariance matches an independent brute-force variation matrix", {
  set.seed(3)
  M <- 120L; S <- 8L
  COVm <- matrix(0, M, S)
  for (p in seq_len(M)) {
    k <- sample(2:S, 1); s <- sample(S, k); COVm[p, s] <- runif(k, 1, 100)
  }
  COV <- Matrix::Matrix(COVm, sparse = TRUE)
  species <- paste0("s", seq_len(S))

  G <- specvec:::.clr_covariance(COV, species, min_copresent = 2L)

  ## brute force: pairwise log-ratio variance over co-present plots, undefined
  ## pairs filled with the largest variation, then explicit double-centering.
  B <- (COVm > 0) * 1
  Tbf <- matrix(NA_real_, S, S)
  for (i in seq_len(S)) for (j in seq_len(S)) {
    co <- which(B[, i] > 0 & B[, j] > 0)
    if (length(co) >= 2L) {
      r <- log(COVm[co, i]) - log(COVm[co, j])
      Tbf[i, j] <- mean(r^2) - mean(r)^2
    }
  }
  diag(Tbf) <- 0
  Tbf[!is.finite(Tbf)] <- max(Tbf[is.finite(Tbf)])
  Tbf <- (Tbf + t(Tbf)) / 2
  Jc <- diag(S) - matrix(1 / S, S, S)
  Gref <- -0.5 * (Jc %*% Tbf %*% Jc)
  Gref <- (Gref + t(Gref)) / 2
  dimnames(Gref) <- list(species, species)

  expect_equal(G, Gref, tolerance = 1e-8)
  expect_true(isSymmetric(unname(G)))
  expect_equal(unname(diag(G)) >= -1e-8, rep(TRUE, S))   # self-variance non-negative
})

test_that("clr without cover falls back to presence PMI with a note", {
  df <- sim_communities(M = 150L, S = 20L, seed = 5L, abundance = FALSE)
  x <- specvec(df, "plot", "species")
  expect_message(
    emb <- species_embedding(x, method = "clr", dim = 6L, min_occurrence = 2L),
    "no cover")
  pmi <- species_embedding(x, method = "pmi", dim = 6L, min_occurrence = 2L)
  expect_equal(emb$V, pmi$V)                              # same operator, same vectors
})

test_that("clr is a registered method usable through the standard verbs", {
  expect_true("clr" %in% specvec_methods())
  df <- sim_communities(M = 200L, S = 25L, seed = 2L, abundance = TRUE)
  x <- specvec(df, "plot", "species", abundance = "cover")
  emb <- species_embedding(x, method = "clr", dim = 8L, min_occurrence = 2L)
  expect_s3_class(emb, "specvec_embedding")
  expect_equal(ncol(emb$V), 8L)
  expect_true(all(is.finite(emb$V)))
  comm <- community_embedding(x, embedding = emb)        # pools the clr species vectors
  expect_equal(nrow(comm$U), length(x$plots))
})
