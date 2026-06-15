test_that("presence pooling equals the mean of present species vectors", {
  df <- sim_communities(M = 300L, S = 30L, seed = 3L, abundance = TRUE)
  x <- specvec(df, "plot", "species", abundance = "cover")
  emb <- species_embedding(x, "pmi", dim = 8L, min_occurrence = 2L)

  comm <- community_embedding(x, embedding = emb, weights = "presence")
  ks <- match(emb$species, x$species)
  P <- x$P[, ks, drop = FALSE]
  U_manual <- as.matrix(P %*% emb$V) / pmax(as.numeric(Matrix::rowSums(P)), 1)
  rownames(U_manual) <- x$plots
  expect_equal(comm$U, U_manual, tolerance = 1e-10)
})

test_that("cover pooling differs from presence and pulls toward the dominant species", {
  ## p1 is overwhelmingly species A. With known species vectors the pooling math
  ## is tested directly: cover pooling should sit nearer A than the presence mean.
  df <- data.frame(
    plot = c("p1","p1","p2","p2","p3","p3"),
    species = c("A","B","A","C","B","C"),
    cover = c(99, 1, 50, 50, 50, 50),
    stringsAsFactors = FALSE
  )
  x <- specvec(df, "plot", "species", abundance = "cover")
  V <- rbind(A = c(1, 0), B = c(0, 1), C = c(-1, 0))
  emb <- structure(list(V = V, species = rownames(V), method = "pmi", dim = 2L,
                        weighting = "ppmi", factorization = "eigen"),
                   class = "specvec_embedding")

  cc <- community_embedding(x, embedding = emb, weights = "cover")
  cp <- community_embedding(x, embedding = emb, weights = "presence")
  expect_false(isTRUE(all.equal(cc$U, cp$U)))

  cos <- function(u, v) sum(u * v) / sqrt(sum(u^2) * sum(v^2))
  expect_gt(cos(cc$U["p1", ], V["A", ]), cos(cp$U["p1", ], V["A", ]))
})

test_that("community pooling falls back to presence without cover", {
  df <- sim_communities(M = 200L, S = 20L, seed = 7L, abundance = FALSE)
  x <- specvec(df, "plot", "species")
  expect_message(
    comm <- community_embedding(x, method = "pmi", dim = 6L, min_occurrence = 1L),
    "uses presence"
  )
  expect_identical(comm$pooling, "presence")
})

test_that("community_novelty scores far communities above near ones", {
  set.seed(11)
  R <- matrix(stats::rnorm(120 * 4), 120, 4)
  U <- rbind(R[1, ] + 1e-3, R[1, ] + 50)
  rownames(U) <- c("near", "far")
  obj <- structure(list(U = U), class = "specvec_community")
  nov <- community_novelty(obj, R, k = 5L)
  expect_lt(nov["near"], nov["far"])
  expect_named(nov, c("near", "far"))
})

test_that("community_similarity self-matrix has unit diagonal (cosine)", {
  df <- sim_communities(M = 150L, S = 18L, seed = 9L, abundance = TRUE)
  x <- specvec(df, "plot", "species", abundance = "cover")
  comm <- community_embedding(x, method = "abund_pmi", dim = 6L, min_occurrence = 1L)
  S <- community_similarity(comm)
  expect_equal(dim(S), c(nrow(comm$U), nrow(comm$U)))
  expect_equal(unname(diag(S)), rep(1, nrow(S)), tolerance = 1e-8)
})
