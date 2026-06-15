test_that("method presets expand to the right weighting and factorization", {
  df <- sim_communities(M = 200L, S = 24L, seed = 7L, abundance = TRUE)
  x <- specvec(df, "plot", "species", abundance = "cover")

  e_ca  <- species_embedding(x, "ca", dim = 6L, min_occurrence = 1L)
  e_pmi <- species_embedding(x, "pmi", dim = 6L, min_occurrence = 1L)
  e_ab  <- species_embedding(x, "abund_pmi", dim = 6L, min_occurrence = 1L)

  expect_identical(c(e_ca$weighting, e_ca$factorization), c("chi_square", "svd"))
  expect_identical(c(e_pmi$weighting, e_pmi$factorization), c("ppmi", "eigen"))
  expect_identical(c(e_ab$weighting, e_ab$factorization), c("abundance_pmi", "eigen"))

  expect_true(all(c("ca","pmi","abund_pmi","glove") %in% specvec_methods()))
})

test_that("incompatible weighting + factorization is rejected by capability check", {
  df <- sim_communities(M = 150L, S = 20L, seed = 8L)
  x <- specvec(df, "plot", "species")
  ## ppmi produces a 'sym' operator; the CA 'svd' factorizer consumes 'implicit'.
  expect_error(
    species_embedding(x, method = NULL, weighting = "ppmi", factorization = "svd",
                      dim = 4L, min_occurrence = 1L),
    "consumes operator kind"
  )
})

test_that("unknown names error with the available list", {
  df <- sim_communities(M = 120L, S = 18L, seed = 9L)
  x <- specvec(df, "plot", "species")
  expect_error(species_embedding(x, method = "nope", dim = 4L), "unknown method")
  expect_error(cooc_matrix(x, "counts", min_occurrence = 999L), "fewer than 2 species")
})

test_that("dim larger than species count is padded, not an error", {
  df <- sim_communities(M = 150L, S = 15L, seed = 10L)
  x <- specvec(df, "plot", "species")
  emb <- species_embedding(x, "pmi", dim = 40L, min_occurrence = 1L)
  expect_equal(ncol(emb$V), 40L)
  expect_true(nrow(emb$V) <= 15L)
})
