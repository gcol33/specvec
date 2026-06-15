test_that("nearest_species ranks niche-neighbours and excludes self", {
  df <- sim_communities(M = 800L, S = 60L, seed = 2L)
  x <- specvec(df, "plot", "species")
  emb <- species_embedding(x, "pmi", dim = 12L, min_occurrence = 3L)
  mu <- attr(df, "mu")

  focal <- emb$species[1]
  nn <- nearest_species(emb, focal, n = 5L)
  expect_false(focal %in% nn$species)             # self excluded
  expect_equal(nrow(nn), 5L)
  expect_true(all(diff(nn$similarity) <= 1e-9))   # ranked descending

  ## Embedding neighbours should overlap latent neighbours far above chance.
  kept <- emb$species
  mu_k <- mu_for(kept, mu)
  fi <- match(focal, kept)
  latent_d <- sqrt(rowSums((mu_k - matrix(mu_k[fi, ], nrow(mu_k), ncol(mu_k), byrow = TRUE))^2))
  latent_near <- kept[order(latent_d)[2:11]]      # 10 nearest in latent space
  emb_near <- nearest_species(emb, focal, n = 10L)$species
  overlap <- length(intersect(emb_near, latent_near))
  expect_gt(overlap, 3L)                           # chance ~10*10/59 < 2
})

test_that("species_similarity is symmetric and self-similarity is 1 (cosine)", {
  df <- sim_communities(M = 400L, S = 30L, seed = 4L)
  x <- specvec(df, "plot", "species")
  emb <- species_embedding(x, "pmi", dim = 8L, min_occurrence = 2L)
  a <- emb$species[1]; b <- emb$species[2]

  expect_equal(species_similarity(emb, a, b), species_similarity(emb, b, a))
  expect_equal(species_similarity(emb, a, a), 1, tolerance = 1e-8)

  full <- species_similarity(emb, a)
  expect_equal(length(full), nrow(emb$V))
  expect_equal(unname(full[a]), 1, tolerance = 1e-8)
  expect_equal(unname(full[b]), species_similarity(emb, a, b), tolerance = 1e-8)
})

test_that("euclidean nearest returns a distance column ranked ascending", {
  df <- sim_communities(M = 300L, S = 24L, seed = 6L)
  x <- specvec(df, "plot", "species")
  emb <- species_embedding(x, "pmi", dim = 6L, min_occurrence = 2L)
  nn <- nearest_species(emb, emb$species[1], n = 4L, metric = "euclidean")
  expect_true("distance" %in% names(nn))
  expect_true(all(diff(nn$distance) >= -1e-9))
})
