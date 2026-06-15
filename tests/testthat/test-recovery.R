## Recovery: species with nearby latent niches co-occur more, so their embedding
## vectors should be closer. We check that pairwise embedding distance tracks
## pairwise latent distance (Spearman), across seeds. This tests the method, not
## just the plumbing.
recovery_spearman <- function(method, seed, abundance = FALSE, dim = 12L) {
  df <- sim_communities(M = 800L, S = 60L, seed = seed, abundance = abundance)
  x <- specvec(df, "plot", "species",
               abundance = if (abundance) "cover" else NULL)
  emb <- species_embedding(x, method = method, dim = dim, min_occurrence = 3L)
  mu <- mu_for(emb$species, attr(df, "mu"))
  d_latent <- as.numeric(dist(mu))
  d_emb    <- as.numeric(dist(emb$V))
  suppressWarnings(cor(d_latent, d_emb, method = "spearman"))
}

test_that("PMI recovers latent niche geometry across seeds", {
  vals <- vapply(1:3, function(s) recovery_spearman("pmi", s), numeric(1))
  expect_true(all(vals > 0.45),
              info = paste("recovery spearman:", paste(round(vals, 3), collapse = ", ")))
})

test_that("AbundPMI recovers latent niche geometry with simulated cover", {
  vals <- vapply(1:3, function(s) recovery_spearman("abund_pmi", s, abundance = TRUE),
                 numeric(1))
  expect_true(all(vals > 0.45),
              info = paste("recovery spearman:", paste(round(vals, 3), collapse = ", ")))
})

test_that("CA also recovers the latent geometry", {
  vals <- vapply(1:3, function(s) recovery_spearman("ca", s), numeric(1))
  expect_true(all(vals > 0.40),
              info = paste("recovery spearman:", paste(round(vals, 3), collapse = ", ")))
})
