test_that("compare_embeddings returns a tidy method x metric summary", {
  df <- sim_communities(M = 600L, S = 80L, seed = 1L, abundance = TRUE)
  x <- specvec(df, "plot", "species", abundance = "cover")
  b <- compare_embeddings(x, methods = c("ca", "pmi", "abund_pmi"),
                          metrics = c("cooc_ppmi", "cooc_raw", "link_auc"),
                          dim = 12L, seeds = 1:2, min_occurrence = 2L)

  expect_s3_class(b, "specvec_benchmark")
  expect_true(all(c("method", "metric", "mean", "sd", "n") %in% names(b$summary)))
  expect_setequal(unique(b$summary$method), c("ca", "pmi", "abund_pmi"))
  expect_setequal(unique(b$summary$metric), c("cooc_ppmi", "cooc_raw", "link_auc"))
  expect_equal(nrow(b$raw), 6L)                       # 3 methods x 2 seeds

  auc <- b$summary[b$summary$metric == "link_auc", ]
  expect_true(all(auc$mean >= 0 & auc$mean <= 1))
})

get_mean <- function(b, m, cc)
  b$summary$mean[b$summary$method == m & b$summary$metric == cc]

test_that("every method recovers held-out co-occurrence well above chance", {
  ## sim_communities draws presence from smooth Gaussian latent niches, the
  ## unimodal-gradient regime CA is built for, so CA leads here; the AbundPMI
  ## advantage over CA is a real-data / high-sparsity result (bakeoff.R, sec 13),
  ## not a property a clean low-dim sim adjudicates. What the sim does support is
  ## that each method reconstructs held-out co-occurrence far above chance.
  df <- sim_communities(M = 800L, S = 90L, seed = 2L, abundance = TRUE)
  x <- specvec(df, "plot", "species", abundance = "cover")
  b <- compare_embeddings(x, methods = c("ca", "pmi", "abund_pmi"),
                          metrics = c("cooc_raw", "link_auc"),
                          dim = 16L, seeds = 1:3, min_occurrence = 3L)

  for (m in c("ca", "pmi", "abund_pmi")) {
    expect_gt(get_mean(b, m, "link_auc"), 0.6)       # AUC 0.5 = chance
    expect_gt(get_mean(b, m, "cooc_raw"), 0)
  }
  expect_true(b$verdict$has_reference)
  expect_true("abund_pmi" %in% names(b$verdict$rows))
})

test_that("abundance weighting adds held-out signal over presence PMI", {
  ## The clean A/B the sim does support: AbundPMI is PMI with sqrt(COV) for
  ## presence, and cover here is a monotone niche signal, so it should not lose.
  df <- sim_communities(M = 800L, S = 90L, seed = 8L, abundance = TRUE)
  x <- specvec(df, "plot", "species", abundance = "cover")
  b <- compare_embeddings(x, methods = c("pmi", "abund_pmi"),
                          metrics = c("cooc_raw", "link_auc"),
                          dim = 16L, seeds = 1:3, min_occurrence = 3L)

  expect_gte(get_mean(b, "abund_pmi", "link_auc"), get_mean(b, "pmi", "link_auc") - 1e-3)
  expect_gte(get_mean(b, "abund_pmi", "cooc_raw"), get_mean(b, "pmi", "cooc_raw") - 1e-3)
})

test_that("benchmark is deterministic given the seeds", {
  df <- sim_communities(M = 400L, S = 50L, seed = 5L, abundance = TRUE)
  x <- specvec(df, "plot", "species", abundance = "cover")
  args <- list(x = x, methods = c("ca", "pmi", "abund_pmi"),
               metrics = c("cooc_raw", "link_auc"), dim = 8L, seeds = 1:2,
               min_occurrence = 2L)
  b1 <- do.call(compare_embeddings, args)
  b2 <- do.call(compare_embeddings, args)
  expect_equal(b1$summary, b2$summary)
  expect_equal(b1$raw, b2$raw)
})

test_that("EUNIS metric activates with plot labels and skips without", {
  df <- sim_communities(M = 600L, S = 70L, seed = 3L, abundance = TRUE)
  ## plot-level habitat tied to the dominant species' cluster
  lab <- tapply(seq_len(nrow(df)), df$plot, function(ix)
    paste0("h", as.integer(sub("^s", "", df$species[ix][which.max(df$cover[ix])])) %% 4L))
  df$habitat <- as.character(lab[df$plot])

  x_lab <- specvec(df, "plot", "species", abundance = "cover", labels = "habitat")
  b <- compare_embeddings(x_lab, methods = c("pmi", "abund_pmi"),
                          metrics = c("link_auc", "eunis"), labels = "habitat",
                          dim = 10L, seeds = 1L, min_occurrence = 2L)
  expect_true("eunis_f1" %in% b$summary$metric)
  f1 <- b$summary$mean[b$summary$metric == "eunis_f1"]
  expect_true(all(is.finite(f1) & f1 >= 0 & f1 <= 1))

  x_nolab <- specvec(df, "plot", "species", abundance = "cover")
  expect_message(
    b2 <- compare_embeddings(x_nolab, methods = "abund_pmi",
                             metrics = c("link_auc", "eunis"),
                             dim = 10L, seeds = 1L, min_occurrence = 2L),
    "skipping"
  )
  expect_false("eunis_f1" %in% b2$summary$metric)
})

test_that("trait metric activates with a species trait table", {
  df <- sim_communities(M = 700L, S = 90L, seed = 4L, abundance = TRUE)
  mu <- attr(df, "mu")
  x <- specvec(df, "plot", "species", abundance = "cover")
  sp <- sort(unique(df$species))
  idx <- as.integer(sub("^s", "", sp))
  set.seed(1)
  traits <- matrix(mu[idx, 1] + stats::rnorm(length(sp), sd = 0.3),
                   ncol = 1, dimnames = list(sp, "niche1"))

  b <- compare_embeddings(x, methods = c("pmi", "abund_pmi"),
                          metrics = c("link_auc", "trait"), traits = traits,
                          dim = 12L, seeds = 1:2, min_occurrence = 3L)
  expect_true("trait_r2" %in% b$summary$metric)
  expect_true(all(is.finite(b$summary$mean[b$summary$metric == "trait_r2"])))
})

test_that("specvec_benchmark is an alias and print renders without error", {
  expect_identical(specvec_benchmark, compare_embeddings)
  df <- sim_communities(M = 300L, S = 40L, seed = 6L, abundance = TRUE)
  x <- specvec(df, "plot", "species", abundance = "cover")
  b <- compare_embeddings(x, methods = c("ca", "abund_pmi"),
                          metrics = c("cooc_raw", "link_auc"),
                          dim = 6L, seeds = 1L, min_occurrence = 2L)
  expect_output(print(b), "specvec_benchmark")
  expect_output(print(b), "verdict")
})

test_that("glove runs through the benchmark when text2vec is available", {
  skip_if_not_installed("text2vec")
  df <- sim_communities(M = 400L, S = 50L, seed = 7L, abundance = TRUE)
  x <- specvec(df, "plot", "species", abundance = "cover")
  b <- compare_embeddings(x, methods = c("abund_pmi", "glove"),
                          metrics = c("link_auc"), dim = 8L, seeds = 1L,
                          min_occurrence = 2L, glove_iter = 8L)
  expect_true("glove" %in% b$summary$method)
})
