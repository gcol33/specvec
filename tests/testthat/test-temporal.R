## Temporal layer (v0.2): windowed embeddings, the fixed-frame trajectory
## projection, and novelty over time. The trajectory test is recovery-first --
## it checks that a focal species whose co-occurrence drifts between two frame
## clusters actually moves from one to the other in the fixed frame.

test_that("time window selects exactly the matching plots", {
  df <- data.frame(
    plot    = c("p1","p1","p2","p2","p3","p3","p4","p4"),
    species = c("A","B","A","C","B","C","A","B"),
    decade  = c(1990,1990,1990,1990,2000,2000,2000,2000)
  )
  x <- specvec(df, "plot", "species", time = "decade")
  rows <- which(x$time == 1990)
  expect_setequal(x$plots[rows], c("p1", "p2"))

  ## windowed co-occurrence counts equal a manual crossprod on the same plots
  Cw <- cooc_matrix(x, "counts", time = 1990, min_occurrence = 1L)
  ks <- which(as.numeric(Matrix::colSums(x$P[rows, , drop = FALSE])) >= 1)
  manual <- as.matrix(Matrix::crossprod(x$P[rows, ks, drop = FALSE]))
  dimnames(manual) <- list(x$species[ks], x$species[ks])
  expect_equal(as.matrix(Cw), manual)
})

test_that("windowed species_embedding fits on the window's plots only", {
  df <- sim_drift(n_windows = 2L, plots_per_window = 60L, per_cluster = 6L, seed = 4L)
  x <- specvec(df, "plot", "species", abundance = "cover", time = "decade")
  emb <- species_embedding(x, "abund_pmi", dim = 4L, time = 1990L, min_occurrence = 2L)
  expect_equal(emb$preprocessing$n_plots, sum(x$time == 1990L))
  expect_equal(emb$preprocessing$time_window, 1990L)
})

test_that("species_embedding errors cleanly when no time is stored", {
  df <- sim_communities(M = 80L, S = 12L, seed = 2L)
  x <- specvec(df, "plot", "species")
  expect_error(species_embedding(x, "pmi", dim = 4L, time = 1L, min_occurrence = 1L),
               "no `time`")
})

test_that("fixed-frame trajectory tracks a focal species drifting A -> B", {
  ## Signed gap g(w) = dist(traj_w, centroid_A) - dist(traj_w, centroid_B).
  ## Early windows: focal near A so g < 0; late windows: near B so g > 0.
  gap_path <- function(seed) {
    df <- sim_drift(n_windows = 4L, plots_per_window = 200L, per_cluster = 10L,
                    seed = seed, abundance = TRUE, sep = 6)
    x <- specvec(df, "plot", "species", abundance = "cover", time = "decade")
    tr <- species_trajectory(x, species = "X", frame = attr(df, "frame"),
                             dim = 8L, weights = "cover", min_occurrence = 5L)
    Vf <- tr$frame$V
    cA <- colMeans(Vf[attr(df, "clusterA"), , drop = FALSE])
    cB <- colMeans(Vf[attr(df, "clusterB"), , drop = FALSE])
    vapply(seq_len(nrow(tr$windows)), function(w) {
      u <- tr$U["X", w, ]
      sqrt(sum((u - cA)^2)) - sqrt(sum((u - cB)^2))
    }, numeric(1))
  }
  for (s in 1:2) {
    g <- gap_path(s)
    expect_lt(g[1], 0)                      # first window sits with cluster A
    expect_gt(g[length(g)], 0)              # last window sits with cluster B
    expect_gt(g[length(g)], g[1])           # net drift A -> B
    info <- paste("gap path:", paste(round(g, 3), collapse = ", "))
    expect_true(g[length(g)] - g[1] > 1, info = info)
  }
})

test_that("trajectory cells with no co-occurrence are NA with zero support", {
  df <- data.frame(
    plot = c("p1","p1","p2","p2","p3","p3","p4","p4",
             "p5","p5","p6","p6","p7","p7","p8","p8"),
    species = c("A","B","A","C","B","C","A","B",
                "A","Z","B","Z","C","Z","A","Z"),
    decade = c(rep(1990L, 8), rep(2000L, 8))
  )
  x <- specvec(df, "plot", "species", time = "decade")
  tr <- species_trajectory(x, species = "Z", frame = c("A","B","C"),
                           dim = 2L, weights = "presence", min_occurrence = 1L)
  expect_equal(unname(tr$support["Z", "1990"]), 0L)
  expect_true(all(is.na(tr$U["Z", "1990", ])))
  expect_gt(tr$support["Z", "2000"], 0L)
  expect_false(any(is.na(tr$U["Z", "2000", ])))

  d <- as.data.frame(tr)
  expect_equal(nrow(d), 2L)                       # 1 focal x 2 windows
  expect_true(all(c("species","window","center","support","d1","d2") %in% names(d)))
  expect_equal(nrow(as.data.frame(tr, na.rm = TRUE)), 1L)
})

test_that("species_trajectory reuses a supplied frame embedding", {
  df <- sim_drift(n_windows = 3L, plots_per_window = 80L, per_cluster = 8L, seed = 5L)
  x <- specvec(df, "plot", "species", abundance = "cover", time = "decade")
  fr <- species_embedding(x, "abund_pmi", dim = 5L, min_occurrence = 3L)
  tr <- species_trajectory(x, species = "X", frame_embedding = fr, weights = "cover")
  expect_identical(tr$frame$V, fr$V)               # frame reused verbatim, not refit
  expect_equal(dim(tr$U), c(1L, nrow(tr$windows), 5L))
})

test_that("community novelty over time rises as communities drift from the reference", {
  df <- sim_drift(n_windows = 4L, plots_per_window = 200L, per_cluster = 10L,
                  seed = 1L, abundance = TRUE, sep = 6)
  x <- specvec(df, "plot", "species", abundance = "cover", time = "decade")
  ct <- community_trajectory(x, dim = 8L, weights = "cover", k = 5L, min_occurrence = 5L)
  nv <- ct$novelty
  expect_equal(nrow(nv), 4L)
  expect_gt(nv$mean_novelty[nrow(nv)], nv$mean_novelty[1])   # later windows more novel
  expect_equal(which.min(nv$mean_novelty), 1L)              # reference window is least novel
  expect_identical(ct$reference, nv$window[1])
})

test_that("by-breaks windowing bins the time axis", {
  df <- sim_drift(n_windows = 4L, plots_per_window = 40L, per_cluster = 6L, seed = 7L)
  x <- specvec(df, "plot", "species", abundance = "cover", time = "decade")
  ## decades are 1990,2000,2010,2020; two bins split at 2005
  tr <- species_trajectory(x, species = "X", frame = attr(df, "frame"),
                           by = c(1985, 2005, 2025), dim = 4L, min_occurrence = 3L)
  expect_equal(nrow(tr$windows), 2L)
  expect_equal(sum(tr$windows$n_plots), nrow(x$P))
})
