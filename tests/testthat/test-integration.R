## Alien integration (v0.3): the focal neophyte's distance to the native-community
## centroid over time, in one fixed frame. The primary test is recovery-first --
## a focal whose co-occurrence drifts from a disturbed pool into a stable native
## community must show that distance falling monotonically toward zero.

test_that("integration distance to the native centroid falls as the focal drifts in", {
  dist_path <- function(seed) {
    df <- sim_integration(n_windows = 4L, plots_per_window = 200L, per_cluster = 10L,
                          seed = seed, abundance = TRUE, sep = 6)
    x <- specvec(df, "plot", "species", abundance = "cover", time = "decade")
    it <- integration_trajectory(x, species = "X", native = attr(df, "native"),
                                 dim = 8L, weights = "cover", min_occurrence = 5L)
    it$distance["X", ]
  }
  for (s in 1:3) {
    d <- dist_path(s)
    info <- paste("distance path:", paste(round(d, 3), collapse = ", "))
    expect_true(all(diff(d) < 0), info = info)          # strictly decreasing
    expect_equal(unname(which.max(d)), 1L, info = info)  # starts farthest from natives
    expect_equal(unname(which.min(d)), length(d), info = info)  # ends closest
    expect_gt(d[1] - d[length(d)], 1.0)                  # large net integration
    expect_lt(d[length(d)], 0.1)                         # effectively integrated
  }
})

test_that("cosine integration distance also falls monotonically", {
  df <- sim_integration(n_windows = 4L, plots_per_window = 200L, per_cluster = 10L,
                        seed = 2L, abundance = TRUE, sep = 6)
  x <- specvec(df, "plot", "species", abundance = "cover", time = "decade")
  it <- integration_trajectory(x, species = "X", native = attr(df, "native"),
                               dim = 8L, weights = "cover", metric = "cosine",
                               min_occurrence = 5L)
  d <- it$distance["X", ]
  expect_identical(it$metric, "cosine")
  expect_true(all(diff(d) < 0), info = paste(round(d, 4), collapse = ", "))
})

test_that("a focal absent in a window has NA distance and zero support", {
  df <- data.frame(
    plot = c("p1","p1","p2","p2","p3","p3","p4","p4",               # 1990: natives only
             "p5","p5","p5","p6","p6","p6","p7","p7","p7","p8","p8","p8"),  # 2000: Z + natives
    species = c("A","B","A","C","B","C","A","B",
                "A","B","Z","B","C","Z","A","C","Z","A","B","Z"),
    decade = c(rep(1990L, 8), rep(2000L, 12)))
  x <- specvec(df, "plot", "species", time = "decade")
  it <- integration_trajectory(x, species = "Z", native = c("A","B","C"),
                               dim = 2L, weights = "presence", min_occurrence = 1L)
  expect_equal(unname(it$support["Z", "1990"]), 0L)
  expect_true(is.na(it$distance["Z", "1990"]))
  expect_gt(it$native_support[["1990"]], 0L)                # native community exists
  expect_gt(it$support["Z", "2000"], 0L)
  expect_false(is.na(it$distance["Z", "2000"]))

  d <- as.data.frame(it)
  expect_equal(nrow(d), 2L)                                 # 1 focal x 2 windows
  expect_true(all(c("species","window","center","n_plots","support",
                    "native_support","distance") %in% names(d)))
  expect_equal(nrow(as.data.frame(it, na.rm = TRUE)), 1L)
})

test_that("a window with no native community yields NA distance even when the focal is present", {
  df <- data.frame(
    plot = c("p1","p1","p2","p2",                           # 1990: focal + ruderal, no natives
             "p3","p3","p4","p4","p5","p5","p6","p6"),       # 2000: natives (+ focal)
    species = c("Z","R","Z","R",
                "A","B","B","C","A","C","Z","A"),
    decade = c(rep(1990L, 4), rep(2000L, 8)))
  x <- specvec(df, "plot", "species", time = "decade")
  it <- integration_trajectory(x, species = "Z", native = c("A","B","C"),
                               dim = 2L, weights = "presence", min_occurrence = 1L)
  expect_equal(it$native_support[["1990"]], 0L)             # no native plots this window
  expect_gt(it$support["Z", "1990"], 0L)                    # focal itself is present
  expect_true(is.na(it$distance["Z", "1990"]))              # but no native target
  expect_gt(it$native_support[["2000"]], 0L)
  expect_false(is.na(it$distance["Z", "2000"]))
})

test_that("native species absent from the fitted frame error cleanly", {
  df <- sim_integration(n_windows = 3L, plots_per_window = 60L, per_cluster = 6L, seed = 3L)
  x <- specvec(df, "plot", "species", abundance = "cover", time = "decade")
  expect_error(
    integration_trajectory(x, species = "X", native = c("ghost1","ghost2"),
                           dim = 4L, weights = "cover", min_occurrence = 3L),
    "no `native` species are in the fitted frame")
})

test_that("integration_trajectory reuses a supplied frame embedding and composes the trajectory", {
  df <- sim_integration(n_windows = 3L, plots_per_window = 80L, per_cluster = 8L, seed = 5L)
  x <- specvec(df, "plot", "species", abundance = "cover", time = "decade")
  fr <- species_embedding(x, "abund_pmi", dim = 5L, min_occurrence = 3L)
  it <- integration_trajectory(x, species = "X", native = attr(df, "native"),
                               frame_embedding = fr, weights = "cover")
  expect_identical(it$frame$V, fr$V)                        # frame reused verbatim, not refit
  expect_identical(it$trajectory$frame$V, fr$V)             # distance read out from this frame
  expect_equal(dim(it$distance), c(1L, nrow(it$windows)))
  expect_s3_class(it$trajectory, "specvec_trajectory")
})

test_that("multiple focal species are traced against the same native centroid", {
  df <- sim_integration(n_windows = 3L, plots_per_window = 80L, per_cluster = 8L, seed = 6L)
  x <- specvec(df, "plot", "species", abundance = "cover", time = "decade")
  it <- integration_trajectory(x, species = c("X", "d1"), native = attr(df, "native"),
                               dim = 6L, weights = "cover", min_occurrence = 3L)
  expect_equal(nrow(it$distance), 2L)
  expect_setequal(rownames(it$distance), c("X", "d1"))
  d <- as.data.frame(it)
  expect_equal(nrow(d), 2L * nrow(it$windows))
})
