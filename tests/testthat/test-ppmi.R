test_that("PPMI matches a hand-computed example and zeroes the diagonal", {
  ## 4 plots: A,B co-occur in 3 plots; C,D co-occur in 1.
  ## occ(A)=occ(B)=3, occ(C)=occ(D)=1, n_plots=4.
  ## PPMI(A,B) = max(log(3 * 4 / (3 * 3)), 0) = log(12/9).
  df <- data.frame(
    plot = c("p1","p1","p2","p2","p3","p3","p4","p4"),
    species = c("A","B","A","B","A","B","C","D"),
    stringsAsFactors = FALSE
  )
  x <- specvec(df, "plot", "species")
  M <- cooc_matrix(x, "ppmi", min_occurrence = 1L)

  expect_equal(rownames(M), c("A","B","C","D"))
  expect_equal(M["A","B"], log(12 / 9), tolerance = 1e-8)
  expect_equal(M["C","D"], log(1 * 4 / (1 * 1)), tolerance = 1e-8)
  expect_true(all(Matrix::diag(M) == 0))
  expect_equal(M["A","B"], M["B","A"])  # symmetric
})

test_that("AbundPMI uses the geometric mean of covers", {
  ## Two species in two plots; AbundPMI off-diagonal A[a,b] = sum_p sqrt(cov_a cov_b).
  df <- data.frame(
    plot = c("p1","p1","p2","p2"),
    species = c("A","B","A","B"),
    cover = c(100, 100, 25, 25),  # proportions 1.0 and 0.25
    stringsAsFactors = FALSE
  )
  x <- specvec(df, "plot", "species", abundance = "cover")
  ## A[A,B] = sqrt(1*1) + sqrt(.25*.25) = 1 + .25 = 1.25
  ## f[A] = f[B] = 1 + .25 = 1.25 ; n_plots = 2
  ## PPMI = max(log(1.25 * 2 / (1.25 * 1.25)), 0) = log(2/1.25)
  M <- cooc_matrix(x, "abundance_pmi", min_occurrence = 1L)
  expect_equal(M["A","B"], log(2 / 1.25), tolerance = 1e-8)
})

test_that("abund_pmi falls back to presence PMI when no abundance is present", {
  df <- data.frame(
    plot = c("p1","p1","p2","p2","p3","p3","p4","p4"),
    species = c("A","B","A","B","A","B","C","D"),
    stringsAsFactors = FALSE
  )
  x <- specvec(df, "plot", "species")  # no abundance
  expect_message(
    emb <- species_embedding(x, method = "abund_pmi", dim = 2L, min_occurrence = 1L),
    "falls back to presence PMI"
  )
  pmi <- species_embedding(x, method = "pmi", dim = 2L, min_occurrence = 1L)
  ## Gram matrices agree: the fallback IS presence PMI.
  expect_equal(emb$V %*% t(emb$V), pmi$V %*% t(pmi$V), tolerance = 1e-6)
})
