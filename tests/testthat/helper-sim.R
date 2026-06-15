## Simulate communities from known latent species niches. Species s sits at a
## latent position mu_s; plot p has environment env_p; presence probability
## decays with squared latent distance, so species with nearby niches co-occur
## more. The known mu is the ground truth a recovery test checks against.
sim_communities <- function(M = 600L, S = 50L, L = 2L, seed = 1L,
                            abundance = FALSE, width = 1.0, spread = 1.6) {
  set.seed(seed)
  mu  <- matrix(stats::rnorm(S * L, sd = spread), S, L)
  env <- matrix(stats::rnorm(M * L, sd = spread), M, L)
  rows <- vector("list", M)
  for (p in seq_len(M)) {
    d2 <- rowSums((mu - matrix(env[p, ], S, L, byrow = TRUE))^2)
    prob <- exp(-d2 / (2 * width^2))
    present <- which(stats::runif(S) < prob)
    if (length(present) < 2L) present <- order(prob, decreasing = TRUE)[1:2]
    cov <- if (abundance) round(100 * prob[present] / max(prob[present]), 1) else NA_real_
    rows[[p]] <- data.frame(plot = paste0("p", p),
                            species = paste0("s", present),
                            cover = cov,
                            stringsAsFactors = FALSE)
  }
  df <- do.call(rbind, rows)
  attr(df, "mu") <- mu
  df
}

## Map kept species ("s12") back to their latent rows in mu.
mu_for <- function(species, mu) {
  idx <- as.integer(sub("^s", "", species))
  mu[idx, , drop = FALSE]
}

## Simulate a focal species whose co-occurrence drifts between two well-separated
## frame clusters across time windows. Cluster A species ("a1".."aK") sit on the
## left, cluster B ("b1".."bK") on the right. In window w, a fraction (w-1)/(W-1)
## of plots are B-affiliated and the rest A-affiliated; the focal "X" occurs in
## every plot, so its co-occurrents shift A -> B over time. A correct fixed-frame
## projection puts X near the A centroid early and the B centroid late. Returns a
## long data frame with a `decade` time column and cluster membership attributes.
sim_drift <- function(n_windows = 4L, plots_per_window = 200L, per_cluster = 10L,
                      seed = 1L, abundance = TRUE, sep = 6, width = 1.0) {
  set.seed(seed)
  K <- per_cluster
  muA <- cbind(stats::rnorm(K, -sep / 2, 0.4), stats::rnorm(K, 0, 0.4))
  muB <- cbind(stats::rnorm(K, +sep / 2, 0.4), stats::rnorm(K, 0, 0.4))
  mu  <- rbind(muA, muB)
  nm  <- c(paste0("a", seq_len(K)), paste0("b", seq_len(K)))
  S <- 2L * K
  rows <- list(); pid <- 0L
  for (w in seq_len(n_windows)) {
    fB <- if (n_windows > 1L) (w - 1) / (n_windows - 1) else 0
    for (i in seq_len(plots_per_window)) {
      pid <- pid + 1L
      center <- if (stats::runif(1) < fB) c(sep / 2, 0) else c(-sep / 2, 0)
      d2 <- rowSums((mu - matrix(center, S, 2L, byrow = TRUE))^2)
      prob <- exp(-d2 / (2 * width^2))
      present <- which(stats::runif(S) < prob)
      if (length(present) < 2L) present <- order(prob, decreasing = TRUE)[1:3]
      sp  <- c(nm[present], "X")
      cov <- if (abundance) round(100 * c(prob[present], max(prob[present])) /
                                    max(prob[present]), 1) else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        plot = paste0("p", pid), species = sp, cover = cov,
        decade = 1990L + 10L * (w - 1L), stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)
  attr(df, "frame")    <- nm
  attr(df, "clusterA") <- paste0("a", seq_len(K))
  attr(df, "clusterB") <- paste0("b", seq_len(K))
  df
}

## The alien-integration scenario (the sim_drift family, recast for v0.3). A
## native community (cluster "n1".."nK", right pole) is present in every window
## as a stable backbone, so the native centroid stays put through time. A
## disturbed/ruderal pool ("d1".."dK", left pole) is the flora the focal arrives
## with. The focal "X" rides disturbed plots early and native plots late: a
## fraction (w-1)/(W-1) of its plots are native-centered in window w, so its
## co-occurrents drift D -> N. A correct integration_trajectory shows the focal's
## distance to the native centroid falling monotonically toward zero. Attributes
## carry the frame and the native/disturbed memberships.
sim_integration <- function(n_windows = 4L, plots_per_window = 200L,
                            per_cluster = 10L, seed = 1L, abundance = TRUE,
                            sep = 6, width = 1.0, native_frac = 0.5) {
  set.seed(seed)
  K <- per_cluster
  muN <- cbind(stats::rnorm(K, +sep / 2, 0.4), stats::rnorm(K, 0, 0.4))
  muD <- cbind(stats::rnorm(K, -sep / 2, 0.4), stats::rnorm(K, 0, 0.4))
  mu  <- rbind(muN, muD)
  nm  <- c(paste0("n", seq_len(K)), paste0("d", seq_len(K)))
  S <- 2L * K
  cN <- c(+sep / 2, 0); cD <- c(-sep / 2, 0)
  gen <- function(center, focal) {
    d2 <- rowSums((mu - matrix(center, S, 2L, byrow = TRUE))^2)
    prob <- exp(-d2 / (2 * width^2))
    present <- which(stats::runif(S) < prob)
    if (length(present) < 2L) present <- order(prob, decreasing = TRUE)[1:3]
    sp  <- nm[present]
    cov <- if (abundance) round(100 * prob[present] / max(prob[present]), 1) else NA_real_
    if (focal) { sp <- c(sp, "X"); if (abundance) cov <- c(cov, max(cov, 100)) }
    list(sp = sp, cov = cov)
  }
  add <- function(rows, pid, dec, g)
    c(rows, list(data.frame(plot = paste0("p", pid), species = g$sp, cover = g$cov,
                            decade = dec, stringsAsFactors = FALSE)))
  rows <- list(); pid <- 0L
  for (w in seq_len(n_windows)) {
    fXn <- if (n_windows > 1L) (w - 1) / (n_windows - 1) else 1
    dec <- 1990L + 10L * (w - 1L)
    n_back <- round(plots_per_window * native_frac)   # native backbone (no focal)
    n_foc  <- plots_per_window - n_back               # focal-carrying plots
    for (i in seq_len(n_back)) { pid <- pid + 1L; rows <- add(rows, pid, dec, gen(cN, FALSE)) }
    for (i in seq_len(n_foc))  {
      pid <- pid + 1L
      center <- if (stats::runif(1) < fXn) cN else cD
      rows <- add(rows, pid, dec, gen(center, TRUE))
    }
  }
  df <- do.call(rbind, rows)
  attr(df, "frame")     <- nm
  attr(df, "native")    <- paste0("n", seq_len(K))
  attr(df, "disturbed") <- paste0("d", seq_len(K))
  df
}
