## Neutral benchmark metrics, transcribed from bakeoff/bakeoff.R. All scoring is
## on a fixed kept-species index and a shared train/test plot split, so every
## method is judged on identical pairs, identical pos/neg edge sets, and one
## community definition (presence-mean pooling). The diagonal is excluded.

## Scored species pairs: all upper-triangular pairs when cheap, else a
## deduplicated sample. Uses the seeded RNG already set by the benchmark loop.
.make_pairs <- function(S, control) {
  full <- S * (S - 1) / 2
  if (full <= control$pair_all_max) {
    ij <- which(upper.tri(matrix(0L, S, S)), arr.ind = TRUE)
    return(list(ip = ij[, 1], jp = ij[, 2]))
  }
  n <- ceiling(control$pair_sample * 1.4)
  i <- sample.int(S, n, TRUE); j <- sample.int(S, n, TRUE)
  k <- i < j; i <- i[k]; j <- j[k]
  d <- !duplicated(i * (S + 1) + j); i <- i[d]; j <- j[d]
  m <- min(control$pair_sample, length(i))
  list(ip = i[seq_len(m)], jp = j[seq_len(m)])
}

## Co-occurrence metrics: Spearman of the embedding gram score against held-out
## raw counts (neutral) and held-out PPMI (home-field, reference), plus the
## rank-based link-prediction AUC over never-seen negatives.
.co_metrics <- function(emb, ev) {
  g <- rowSums(emb[ev$ip, , drop = FALSE] * emb[ev$jp, , drop = FALSE])
  ppmi <- suppressWarnings(stats::cor(g, ev$test_ppmi, method = "spearman"))
  raw  <- suppressWarnings(stats::cor(g, ev$test_raw,  method = "spearman"))
  if (length(ev$pos) >= 20L && length(ev$neg) >= 20L) {
    sc <- c(g[ev$pos], g[ev$neg]); rk <- rank(sc)
    np <- as.numeric(length(ev$pos)); nn <- as.numeric(length(ev$neg))
    auc <- (sum(rk[seq_len(np)]) - np * (np + 1) / 2) / (np * nn)
  } else auc <- NA_real_
  c(cooc_ppmi = ppmi, cooc_raw = raw, link_auc = auc)
}

## Held-out co-occurrence evidence for one seed's split, shared by every method.
.make_eval <- function(BIN, tr_plot, te_plot, control) {
  Ctr <- methods::as(Matrix::crossprod(BIN[tr_plot, , drop = FALSE]), "CsparseMatrix")
  Cte <- methods::as(Matrix::crossprod(BIN[te_plot, , drop = FALSE]), "CsparseMatrix")
  pr  <- .make_pairs(ncol(BIN), control)
  fte <- Matrix::diag(Cte)
  te_raw <- Cte[cbind(pr$ip, pr$jp)]; tr_raw <- Ctr[cbind(pr$ip, pr$jp)]
  den <- fte[pr$ip] * fte[pr$jp]
  te_ppmi <- pmax(ifelse(te_raw > 0 & den > 0, log(te_raw * length(te_plot) / den), 0), 0)
  pos <- which(te_raw > 0); neg <- which(te_raw == 0 & tr_raw == 0)
  if (length(pos) > control$auc_cap) pos <- sample(pos, control$auc_cap)
  if (length(neg) > control$auc_cap) neg <- sample(neg, control$auc_cap)
  list(ip = pr$ip, jp = pr$jp, test_ppmi = te_ppmi, test_raw = te_raw, pos = pos, neg = neg)
}

## kNN over Euclidean distance. FNN (Suggests) when present; a dependency-light
## brute-force fallback otherwise, bounded by the EUNIS reference/query caps.
.knn_index <- function(ref, query, k) {
  k <- min(as.integer(k), nrow(ref))
  if (requireNamespace("FNN", quietly = TRUE))
    return(FNN::knnx.index(ref, query, k = k))
  D <- .pairwise_euclid(query, ref)
  t(apply(D, 1L, function(row) order(row)[seq_len(k)]))
}

.knn_vote <- function(ref, query, y, k) {
  idx <- .knn_index(ref, query, k)
  apply(idx, 1L, function(ix) { tb <- table(y[ix]); names(tb)[which.max(tb)] })
}

.knn_reg <- function(ref, query, y, k) {
  idx <- .knn_index(ref, query, k)
  apply(idx, 1L, function(ix) mean(y[ix]))
}

## Standardize a train/test feature block by the training mean/sd.
.scale_tt <- function(tr, te) {
  m <- colMeans(tr); s <- apply(tr, 2L, stats::sd)
  s[!is.finite(s) | s == 0] <- 1
  list(tr = sweep(sweep(tr, 2L, m), 2L, s, "/"),
       te = sweep(sweep(te, 2L, m), 2L, s, "/"))
}

.macro_f1 <- function(pred, truth) {
  pred <- as.character(pred); truth <- as.character(truth)
  classes <- sort(unique(truth))
  mean(vapply(classes, function(k) {
    tp <- sum(pred == k & truth == k)
    fp <- sum(pred == k & truth != k)
    fn <- sum(pred != k & truth == k)
    p <- if (tp + fp == 0) 0 else tp / (tp + fp)
    r <- if (tp + fn == 0) 0 else tp / (tp + fn)
    if (p + r == 0) 0 else 2 * p * r / (p + r)
  }, numeric(1)))
}

## EUNIS habitat recovery: community embeddings (presence-mean pooling) -> kNN
## vote, macro-F1 + accuracy against the majority baseline.
.eunis_metric <- function(emb, BIN, eun_map, tr_plot, te_plot, control) {
  na_out <- c(eunis_f1 = NA_real_, eunis_acc = NA_real_, eunis_base = NA_real_)
  trp <- tr_plot[!is.na(eun_map[tr_plot])]
  tep <- te_plot[!is.na(eun_map[te_plot])]
  if (length(trp) < control$knn_k + 1L || length(tep) < 1L) return(na_out)
  if (length(trp) > control$eunis_ref_cap) trp <- sample(trp, control$eunis_ref_cap)
  if (length(tep) > control$eunis_q_cap)   tep <- sample(tep, control$eunis_q_cap)
  ytr <- eun_map[trp]
  tep <- tep[eun_map[tep] %in% unique(ytr)]; yte <- eun_map[tep]
  if (length(tep) < 1L) return(na_out)
  sc <- .scale_tt(.pool_rows(BIN[trp, , drop = FALSE], emb),
                  .pool_rows(BIN[tep, , drop = FALSE], emb))
  pred <- .knn_vote(sc$tr, sc$te, ytr, control$knn_k)
  maj <- names(sort(table(ytr), decreasing = TRUE))[1]
  c(eunis_f1 = .macro_f1(pred, yte),
    eunis_acc = mean(as.character(pred) == yte),
    eunis_base = mean(yte == maj))
}

## Trait recovery: species embeddings -> kNN regression R^2, averaged over traits.
.trait_metric <- function(emb, TRAITM, sp_tr, sp_te, control) {
  r2 <- numeric(0)
  for (t in colnames(TRAITM)) {
    y <- TRAITM[, t]
    itr <- sp_tr[is.finite(y[sp_tr])]; ite <- sp_te[is.finite(y[sp_te])]
    if (length(itr) < control$knn_k + 5L || length(ite) < 5L) next
    sc <- .scale_tt(emb[itr, , drop = FALSE], emb[ite, , drop = FALSE])
    pred <- .knn_reg(sc$tr, sc$te, y[itr], control$knn_k)
    sse <- sum((y[ite] - pred)^2); sst <- sum((y[ite] - mean(y[itr]))^2)
    r2 <- c(r2, if (sst > 0) 1 - sse / sst else NA_real_)
  }
  if (length(r2) == 0L) NA_real_ else mean(r2, na.rm = TRUE)
}
