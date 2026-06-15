# SpecVec bake-off v0 -- presence/absence, 10k-plot EVA sample
# =============================================================
# WIN CRITERION (fixed BEFORE looking at any result):
#   Primary   : held-out co-occurrence Spearman (species-vector gram vs test PPMI)
#   Secondary : EUNIS-lvl2 macro-F1 (community embedding -> kNN classify)
#               mean trait R^2     (species embedding   -> kNN regression)
#   A method "beats CA" iff, across seeds, it exceeds CA's mean on the PRIMARY
#   metric by more than 2 pooled SDs, OR exceeds CA on >=2 of the 3 metrics each
#   by more than 2 pooled SDs.
#   If nothing clears that bar -> verdict: "no advantage over CA at presence/
#   absence; SpecVec's case rests on abundance + community tooling + scale."
#   The result is reported regardless of direction.
# Fairness notes:
#   - All methods produce 64-dim species vectors; CA, PMI-SVD, GloVe judged identically.
#   - Community embedding defined identically for every method = mean of present
#     species vectors (so CA gets no special plot-score advantage).
#   - Same train/test plot split and same trait species split reused across methods
#     within a seed.
#   - Diagonal (self co-occurrence) zeroed everywhere.

suppressWarnings(suppressMessages({
  library(data.table)
  library(Matrix)
  library(RSpectra)
  library(text2vec)
  library(FNN)
}))
try(lgr::get_logger("text2vec")$set_threshold("warn"), silent = TRUE)
try(lgr::threshold("warn"), silent = TRUE)

## ---- constants -------------------------------------------------------------
DIM       <- 64L
MIN_OCC   <- 10L      # keep species occurring in >= MIN_OCC plots (global)
KNN_K     <- 15L
SEEDS     <- 1:3
TEST_FRAC <- 0.20
GLOVE_ITER<- 20L
TRAITS    <- c("SLA","PlantHeight","SeedMass","LDMC")
DATA      <- "J:/Phd Local/Gilles_paper_resolve/data"
SP_FILE   <- file.path(DATA, "species_preprocessed_sample10000.csv")
HD_FILE   <- file.path(DATA, "header_preprocessed_sample10000.csv")

cat(sprintf("[cfg] DIM=%d MIN_OCC=%d KNN_K=%d seeds=%s test_frac=%.2f glove_iter=%d\n",
            DIM, MIN_OCC, KNN_K, paste(SEEDS, collapse=","), TEST_FRAC, GLOVE_ITER))

## ---- load ------------------------------------------------------------------
sp <- fread(SP_FILE, showProgress = FALSE)
setnames(sp, names(sp), gsub("^'|'$", "", names(sp)))          # strip quoted trait names
hd <- fread(HD_FILE, select = c("PlotObservationID","Eunis_lvl2"), showProgress = FALSE)

stopifnot(all(c("PlotObservationID","WFO_TAXON","WFO_TAXON_RANK","Cover %") %in% names(sp)))
setnames(sp, c("PlotObservationID","WFO_TAXON","Cover %"), c("plot","species","cover"))

cat(sprintf("[load] species rows=%d  header rows=%d  raw distinct plots=%d  raw distinct species=%d\n",
            nrow(sp), nrow(hd), uniqueN(sp$plot), uniqueN(sp$species)))

## ---- filter ----------------------------------------------------------------
n0 <- nrow(sp)
sp <- sp[WFO_TAXON_RANK == "species" & !is.na(species) & species != ""]
if ("EVA_TAXON_GROUP" %in% names(sp)) sp <- sp[EVA_TAXON_GROUP == "Vascular plant"]
cat(sprintf("[filter] species-rank + vascular: %d -> %d records\n", n0, nrow(sp)))

freq <- sp[, .(f = uniqueN(plot)), by = species]
keep_sp <- freq[f >= MIN_OCC, species]
sp <- sp[species %in% keep_sp]
cat(sprintf("[filter] min_occ>=%d: kept %d / %d species (dropped %d rare)\n",
            MIN_OCC, length(keep_sp), nrow(freq), nrow(freq) - length(keep_sp)))

## ---- global index + matrices ----------------------------------------------
plots   <- sort(unique(sp$plot))
species <- sort(unique(sp$species))
pi_     <- match(sp$plot,    plots)
si_     <- match(sp$species, species)
M <- length(plots); S <- length(species)
cat(sprintf("[index] plots=%d  species=%d\n", M, S))
stopifnot(M > 500, S >= 200, S > DIM + 5)

# presence (dedup plot x species) and cover (max cover per plot x species, in [0,1])
pa  <- unique(data.table(p = pi_, s = si_))
BIN <- sparseMatrix(i = pa$p, j = pa$s, x = 1, dims = c(M, S))
covdt <- sp[, .(cov = suppressWarnings(max(cover, na.rm = TRUE))), by = .(p = pi_, s = si_)]
covdt[!is.finite(cov), cov := NA_real_]
med_cov <- median(covdt$cov, na.rm = TRUE); if (!is.finite(med_cov)) med_cov <- 1
covdt[is.na(cov), cov := med_cov]
COV <- sparseMatrix(i = covdt$p, j = covdt$s, x = pmin(pmax(covdt$cov,0),100)/100, dims = c(M, S))

# eunis per plot (aligned to `plots`); traits per species (aligned to `species`)
eun_map <- hd[match(plots, PlotObservationID), Eunis_lvl2]
eun_map[eun_map %in% c("", "~", "NA")] <- NA
tr <- sp[, lapply(.SD, function(x){ v <- suppressWarnings(as.numeric(x)); m <- median(v[is.finite(v)]); if (is.finite(m)) m else NA_real_ }),
         by = species, .SDcols = TRAITS]
setkey(tr, species)
TRAITM <- as.matrix(tr[species, ..TRAITS])               # S x nTrait, aligned to `species`
cat(sprintf("[aux] plots with EUNIS=%d/%d (classes=%d); species with all traits=%d\n",
            sum(!is.na(eun_map)), M, length(unique(na.omit(eun_map))),
            sum(rowSums(is.finite(TRAITM)) == length(TRAITS))))

## ---- helpers ---------------------------------------------------------------
ppmi_from_C <- function(C, nplots) {                     # C: species x species co-occ counts
  C <- as(C, "TsparseMatrix")
  i <- C@i + 1L; j <- C@j + 1L; x <- C@x
  f <- Matrix::diag(as(C, "CsparseMatrix"))
  off <- i != j & x > 0 & f[i] > 0 & f[j] > 0
  i <- i[off]; j <- j[off]
  pmi <- log(x[off] * nplots / (f[i] * f[j]))
  pp  <- pmax(pmi, 0)
  keep <- pp > 0                                          # drop explicit zeros
  sparseMatrix(i = i[keep], j = j[keep], x = pp[keep], dims = dim(C))
}

emb_pmisvd <- function(C, nplots) {
  PP <- ppmi_from_C(C, nplots)
  if (length(PP@x) == 0) return(matrix(0, S, DIM))
  PPd <- as.matrix(PP); PPd <- (PPd + t(PPd)) / 2         # exact symmetry
  e <- eigen(PPd, symmetric = TRUE)                        # values descending
  k <- min(DIM, sum(e$values > 1e-8))
  if (k < 1) return(matrix(0, S, DIM))
  emb <- e$vectors[, seq_len(k), drop = FALSE] %*% diag(sqrt(e$values[seq_len(k)]), k, k)
  if (k < DIM) emb <- cbind(emb, matrix(0, S, DIM - k))
  emb
}

emb_ca <- function(Btr) {                                # Btr: plots x species binary (train)
  P <- as.matrix(Btr); tot <- sum(P); P <- P / tot
  r <- rowSums(P); c <- colSums(P)
  cs <- c; cs[cs == 0] <- NA
  E <- outer(r, c)
  Sresid <- (P - E) / sqrt(E)
  Sresid[!is.finite(Sresid)] <- 0
  sv <- RSpectra::svds(Sresid, k = DIM, nu = 0, nv = DIM)
  V  <- sv$v                                             # species x DIM
  emb <- sweep(V, 1, sqrt(cs), "/")
  emb <- sweep(emb, 2, sv$d, "*")
  emb[!is.finite(emb)] <- 0
  emb
}

emb_glove <- function(C) {                               # C: species x species co-occ counts
  tcm <- as(C, "TsparseMatrix"); Matrix::diag(tcm) <- 0
  dimnames(tcm) <- list(species, species)
  gv <- GlobalVectors$new(rank = DIM, x_max = 10)
  wv <- suppressMessages(gv$fit_transform(tcm, n_iter = GLOVE_ITER,
                                          convergence_tol = 0.01, progressbar = FALSE))
  emb <- wv + t(gv$components)
  emb[species, , drop = FALSE]
}

community_emb <- function(B, emb) {                      # B: plots x species binary; emb: species x DIM
  cnt <- pmax(Matrix::rowSums(B), 1)
  as.matrix(B %*% emb) / cnt
}

scale_tt <- function(tr, te) {
  m <- colMeans(tr); s <- apply(tr, 2, sd); s[!is.finite(s) | s == 0] <- 1
  list(tr = sweep(sweep(tr, 2, m), 2, s, "/"),
       te = sweep(sweep(te, 2, m), 2, s, "/"))
}

macroF1 <- function(pred, truth) {
  pred <- as.character(pred); truth <- as.character(truth)
  cls <- sort(unique(truth))
  mean(vapply(cls, function(k) {
    tp <- sum(pred == k & truth == k); fp <- sum(pred == k & truth != k); fn <- sum(pred != k & truth == k)
    p <- if (tp + fp == 0) 0 else tp / (tp + fp)
    r <- if (tp + fn == 0) 0 else tp / (tp + fn)
    if (p + r == 0) 0 else 2 * p * r / (p + r)
  }, numeric(1)))
}

cooc_spearman <- function(emb, test_ppmi_dense) {
  G <- tcrossprod(emb)
  ut <- upper.tri(G)
  suppressWarnings(cor(G[ut], test_ppmi_dense[ut], method = "spearman"))
}

eunis_metric <- function(emb, tr_plots, te_plots) {
  ok_tr <- tr_plots[!is.na(eun_map[tr_plots])]
  ok_te <- te_plots[!is.na(eun_map[te_plots])]
  ytr <- eun_map[ok_tr]; yte <- eun_map[ok_te]
  ok_te <- ok_te[yte %in% unique(ytr)]; yte <- eun_map[ok_te]   # drop unseen classes
  ctr <- community_emb(BIN[ok_tr, , drop = FALSE], emb)
  cte <- community_emb(BIN[ok_te, , drop = FALSE], emb)
  sc <- scale_tt(ctr, cte)
  idx <- FNN::knnx.index(sc$tr, sc$te, k = KNN_K)
  pred <- apply(idx, 1, function(ix) { tb <- table(ytr[ix]); names(tb)[which.max(tb)] })
  maj <- names(sort(table(ytr), decreasing = TRUE))[1]
  list(f1 = macroF1(pred, yte), acc = mean(as.character(pred) == yte),
       base = mean(yte == maj))
}

trait_metric <- function(emb, sp_tr, sp_te) {
  r2 <- c()
  for (t in TRAITS) {
    y <- TRAITM[, t]
    itr <- sp_tr[is.finite(y[sp_tr])]; ite <- sp_te[is.finite(y[sp_te])]
    if (length(itr) < KNN_K + 5 || length(ite) < 5) next
    sc <- scale_tt(emb[itr, , drop = FALSE], emb[ite, , drop = FALSE])
    pred <- FNN::knn.reg(sc$tr, sc$te, y[itr], k = KNN_K)$pred
    sse <- sum((y[ite] - pred)^2); sst <- sum((y[ite] - mean(y[itr]))^2)
    r2 <- c(r2, if (sst > 0) 1 - sse / sst else NA_real_)
  }
  mean(r2, na.rm = TRUE)
}

## ---- run -------------------------------------------------------------------
res <- list()
for (seed in SEEDS) {
  set.seed(seed)
  te_plot <- sort(sample.int(M, floor(M * TEST_FRAC)))
  tr_plot <- setdiff(seq_len(M), te_plot)
  sp_te   <- sort(sample.int(S, floor(S * TEST_FRAC)))
  sp_tr   <- setdiff(seq_len(S), sp_te)

  Btr <- BIN[tr_plot, , drop = FALSE]; Bte <- BIN[te_plot, , drop = FALSE]
  Ctr <- as(crossprod(Btr), "CsparseMatrix")
  test_ppmi <- as.matrix(ppmi_from_C(crossprod(Bte), length(te_plot)))

  embs <- list()
  embs$CA      <- emb_ca(Btr)
  embs$PMISVD  <- emb_pmisvd(Ctr, length(tr_plot))
  embs$GloVe   <- emb_glove(Ctr)
  embs$AbundPMI <- tryCatch({                            # exploratory, non-blocking
    Wtr <- sqrt(COV[tr_plot, , drop = FALSE])
    Cw  <- as(crossprod(Wtr), "CsparseMatrix")
    emb_pmisvd(Cw, length(tr_plot))
  }, error = function(e) { cat("[warn] AbundPMI failed:", conditionMessage(e), "\n"); NULL })

  for (m in names(embs)) {
    emb <- embs[[m]]; if (is.null(emb)) next
    co <- cooc_spearman(emb, test_ppmi)
    eu <- eunis_metric(emb, tr_plot, te_plot)
    trr <- trait_metric(emb, sp_tr, sp_te)
    res[[length(res) + 1]] <- data.table(seed = seed, method = m,
                                          cooc = co, eunis_f1 = eu$f1, eunis_acc = eu$acc,
                                          eunis_base = eu$base, trait_r2 = trr)
    cat(sprintf("[seed %d] %-9s cooc=%.3f  eunisF1=%.3f (acc=%.3f base=%.3f)  traitR2=%.3f\n",
                seed, m, co, eu$f1, eu$acc, eu$base, trr))
  }
}
R <- rbindlist(res)

## ---- summarise + verdict ---------------------------------------------------
agg <- R[, .(cooc_m = mean(cooc), cooc_sd = sd(cooc),
             f1_m = mean(eunis_f1), f1_sd = sd(eunis_f1),
             r2_m = mean(trait_r2), r2_sd = sd(trait_r2)), by = method]
cat("\n================ SUMMARY (mean over seeds) ================\n")
print(agg[order(-cooc_m)], digits = 3)
fwrite(R, "C:/GillesC/documents/dev/specvec/bakeoff/bakeoff_10k_results.csv")

ca <- agg[method == "CA"]
cmp <- agg[method != "CA"]
beats2sd <- function(m, sd_self, ref_m, ref_sd) isTRUE((m - ref_m) > 2 * sqrt((sd_self^2 + ref_sd^2) / 2))
cat("\n================ VERDICT vs CA ================\n")
for (i in seq_len(nrow(cmp))) {
  row <- cmp[i]
  b_cooc <- beats2sd(row$cooc_m, row$cooc_sd, ca$cooc_m, ca$cooc_sd)
  b_f1   <- beats2sd(row$f1_m,   row$f1_sd,   ca$f1_m,   ca$f1_sd)
  b_r2   <- beats2sd(row$r2_m,   row$r2_sd,   ca$r2_m,   ca$r2_sd)
  nbeat  <- sum(b_cooc, b_f1, b_r2)
  verdict <- if (b_cooc || nbeat >= 2) "BEATS CA" else "no advantage over CA"
  cat(sprintf("%-9s  cooc%+.3f[%s]  f1%+.3f[%s]  r2%+.3f[%s]  -> %s\n",
              row$method, row$cooc_m - ca$cooc_m, ifelse(b_cooc,"Y","n"),
              row$f1_m - ca$f1_m, ifelse(b_f1,"Y","n"),
              row$r2_m - ca$r2_m, ifelse(b_r2,"Y","n"), verdict))
}
cat("\n[done]\n")
