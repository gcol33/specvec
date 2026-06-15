# SpecVec bake-off -- presence/absence vs abundance, neutral co-occurrence scoring
# ================================================================================
# Usage: Rscript bakeoff.R [10000|50000|200000]   (default 10000)
#
# QUESTION (fixed): does the abundance-weighted PMI advantage over CA survive
#   NEUTRAL evaluation (not scored on PMI's own objective) and larger scale?
#
# METRICS (all methods produce 64-dim species vectors, scored on identical pairs):
#   cooc_ppmi : Spearman(gram, held-out PPMI)  -- HOME-FIELD for PMI, reference only
#   cooc_raw  : Spearman(gram, held-out raw co-occurrence counts)   -- NEUTRAL
#   link_auc  : AUC ranking co-occurring vs never-co-occurring pairs -- NEUTRAL
#   eunis_f1  : EUNIS-lvl2 macro-F1 (community embedding -> kNN)      -- NEUTRAL downstream
#   trait_r2  : mean trait R^2 (kept for completeness; not under investigation)
#
# WIN CRITERION (fixed before looking): a method "beats CA" iff it exceeds CA's
#   across-seed mean on BOTH neutral co-occurrence metrics (cooc_raw AND link_auc)
#   by more than 2 pooled SDs. If only cooc_ppmi separates them, the advantage was
#   objective alignment, not a better representation. Reported regardless.
# Fairness: identical train/test plot split, identical scored pair set, identical
#   pos/neg edge sets, and identical community-embedding definition (mean of present
#   species vectors) across every method within a seed. Diagonal zeroed everywhere.

suppressWarnings(suppressMessages({
  library(data.table); library(Matrix); library(RSpectra); library(text2vec); library(FNN)
}))
try(lgr::get_logger("text2vec")$set_threshold("warn"), silent = TRUE)
try(lgr::threshold("warn"), silent = TRUE)

## ---- args / constants ------------------------------------------------------
ARGS       <- commandArgs(trailingOnly = TRUE)
SAMPLE     <- if (length(ARGS) >= 1L) ARGS[1] else "10000"
DIM        <- 64L
MIN_OCC    <- if (length(ARGS) >= 2L) as.integer(ARGS[2]) else 10L   # second arg overrides
KNN_K      <- 15L
SEEDS      <- 1:3
TEST_FRAC  <- 0.20
GLOVE_ITER <- 20L
TRAITS     <- c("SLA","PlantHeight","SeedMass","LDMC")
PAIR_ALL_MAX <- 2.5e6      # score all upper-tri pairs below this, else sample
PAIR_SAMPLE  <- 2.0e6
AUC_CAP      <- 50000L     # max pos / neg edges for link AUC
EUNIS_REF_CAP<- 30000L     # max train plots used as kNN reference (scaling guard)
EUNIS_Q_CAP  <- 10000L     # max test plots scored for EUNIS (scaling guard)
DATA       <- "J:/Phd Local/Gilles_paper_resolve/data"
SP_FILE    <- file.path(DATA, sprintf("species_preprocessed_sample%s.csv", SAMPLE))
HD_FILE    <- file.path(DATA, sprintf("header_preprocessed_sample%s.csv", SAMPLE))
OUT_CSV    <- sprintf("C:/GillesC/documents/dev/specvec/bakeoff/bakeoff_%s_mo%d_results.csv", SAMPLE, MIN_OCC)

cat(sprintf("[cfg] SAMPLE=%s DIM=%d MIN_OCC=%d KNN_K=%d seeds=%s test=%.2f glove_iter=%d\n",
            SAMPLE, DIM, MIN_OCC, KNN_K, paste(SEEDS, collapse=","), TEST_FRAC, GLOVE_ITER))

## ---- load + filter ---------------------------------------------------------
sp <- fread(SP_FILE, showProgress = FALSE)
setnames(sp, names(sp), gsub("^'|'$", "", names(sp)))
hd <- fread(HD_FILE, select = c("PlotObservationID","Eunis_lvl2"), showProgress = FALSE)
stopifnot(all(c("PlotObservationID","WFO_TAXON","WFO_TAXON_RANK","Cover %") %in% names(sp)))
setnames(sp, c("PlotObservationID","WFO_TAXON","Cover %"), c("plot","species","cover"))
cat(sprintf("[load] species rows=%d  header rows=%d  distinct plots=%d  distinct species=%d\n",
            nrow(sp), nrow(hd), uniqueN(sp$plot), uniqueN(sp$species)))

n0 <- nrow(sp)
sp <- sp[WFO_TAXON_RANK == "species" & !is.na(species) & species != ""]
if ("EVA_TAXON_GROUP" %in% names(sp)) sp <- sp[EVA_TAXON_GROUP == "Vascular plant"]
freq <- sp[, .(f = uniqueN(plot)), by = species]
keep_sp <- freq[f >= MIN_OCC, species]
sp <- sp[species %in% keep_sp]
cat(sprintf("[filter] %d -> %d records; species kept %d/%d (min_occ>=%d)\n",
            n0, nrow(sp), length(keep_sp), nrow(freq), MIN_OCC))

## ---- global index + matrices ----------------------------------------------
plots   <- sort(unique(sp$plot)); species <- sort(unique(sp$species))
pi_ <- match(sp$plot, plots); si_ <- match(sp$species, species)
M <- length(plots); S <- length(species)
cat(sprintf("[index] plots=%d  species=%d\n", M, S))
stopifnot(M > 500, S >= 200, S > DIM + 5)

pa  <- unique(data.table(p = pi_, s = si_))
BIN <- sparseMatrix(i = pa$p, j = pa$s, x = 1, dims = c(M, S))
covdt <- sp[, .(cov = suppressWarnings(max(cover, na.rm = TRUE))), by = .(p = pi_, s = si_)]
covdt[!is.finite(cov), cov := NA_real_]
mc <- median(covdt$cov, na.rm = TRUE); if (!is.finite(mc)) mc <- 1
covdt[is.na(cov), cov := mc]
COV <- sparseMatrix(i = covdt$p, j = covdt$s, x = pmin(pmax(covdt$cov,0),100)/100, dims = c(M, S))

eun_map <- hd[match(plots, PlotObservationID), Eunis_lvl2]
eun_map[eun_map %in% c("", "~", "NA")] <- NA
tr <- sp[, lapply(.SD, function(x){ v <- suppressWarnings(as.numeric(x)); m <- median(v[is.finite(v)]); if (is.finite(m)) m else NA_real_ }),
         by = species, .SDcols = TRAITS]
setkey(tr, species); TRAITM <- as.matrix(tr[species, ..TRAITS])
cat(sprintf("[aux] EUNIS plots=%d/%d (classes=%d); full-trait species=%d\n",
            sum(!is.na(eun_map)), M, length(unique(na.omit(eun_map))),
            sum(rowSums(is.finite(TRAITM)) == length(TRAITS))))

## ---- embedding methods -----------------------------------------------------
ppmi_sparse <- function(C, nplots) {                       # cleaned sparse PPMI
  C <- as(C, "TsparseMatrix"); i <- C@i + 1L; j <- C@j + 1L; x <- C@x
  f <- Matrix::diag(as(C, "CsparseMatrix"))
  off <- i != j & x > 0 & f[i] > 0 & f[j] > 0
  i <- i[off]; j <- j[off]
  pp <- pmax(log(x[off] * nplots / (f[i] * f[j])), 0)
  keep <- pp > 0
  sparseMatrix(i = i[keep], j = j[keep], x = pp[keep], dims = dim(C))
}
emb_pmisvd <- function(C, nplots) {
  PP <- ppmi_sparse(C, nplots); if (length(PP@x) == 0) return(matrix(0, S, DIM))
  e <- tryCatch(RSpectra::eigs_sym(PP, k = DIM, which = "LA"), error = function(err) NULL)
  if (is.null(e) || all(!is.finite(e$values)) || max(e$values, na.rm = TRUE) <= 0) {
    PPd <- as.matrix(PP); PPd <- (PPd + t(PPd))/2; e <- eigen(PPd, symmetric = TRUE)
  }
  k <- min(DIM, sum(e$values > 1e-8)); if (k < 1) return(matrix(0, S, DIM))
  emb <- e$vectors[, seq_len(k), drop = FALSE] %*% diag(sqrt(e$values[seq_len(k)]), k, k)
  if (k < DIM) emb <- cbind(emb, matrix(0, S, DIM - k)); emb
}
emb_ca_dense <- function(Btr) {
  P <- as.matrix(Btr); tot <- sum(P); P <- P/tot
  r <- rowSums(P); c <- colSums(P); cs <- ifelse(c > 0, c, NA)
  E <- outer(r, c); Sr <- (P - E)/sqrt(E); Sr[!is.finite(Sr)] <- 0
  sv <- RSpectra::svds(Sr, k = DIM, nu = 0, nv = DIM)
  emb <- sweep(sv$v, 1, sqrt(cs), "/"); emb <- sweep(emb, 2, sv$d, "*"); emb[!is.finite(emb)] <- 0; emb
}
emb_ca <- function(Btr) {                                  # implicit matvec, no densify
  N <- Btr; tot <- sum(N)
  r <- as.numeric(Matrix::rowSums(N))/tot; c <- as.numeric(Matrix::colSums(N))/tot
  rs <- ifelse(r > 0, r, NA); cs <- ifelse(c > 0, c, NA); Mn <- nrow(N); Sn <- ncol(N)
  Af  <- function(x, args){ w <- x/sqrt(cs); w[!is.finite(w)] <- 0; pw <- as.numeric(N %*% w)/tot
                            o <- (pw - r*sum(c*w))/sqrt(rs); o[!is.finite(o)] <- 0; o }
  Atf <- function(x, args){ z <- x/sqrt(rs); z[!is.finite(z)] <- 0; ptz <- as.numeric(Matrix::crossprod(N, z))/tot
                            o <- (ptz - c*sum(r*z))/sqrt(cs); o[!is.finite(o)] <- 0; o }
  sv <- RSpectra::svds(Af, k = DIM, nu = 0, nv = DIM, dim = c(Mn, Sn), Atrans = Atf, args = NULL)
  emb <- sweep(sv$v, 1, sqrt(cs), "/"); emb <- sweep(emb, 2, sv$d, "*"); emb[!is.finite(emb)] <- 0; emb
}
emb_glove <- function(C) {
  tcm <- as(C, "TsparseMatrix"); Matrix::diag(tcm) <- 0; dimnames(tcm) <- list(species, species)
  gv <- GlobalVectors$new(rank = DIM, x_max = 10)
  wv <- suppressMessages(gv$fit_transform(tcm, n_iter = GLOVE_ITER, convergence_tol = 0.01, progressbar = FALSE))
  (wv + t(gv$components))[species, , drop = FALSE]
}

## ---- metrics ---------------------------------------------------------------
community_emb <- function(B, emb) as.matrix(B %*% emb) / pmax(Matrix::rowSums(B), 1)
scale_tt <- function(tr, te) { m <- colMeans(tr); s <- apply(tr, 2, sd); s[!is.finite(s) | s == 0] <- 1
  list(tr = sweep(sweep(tr, 2, m), 2, s, "/"), te = sweep(sweep(te, 2, m), 2, s, "/")) }
macroF1 <- function(pred, truth) { pred <- as.character(pred); truth <- as.character(truth)
  mean(vapply(sort(unique(truth)), function(k){ tp <- sum(pred==k & truth==k); fp <- sum(pred==k & truth!=k); fn <- sum(pred!=k & truth==k)
    p <- if (tp+fp==0) 0 else tp/(tp+fp); r <- if (tp+fn==0) 0 else tp/(tp+fn); if (p+r==0) 0 else 2*p*r/(p+r) }, numeric(1))) }

co_metrics <- function(emb, ev) {
  g <- rowSums(emb[ev$ip, , drop=FALSE] * emb[ev$jp, , drop=FALSE])
  ppmi <- suppressWarnings(cor(g, ev$test_ppmi, method = "spearman"))
  raw  <- suppressWarnings(cor(g, ev$test_raw,  method = "spearman"))
  if (length(ev$pos) >= 20 && length(ev$neg) >= 20) {
    sc <- c(g[ev$pos], g[ev$neg]); rk <- rank(sc)
    np <- as.numeric(length(ev$pos)); nn <- as.numeric(length(ev$neg))   # avoid 32-bit overflow
    auc <- (sum(rk[seq_len(np)]) - np*(np+1)/2)/(np*nn)
  } else auc <- NA_real_
  c(ppmi = ppmi, raw = raw, auc = auc)
}
eunis_metric <- function(emb, tr_plot, te_plot) {
  trp <- tr_plot[!is.na(eun_map[tr_plot])]; tep <- te_plot[!is.na(eun_map[te_plot])]
  if (length(trp) > EUNIS_REF_CAP) trp <- sample(trp, EUNIS_REF_CAP)
  if (length(tep) > EUNIS_Q_CAP)   tep <- sample(tep, EUNIS_Q_CAP)
  ytr <- eun_map[trp]; tep <- tep[eun_map[tep] %in% unique(ytr)]; yte <- eun_map[tep]
  sc <- scale_tt(community_emb(BIN[trp, , drop=FALSE], emb), community_emb(BIN[tep, , drop=FALSE], emb))
  idx <- FNN::knnx.index(sc$tr, sc$te, k = KNN_K)
  pred <- apply(idx, 1, function(ix){ tb <- table(ytr[ix]); names(tb)[which.max(tb)] })
  maj <- names(sort(table(ytr), decreasing = TRUE))[1]
  c(f1 = macroF1(pred, yte), acc = mean(as.character(pred) == yte), base = mean(yte == maj))
}
trait_metric <- function(emb, sp_tr, sp_te) {
  r2 <- c()
  for (t in TRAITS) { y <- TRAITM[, t]
    itr <- sp_tr[is.finite(y[sp_tr])]; ite <- sp_te[is.finite(y[sp_te])]
    if (length(itr) < KNN_K+5 || length(ite) < 5) next
    sc <- scale_tt(emb[itr,,drop=FALSE], emb[ite,,drop=FALSE])
    pred <- FNN::knn.reg(sc$tr, sc$te, y[itr], k = KNN_K)$pred
    sse <- sum((y[ite]-pred)^2); sst <- sum((y[ite]-mean(y[itr]))^2)
    r2 <- c(r2, if (sst > 0) 1 - sse/sst else NA_real_) }
  mean(r2, na.rm = TRUE)
}
make_pairs <- function() {
  full <- S*(S-1)/2
  if (full <= PAIR_ALL_MAX) { ij <- which(upper.tri(matrix(0L, S, S)), arr.ind = TRUE); return(list(ip = ij[,1], jp = ij[,2])) }
  n <- ceiling(PAIR_SAMPLE * 1.4); i <- sample.int(S, n, TRUE); j <- sample.int(S, n, TRUE)
  k <- i < j; i <- i[k]; j <- j[k]; d <- !duplicated(i * (S + 1) + j)
  i <- i[d]; j <- j[d]; m <- min(PAIR_SAMPLE, length(i)); list(ip = i[seq_len(m)], jp = j[seq_len(m)])
}

## ---- CA dense-vs-implicit equivalence guard (seed 1 train split) -----------
set.seed(1); .te1 <- sort(sample.int(M, floor(M*TEST_FRAC))); .tr1 <- setdiff(seq_len(M), .te1)
if (as.numeric(SAMPLE) <= 50000) {
  ca_d <- emb_ca_dense(BIN[.tr1, , drop=FALSE]); ca_i <- emb_ca(BIN[.tr1, , drop=FALSE])
  pr <- make_pairs()
  Cte1 <- as(crossprod(BIN[.te1,,drop=FALSE]), "CsparseMatrix")
  fte <- Matrix::diag(Cte1); tr_raw1 <- rep(0, length(pr$ip)); te_raw1 <- Cte1[cbind(pr$ip, pr$jp)]
  ppmi1 <- pmax(ifelse(te_raw1 > 0 & fte[pr$ip]*fte[pr$jp] > 0, log(te_raw1*length(.te1)/(fte[pr$ip]*fte[pr$jp])), 0), 0)
  ev1 <- list(ip = pr$ip, jp = pr$jp, test_ppmi = ppmi1, test_raw = te_raw1, pos = integer(0), neg = integer(0))
  cat(sprintf("[ca-check] dense cooc_raw=%.4f  implicit cooc_raw=%.4f  (should match)\n",
              co_metrics(ca_d, ev1)["raw"], co_metrics(ca_i, ev1)["raw"]))
}

## ---- run -------------------------------------------------------------------
res <- list()
for (seed in SEEDS) {
  set.seed(seed)
  te_plot <- sort(sample.int(M, floor(M*TEST_FRAC))); tr_plot <- setdiff(seq_len(M), te_plot)
  sp_te <- sort(sample.int(S, floor(S*TEST_FRAC)));   sp_tr  <- setdiff(seq_len(S), sp_te)
  Btr <- BIN[tr_plot,,drop=FALSE]; Bte <- BIN[te_plot,,drop=FALSE]
  Ctr <- as(crossprod(Btr), "CsparseMatrix"); Cte <- as(crossprod(Bte), "CsparseMatrix")

  pr <- make_pairs(); fte <- Matrix::diag(Cte)
  te_raw <- Cte[cbind(pr$ip, pr$jp)]; tr_raw <- Ctr[cbind(pr$ip, pr$jp)]
  den <- fte[pr$ip]*fte[pr$jp]
  te_ppmi <- pmax(ifelse(te_raw > 0 & den > 0, log(te_raw*length(te_plot)/den), 0), 0)
  pos <- which(te_raw > 0); neg <- which(te_raw == 0 & tr_raw == 0)
  if (length(pos) > AUC_CAP) pos <- sample(pos, AUC_CAP)
  if (length(neg) > AUC_CAP) neg <- sample(neg, AUC_CAP)
  ev <- list(ip = pr$ip, jp = pr$jp, test_ppmi = te_ppmi, test_raw = te_raw, pos = pos, neg = neg)
  cat(sprintf("[seed %d] scored pairs=%d  pos=%d neg=%d\n", seed, length(pr$ip), length(pos), length(neg)))

  embs <- list(CA = emb_ca(Btr), PMISVD = emb_pmisvd(Ctr, length(tr_plot)), GloVe = emb_glove(Ctr),
               AbundPMI = tryCatch(emb_pmisvd(as(crossprod(sqrt(COV[tr_plot,,drop=FALSE])), "CsparseMatrix"), length(tr_plot)),
                                   error = function(e){ cat("[warn] AbundPMI:", conditionMessage(e), "\n"); NULL }))
  for (m in names(embs)) { emb <- embs[[m]]; if (is.null(emb)) next
    co <- co_metrics(emb, ev); eu <- eunis_metric(emb, tr_plot, te_plot); trr <- trait_metric(emb, sp_tr, sp_te)
    res[[length(res)+1]] <- data.table(seed=seed, method=m, cooc_ppmi=co["ppmi"], cooc_raw=co["raw"],
                                        link_auc=co["auc"], eunis_f1=eu["f1"], eunis_acc=eu["acc"],
                                        eunis_base=eu["base"], trait_r2=trr)
    cat(sprintf("[seed %d] %-9s ppmi=%.3f raw=%.3f auc=%.3f | eunisF1=%.3f(acc=%.3f base=%.3f) traitR2=%.3f\n",
                seed, m, co["ppmi"], co["raw"], co["auc"], eu["f1"], eu["acc"], eu["base"], trr)) }
}
R <- rbindlist(res); fwrite(R, OUT_CSV)

## ---- summarise + verdict ---------------------------------------------------
agg <- R[, .(ppmi_m=mean(cooc_ppmi), raw_m=mean(cooc_raw), raw_sd=sd(cooc_raw),
             auc_m=mean(link_auc), auc_sd=sd(link_auc), f1_m=mean(eunis_f1), f1_sd=sd(eunis_f1),
             r2_m=mean(trait_r2)), by=method]
cat("\n================ SUMMARY (mean over seeds) ================\n"); print(agg[order(-auc_m)], digits=3)
ca <- agg[method=="CA"]; cmp <- agg[method!="CA"]
beats <- function(m, sds, rm, rsd) isTRUE((m-rm) > 2*sqrt((sds^2+rsd^2)/2))
cat("\n================ VERDICT vs CA (neutral metrics) ================\n")
for (i in seq_len(nrow(cmp))) { row <- cmp[i]
  br <- beats(row$raw_m, row$raw_sd, ca$raw_m, ca$raw_sd); ba <- beats(row$auc_m, row$auc_sd, ca$auc_m, ca$auc_sd)
  v <- if (br && ba) "BEATS CA (both neutral)" else if (br || ba) "partial (one neutral metric)" else "no neutral advantage"
  cat(sprintf("%-9s rawSpear%+.3f[%s] linkAUC%+.3f[%s]  (ppmi%+.3f home-field)  -> %s\n",
              row$method, row$raw_m-ca$raw_m, ifelse(br,"Y","n"), row$auc_m-ca$auc_m, ifelse(ba,"Y","n"),
              row$ppmi_m-ca$ppmi_m, v)) }
cat("\n[done]\n")
writeLines(c(sprintf("done %s", SAMPLE), capture.output(print(agg[order(-auc_m)], digits = 3))),
           sub("\\.csv$", ".done", OUT_CSV))
