# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`specvec` is an R package that turns a long `plot x species (x abundance)` table into abundance-aware species embeddings and first-class community embeddings, plus a benchmark harness that scores embedding methods on neutral ecological tasks. Current version 0.3.0 (all milestones 1-7 complete).

`IMPLEMENTATION_PLAN.md` is the design source of truth: method definitions (section 6), the registry architecture (section 3), the benchmark protocol (6.7), the roadmap (section 10), and the full validated evidence tables (section 13). Read it before changing any method math or the dispatch machinery; it records decisions already settled by benchmark (e.g. `abund_pmi` default, `min_occurrence = 5`, compositional A-vs-B resolved in favour of A).

## Paths

The repo is reachable at two paths for the same directory: `C:\Users\Gilles Colling\Documents\dev\specvec` and the space-free junction `C:/GillesC/Documents/dev/specvec`. All `dev_notes/*.R` scripts hardcode the junction path. Prefer the junction in any R/shell command.

## Development commands

Run R via Windows Rscript (resolve the latest: `RSCRIPT="$(ls -d '/c/Program Files/R/R-'*/bin/Rscript.exe | sort -V | tail -1)"`). Never run R with `Rscript.exe -e '...'` inline (segfaults on Windows with complex code); the workflow is driven by scripts in `dev_notes/`:

- **Document + load + test:** `Rscript dev_notes/run_tests.R` — roxygenise, `load_all`, run the full testthat suite, print pass/fail totals. The main loop while developing.
- **Run one test file:** edit `run_tests.R` to call `testthat::test_file(file.path(pkg, "tests/testthat/test-recovery.R"))`, or after `load_all` call `testthat::test_file(...)` directly. Test files: `test-ppmi`, `test-recovery`, `test-equivalence`, `test-dispatch`, `test-nearest`, `test-community`, `test-benchmark`, `test-temporal`, `test-integration`, `test-compositional` (shared sim helpers in `helper-sim.R`).
- **Full R CMD check:** `Rscript dev_notes/check.R` — `devtools::check(--as-cran-ish)`, prints errors/warnings/notes. Must stay 0/0/0 before any release.
- **Knit vignettes (quick error catch):** `Rscript dev_notes/knit_check.R` — installs the package and renders all vignettes to a temp dir.
- **Build pkgdown site:** `Rscript dev_notes/build_site.R` (wraps `~/.R/build_pkgdown.R`).
- **Bake-off (development data, not shipped):** `cd bakeoff && Rscript bakeoff.R <sample_size> <min_occurrence>` (e.g. `Rscript bakeoff.R 10000 5`); standalone reference implementation, writes `bakeoff_<sample>_mo<min_occ>_results.csv`. Other `dev_notes/run*.R` and `run_AvsB_10k.R` drive the package-native `compare_embeddings` over the EVA samples.

## Architecture

### Registry-based method dispatch (the core design)

A "method" is a registered `(weighting, factorization)` pair plus a capability descriptor. New methods are added by registration, never by copying a function (O(1) feature-add). This is the invariant to preserve when extending the package: do not write a `method_variant()` function; register a weighting or factorizer and add a method preset.

Three registries live in `R/registry.R` (weighting, factorization, method). `R/cover-scale.R` holds a fourth (cover scales). All four are populated in `R/zzz.R`'s `.onLoad`. To add a method end to end: write the weighting builder and/or factorizer (`R/weighting.R`, `R/factorization.R`), then add `register_*` calls in `zzz.R`.

The pipeline:

```
plot x species (+ abundance, + time)
   -> weighting   (counts | ppmi | abundance_pmi | chi_square | clr)   builds the operator
   -> factorization (eigen | svd | glove)                              consumes the operator
   -> species embedding -> community embedding (pooled species vectors)
```

### Operator kinds and the capability contract

A weighting returns an operator with a `kind` field, and a factorization declares the kind it `accepts`. `.fit_embedding()` (in `R/species-embedding.R`) checks they match and errors otherwise. The three kinds:

- `"sym"` — sparse symmetric species x species matrix (`$M`), consumed by `eigen`. Used by `ppmi`, `abundance_pmi`, `clr`.
- `"counts"` — sparse species x species co-occurrence counts (`$M`), consumed by `glove`.
- `"implicit"` — matvec closures (`$Af`/`$Atf`) plus a `$readout`, no densify, consumed by `svd`. Used by `ca` (the chi-square plot x species operand).

The capability descriptor on a method (`input`, `native_output`, `supports_abundance`, `default_dim`) drives routing and validation. `input` is `"species_species"` for PMI-family methods (factorize an association, read out species vectors, pool communities from them) or `"plot_species"` for CA (factorizes the contingency natively). The community readout is uniform across methods (pooled species vectors) so plot vectors stay comparable.

### Shared engine, no copy-paste

`.fit_embedding(data, ks, w_name, f_name, dim, n_plots, ...)` is the single weighting-then-factorize-then-sign-orient path. Both `species_embedding()` (global rare-species filter via `.kept_species`) and the benchmark (kept species fixed globally, data restricted to training plots) call it, so the two never duplicate the plumbing. When touching weighting/factorization logic, change it here, not per caller.

`.ppmi_sparse()` in `R/weighting.R` is shared by plain PMI and AbundPMI: AbundPMI is exactly PMI with the presence matrix replaced by `sqrt(COV)` (the operand becomes `A[a,b] = sum_p sqrt(cov[p,a]*cov[p,b])`, the geometric mean of covers). This is why abundance is a clean A/B and not a fork.

### Temporal and integration layers built on one primitive

`R/temporal.R` adds `time =` windowing to the core verbs (via `.time_rows`) and the fixed-frame projection `species_trajectory()`: fit one frame embedding, then place each focal species per window at the cover-weighted centroid of its co-occurrents in that fixed frame (no per-window refit, so no cross-time alignment). `community_trajectory()` is novelty over time. `R/integration.R`'s flagship `integration_trajectory()` composes the same fixed-frame primitive with `.pool_rows()` to read a neophyte's distance to the native-community centroid across windows. Build new diachronic features on the fixed-frame primitive, not by refitting per window.

### Return objects

S3 classes, each self-describing (ids + method metadata + fitted matrices + the preprocessing that produced them): `specvec_data`, `specvec_embedding`, `specvec_community`, `specvec_benchmark`, `specvec_trajectory`, `specvec_community_trajectory`, `specvec_integration`. Print/`as.data.frame` methods live in `R/print.R` (and trajectory/integration `as.data.frame` for tidying). Embeddings are sign-oriented (each column flipped so its largest-magnitude entry is positive) and ids are stored sorted, so refits are byte-identical.

## Conventions specific to this package

- **External calls are fully namespace-qualified** (`Matrix::`, `RSpectra::`, `data.table::`, `methods::`). The NAMESPACE stays minimal; the only `importFrom` is the `data.table` NSE tokens (`:=`, `.N`), with `utils::globalVariables` in `R/specvec-package.R` covering `cov`, `N`, `p`, `s`.
- **Tests are recovery-first.** Shape/dispatch tests prove plumbing only. The bar is parameter recovery against simulated truth (Procrustes/distance-correlation above tolerance across seeds), plus equivalence guards (dense CA == implicit-matvec CA; abundance vs presence pooling). When adding a fitter or method, add a recovery test, not just a shape test.
- **Reference data is not package data.** The EVA/ASAAS samples on `J:\Phd Local\...` are development inputs. The package ships no large `data/` and no `J:`-path dependency; `compare_embeddings` takes a user-supplied `specvec` object.
- **`text2vec` (GloVe) and `FNN` (kNN) are Suggests**, reached opportunistically with a clean error or brute-force fallback. Core install needs only `Matrix` + `RSpectra` + `data.table`. Keep it that way; write compositional/numeric math natively rather than pulling a stack.
- **Dependency direction:** specvec never depends on RESOLVE; RESOLVE may later consume specvec embeddings.
