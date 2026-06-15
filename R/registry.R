## Capability-declaring registries. A "method" is a registered
## (weighting, factorization) pair plus a capability descriptor; new methods
## are added by registration, never by copying a function (O(1) feature-add).

.specvec_registry <- new.env(parent = emptyenv())
.specvec_registry$weighting     <- new.env(parent = emptyenv())
.specvec_registry$factorization <- new.env(parent = emptyenv())
.specvec_registry$method        <- new.env(parent = emptyenv())

#' Register a weighting
#'
#' A weighting builds the operator a factorization consumes. `fn` takes
#' `(data, kept_species, n_plots, min_cooccurrence)` and returns an operator
#' list with a `kind` field (`"sym"`, `"counts"`, or `"implicit"`).
#'
#' @param name Weighting name.
#' @param fn Operator-building function.
#' @param input `"species_species"` or `"plot_species"`.
#' @param supports_abundance Logical; whether the weighting uses cover.
#' @return Invisibly `NULL`; called for the side effect of registration.
#' @export
register_weighting <- function(name, fn, input = c("species_species", "plot_species"),
                               supports_abundance = FALSE) {
  input <- match.arg(input)
  stopifnot(is.function(fn))
  assign(name, list(fn = fn, input = input, supports_abundance = supports_abundance),
         envir = .specvec_registry$weighting)
  invisible(NULL)
}

#' Register a factorization
#'
#' @param name Factorization name.
#' @param fn Factorization function `(operator, dim, ...)` returning a
#'   species x dim matrix with species row names.
#' @param kind `"matrix"` (svd/eigen; composes with a symmetric weighting) or
#'   `"objective"` (brings its own loss and subsumes its weighting).
#' @param accepts Operator kind this factorization consumes
#'   (`"sym"`, `"implicit"`, or `"counts"`).
#' @return Invisibly `NULL`; called for the side effect of registration.
#' @export
register_factorization <- function(name, fn, kind = c("matrix", "objective"),
                                   accepts = c("sym", "implicit", "counts")) {
  kind <- match.arg(kind)
  accepts <- match.arg(accepts)
  stopifnot(is.function(fn))
  assign(name, list(fn = fn, kind = kind, accepts = accepts),
         envir = .specvec_registry$factorization)
  invisible(NULL)
}

#' Register a method preset
#'
#' Maps a friendly method name to a `(weighting, factorization)` pair plus a
#' capability descriptor used for dispatch and validation.
#'
#' @param name Method name (e.g. `"abund_pmi"`).
#' @param weighting Registered weighting name.
#' @param factorization Registered factorization name.
#' @param input `"species_species"` or `"plot_species"`.
#' @param native_output `"species"` or `"species_and_community"`.
#' @param supports_abundance Logical.
#' @param default_dim Default embedding dimension.
#' @return Invisibly `NULL`; called for the side effect of registration.
#' @export
register_method <- function(name, weighting, factorization,
                            input = c("species_species", "plot_species"),
                            native_output = c("species", "species_and_community"),
                            supports_abundance = FALSE, default_dim = 64L) {
  input <- match.arg(input)
  native_output <- match.arg(native_output)
  assign(name, list(weighting = weighting, factorization = factorization,
                    input = input, native_output = native_output,
                    supports_abundance = supports_abundance,
                    default_dim = as.integer(default_dim)),
         envir = .specvec_registry$method)
  invisible(NULL)
}

.get_weighting <- function(name) {
  if (!exists(name, envir = .specvec_registry$weighting, inherits = FALSE))
    stop(sprintf("unknown weighting '%s'; available: %s", name,
                 paste(.list_weightings(), collapse = ", ")), call. = FALSE)
  get(name, envir = .specvec_registry$weighting, inherits = FALSE)
}
.get_factorization <- function(name) {
  if (!exists(name, envir = .specvec_registry$factorization, inherits = FALSE))
    stop(sprintf("unknown factorization '%s'; available: %s", name,
                 paste(.list_factorizations(), collapse = ", ")), call. = FALSE)
  get(name, envir = .specvec_registry$factorization, inherits = FALSE)
}
.get_method <- function(name) {
  if (!exists(name, envir = .specvec_registry$method, inherits = FALSE))
    stop(sprintf("unknown method '%s'; available: %s", name,
                 paste(.list_methods(), collapse = ", ")), call. = FALSE)
  get(name, envir = .specvec_registry$method, inherits = FALSE)
}

.list_weightings     <- function() sort(ls(.specvec_registry$weighting))
.list_factorizations <- function() sort(ls(.specvec_registry$factorization))
.list_methods        <- function() sort(ls(.specvec_registry$method))

#' List registered methods
#'
#' @return A character vector of registered method names.
#' @export
#' @examples
#' specvec_methods()
specvec_methods <- function() .list_methods()
