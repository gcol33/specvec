## Cover-scale registry: map an abundance column to cover-proportion in [0,1].
## Ordinal vegetation scales are first-class rather than a pre-processing chore
## the caller must get right. Adding a scale (Domin, Londo, van der Maarel) is
## an O(1) registration, not a new function.

.specvec_scales <- new.env(parent = emptyenv())

#' Register a cover scale
#'
#' @param name Scale name.
#' @param fn Function `(x, mapping)` returning cover-proportion in `[0,1]`,
#'   with `NA` for missing or unrecognized values (imputed downstream).
#' @return Invisibly `NULL`; called for the side effect of registration.
#' @export
register_cover_scale <- function(name, fn) {
  stopifnot(is.function(fn))
  assign(name, fn, envir = .specvec_scales)
  invisible(NULL)
}

.get_cover_scale <- function(name) {
  if (!exists(name, envir = .specvec_scales, inherits = FALSE))
    stop(sprintf("unknown cover scale '%s'; available: %s", name,
                 paste(sort(ls(.specvec_scales)), collapse = ", ")), call. = FALSE)
  get(name, envir = .specvec_scales, inherits = FALSE)
}

#' Convert an abundance column to cover-proportion
#'
#' @param x Abundance values: numeric percent, numeric proportion, or ordinal
#'   cover-abundance codes.
#' @param scale One of `"percent"`, `"proportion"`, `"braun_blanquet"`, or any
#'   registered scale.
#' @param mapping Optional named numeric vector overriding the default lookup
#'   (used by ordinal scales such as Braun-Blanquet).
#' @return Numeric cover-proportion in `[0,1]`, `NA` where missing or
#'   unrecognized.
#' @export
#' @examples
#' cover_from_scale(c(0, 50, 100), scale = "percent")
#' cover_from_scale(c("r", "+", "2", "5"), scale = "braun_blanquet")
cover_from_scale <- function(x, scale = c("percent", "proportion", "braun_blanquet"),
                             mapping = NULL) {
  if (length(scale) > 1L) scale <- match.arg(scale)
  .get_cover_scale(scale)(x, mapping = mapping)
}

.scale_percent <- function(x, mapping = NULL) {
  v <- suppressWarnings(as.numeric(x))
  pmin(pmax(v, 0), 100) / 100
}

.scale_proportion <- function(x, mapping = NULL) {
  v <- suppressWarnings(as.numeric(x))
  pmin(pmax(v, 0), 1)
}

## Conventional Braun-Blanquet cover-class midpoints (percent). Conventions
## genuinely vary by author and the split classes (2a/2b/2m) are not universal,
## so this default is documented and overridable via `mapping`; it is not
## asserted as the only standard.
.braun_blanquet_default <- c(
  "r"  = 0.1,  "+"  = 0.5,
  "1"  = 2.5,  "2"  = 15,
  "2m" = 5,    "2a" = 10,   "2b" = 20,
  "3"  = 37.5, "4"  = 62.5, "5"  = 87.5
)

.scale_braun_blanquet <- function(x, mapping = NULL) {
  map <- if (is.null(mapping)) .braun_blanquet_default else mapping
  names(map) <- tolower(trimws(names(map)))
  key <- tolower(trimws(as.character(x)))
  out <- unname(map[key])
  unrec <- !is.na(key) & key != "" & is.na(out)
  if (any(unrec)) {
    eg <- paste(utils::head(unique(key[unrec]), 3L), collapse = ", ")
    warning(sprintf("braun_blanquet: %d unrecognized code(s) set to NA (e.g. %s)",
                    sum(unrec), eg), call. = FALSE)
  }
  pmin(pmax(out, 0), 100) / 100
}
