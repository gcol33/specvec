#' @keywords internal
"_PACKAGE"

## Operators flow weighting -> factorization. External calls are fully
## namespace-qualified (Matrix::, RSpectra::, data.table::) so the NAMESPACE
## stays minimal; the data.table NSE tokens below are the exception.
#' @importFrom data.table := .N
NULL

utils::globalVariables(c("cov", "N", "p", "s"))
