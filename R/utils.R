#' @title Internal Utilities
#' @name utils
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @importFrom utils head
#' @useDynLib mainger, .registration = TRUE
NULL

# Null-coalescing operator. Returns `a` if non-NULL, otherwise `b`.
# Internal helper; not exported.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Suppress R CMD check NOTE about unbound globals from ggplot2 NSE in
# plot_advantage(), plot_diagnostics(), plot_eta_curve(). These are
# column names referenced inside aes() and are not package-level
# variables.
utils::globalVariables(c("eta", "advantage", "index", "k", "criterion"))
