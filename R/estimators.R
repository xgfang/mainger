#' @title Core Estimators for Mahalanobis-Integrated Regression
#' @name estimators
#' @description Closed-form integrated estimators for each data-sharing scenario.
NULL

#' Full data-sharing estimator
#'
#' Computes the matrix-weighted average:
#' \eqn{(Sigma_int + eta * Sigma_ext)^{-1}(Sigma_int * beta_int + eta * Sigma_ext * beta_ext)}
#'
#' @param beta_int Internal point estimate (p-vector)
#' @param Sigma_int Internal Gram matrix (p x p)
#' @param beta_ext External point estimate (p-vector)
#' @param Sigma_ext External Gram matrix (p x p)
#' @param eta Tuning parameter (scalar >= 0)
#' @return Integrated coefficient vector (p-vector)
#' @export
est_full <- function(beta_int, Sigma_int, beta_ext, Sigma_ext, eta) {
  p <- length(beta_int)
  M <- Sigma_int + eta * Sigma_ext
  rhs <- Sigma_int %*% beta_int + eta * Sigma_ext %*% beta_ext
  as.numeric(solve(M, rhs))
}

#' Partial data-sharing estimator
#'
#' Computes the scalar-weighted average:
#' \eqn{beta_int / (1 + eta) + eta * beta_ext / (1 + eta)}
#'
#' @param beta_int Internal point estimate (p-vector)
#' @param beta_ext External point estimate (p-vector)
#' @param eta Tuning parameter (scalar >= 0)
#' @return Integrated coefficient vector (p-vector)
#' @export
est_partial <- function(beta_int, beta_ext, eta) {
  (beta_int + eta * beta_ext) / (1 + eta)
}

#' Restricted data-sharing estimator
#'
#' Computes the summary-statistic-based estimator with optional ridge penalty:
#' First obtains internal baseline via \eqn{(Sigma_ref + lambda * I)^{-1} r_int},
#' then shrinks toward external.
#'
#' @param r_int Marginal correlations X'Y/n (p-vector)
#' @param Sigma_ref Reference LD matrix (p x p)
#' @param beta_ext External point estimate (p-vector)
#' @param eta Tuning parameter (scalar >= 0)
#' @param lambda Ridge penalty (scalar >= 0, default 0)
#' @return Integrated coefficient vector (p-vector)
#' @export
est_restricted <- function(r_int, Sigma_ref, beta_ext, eta, lambda = 0) {
  p <- length(r_int)
  M <- Sigma_ref + lambda * diag(p)
  beta_check <- as.numeric(solve(M, r_int))
  est_partial(beta_check, beta_ext, eta)
}

#' High-dimensional ridge-integrated estimator
#'
#' Computes \eqn{(Sigma_int + lambda * I + eta * Sigma_ext)^{-1}(r_int + eta * Sigma_ext * beta_ext)}
#'
#' @param r_int Marginal correlations or X'Y/n (p-vector)
#' @param Sigma_int Internal Gram matrix (p x p)
#' @param beta_ext External point estimate (p-vector)
#' @param Sigma_ext External Gram/penalty matrix (p x p)
#' @param eta Tuning parameter (scalar >= 0)
#' @param lambda Ridge penalty (scalar > 0)
#' @return Integrated coefficient vector (p-vector)
#' @export
est_ridge <- function(r_int, Sigma_int, beta_ext, Sigma_ext, eta, lambda) {
  p <- length(r_int)
  M <- Sigma_int + lambda * diag(p) + eta * Sigma_ext
  rhs <- r_int + eta * Sigma_ext %*% beta_ext
  as.numeric(solve(M, rhs))
}

#' Competitors: Empirical Bayes (EB) estimator
#'
#' @param beta_int Internal OLS estimate
#' @param beta_ext External point estimate
#' @param V_int Internal covariance matrix of beta_int
#' @param p Number of parameters
#' @return EB-integrated coefficient vector
#' @export
est_eb <- function(beta_int, beta_ext, V_int, p) {
  delta <- beta_ext - beta_int
  V_reg <- V_int + 1e-6 * diag(p)
  dist_sq <- as.numeric(t(delta) %*% MASS::ginv(V_reg) %*% delta)
  w <- max(0, 1 - (p - 2) / dist_sq)
  w * beta_ext + (1 - w) * beta_int
}

#' Competitors: Positive-part James-Stein (JS+) estimator
#'
#' @param beta_int Internal OLS estimate
#' @param beta_ext External point estimate
#' @param X_int Internal design matrix
#' @param sigma2_int Internal error variance
#' @param p Number of parameters
#' @return JS+-integrated coefficient vector
#' @export
est_js_plus <- function(beta_int, beta_ext, X_int, sigma2_int, p) {
  diff_pred <- X_int %*% (beta_int - beta_ext)
  js_denom <- sum(diff_pred^2)
  w <- max(0, 1 - (p - 2) * sigma2_int / js_denom)
  w * beta_int + (1 - w) * beta_ext
}
