#' @title Covariate Space Transformations
#' @name transform
#' @description Utilities for handling heterogeneous parameter spaces.
NULL

#' Zero-pad an external estimate to match internal dimension
#'
#' @param beta_ext External estimate (p2-vector)
#' @param p_int Internal dimension (p1 >= p2)
#' @param shared_idx Indices of shared covariates in the internal model
#'   (default: 1:length(beta_ext), i.e., first p2 covariates are shared)
#' @return Zero-padded p1-vector
#' @export
zero_pad <- function(beta_ext, p_int, shared_idx = NULL) {
  p_ext <- length(beta_ext)
  if (p_ext >= p_int) return(beta_ext[1:p_int])
  if (is.null(shared_idx)) shared_idx <- 1:p_ext

  out <- numeric(p_int)
  out[shared_idx] <- beta_ext
  out
}

#' Regression-based imputation for covariate space transformation
#'
#' Constructs the CLS transformation: regresses exclusive internal covariates
#' on shared covariates, then projects the external estimate into the full
#' internal space.
#'
#' @param X_int Internal design matrix (n x p1)
#' @param beta_ext External estimate (p2-vector)
#' @param shared_idx Indices of shared covariates (default: 1:length(beta_ext))
#' @return List with transformed estimate (p1-vector) and projection matrix T
#' @export
regression_impute <- function(X_int, beta_ext, shared_idx = NULL) {
  p_ext <- length(beta_ext)
  p_int <- ncol(X_int)
  if (is.null(shared_idx)) shared_idx <- 1:p_ext
  excl_idx <- setdiff(1:p_int, shared_idx)

  if (length(excl_idx) == 0) {
    return(list(beta = beta_ext, T_mat = diag(p_int)))
  }

  X_S <- X_int[, shared_idx, drop = FALSE]
  X_E <- X_int[, excl_idx, drop = FALSE]

  # B = (X_S'X_S)^{-1} X_S' X_E
  B_hat <- MASS::ginv(crossprod(X_S)) %*% crossprod(X_S, X_E)

  # C = [I_p2, B_hat]
  # T = C' (C C')^{-1}
  # theta_tilde = T * beta_ext
  theta_tilde <- numeric(p_int)
  theta_tilde[shared_idx] <- beta_ext
  theta_tilde[excl_idx] <- as.numeric(t(B_hat) %*% beta_ext)

  # Build T matrix for reference
  C_mat <- matrix(0, p_ext, p_int)
  C_mat[, shared_idx] <- diag(p_ext)
  C_mat[, excl_idx] <- B_hat
  T_mat <- t(C_mat) %*% MASS::ginv(C_mat %*% t(C_mat))

  list(beta = theta_tilde, T_mat = T_mat)
}
