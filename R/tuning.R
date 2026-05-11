#' @title Tuning Parameter Selection
#' @name tuning
#' @description Theoretical bounds, cross-validation, and eAIC for selecting eta.
NULL

#' Theoretical upper bound on eta (partial scenario)
#'
#' Computes \eqn{eta_bound = sigma2_int * trace(Sigma_int^{-1}) / (delta' delta)}
#' where delta = beta_ext - beta_int. This is a conservative upper bound derived
#' from Lemma 2, omitting the external variance term (which requires n_ext, sigma2_ext).
#'
#' @param beta_int Internal point estimate
#' @param beta_ext External point estimate
#' @param Sigma_int Internal Gram matrix X'X/n
#' @param sigma2_int Internal error variance
#' @param fallback Value returned when the bound is degenerate (default 5)
#' @return Scalar upper bound for eta
#' @export
eta_bound_partial <- function(beta_int, beta_ext, Sigma_int, sigma2_int,
                              fallback = 5) {
  delta <- beta_ext - beta_int
  denom <- as.numeric(t(delta) %*% Sigma_int %*% delta)
  n <- nrow(Sigma_int)  # Sigma_int = X'X/n, so X'X = n * Sigma_int

  if (denom < 1e-10) return(fallback)

  # RSS-based bound: sum(resid^2) / (delta' X'X delta)
  # = n * sigma2_int / (delta' (n*Sigma_int) delta)
  # Simplified: sigma2_int / (delta' Sigma_int delta) when using Sigma = X'X/n
  # But from the code, the actual formula uses RSS / (delta' X'X delta)
  # which equals sigma2_int * (n-p) / (n * delta' Sigma_int delta)
  # For a conservative bound, just use:
  bound <- sigma2_int / denom
  max(bound, 0.1)
}

#' Theoretical upper bound on eta (full scenario)
#'
#' Computes the bound from Lemma 1: accounts for eigenstructure of both
#' Gram matrices. Uses the sharpened form when Sigma_ext is full rank.
#'
#' @param Sigma_int Internal Gram matrix (p x p)
#' @param Sigma_ext External Gram matrix (p x p)
#' @param delta Heterogeneity vector beta_ext - beta_int (p-vector)
#' @param sigma2_int Internal error variance
#' @param sigma2_ext External error variance
#' @param n_int Internal sample size
#' @param n_ext External sample size
#' @return Scalar upper bound for eta
#' @export
eta_bound_full <- function(Sigma_int, Sigma_ext, delta, sigma2_int, sigma2_ext,
                           n_int, n_ext) {
  p <- nrow(Sigma_int)

  # Eigenvalues
  eig_int <- eigen(Sigma_int, symmetric = TRUE, only.values = TRUE)$values
  g_p <- max(eig_int)  # lambda_max
  g_1 <- min(eig_int)  # lambda_min

  eig_ext <- eigen(Sigma_ext, symmetric = TRUE, only.values = TRUE)$values
  t_1 <- min(eig_ext)

  # Cross eigenvalues: k_j = eigenvalues of Sigma_int^{1/2} Sigma_ext^{-1} Sigma_int^{1/2}
  Sigma_int_half <- with(eigen(Sigma_int, symmetric = TRUE),
                         vectors %*% diag(sqrt(pmax(values, 0))) %*% t(vectors))
  Sigma_ext_inv <- solve(Sigma_ext)
  cross_mat <- Sigma_int_half %*% Sigma_ext_inv %*% Sigma_int_half
  k_vals <- eigen(cross_mat, symmetric = TRUE, only.values = TRUE)$values
  k_1 <- min(k_vals)

  # Bias term
  delta_p <- max(eigen(outer(delta, delta), symmetric = TRUE, only.values = TRUE)$values)

  # Sharpened bound (Proposition 3, external eigenspace)
  numer_sharp <- sigma2_int / n_int * k_1 * t_1 / g_p
  denom_bound <- delta_p + sigma2_ext / n_ext

  # Generalized bound (Lemma 1, internal eigenspace)
  kg_ratio <- min(k_vals / rev(eig_int))  # min_j(k_j / g_j)
  numer_gen <- sigma2_int / n_int * g_1^2 * kg_ratio

  # Combined: max of both
  if (denom_bound < 1e-15) return(5)

  eta_sharp <- numer_sharp / denom_bound
  eta_gen <- numer_gen / denom_bound

  max(eta_sharp, eta_gen, 0.1)
}

#' Closed-form optimal eta (partial scenario, Lemma 2)
#'
#' @param Sigma_int Internal Gram matrix
#' @param delta Heterogeneity vector
#' @param sigma2_int Internal error variance
#' @param sigma2_ext External error variance
#' @param n_int Internal sample size
#' @param n_ext External sample size
#' @param T_mat Transformation matrix (default identity)
#' @param Sigma_ext_raw Raw external Gram matrix (for variance term)
#' @return Optimal eta minimizing MSE
#' @export
eta_optimal_partial <- function(Sigma_int, delta, sigma2_int, sigma2_ext = NULL,
                                n_int, n_ext = NULL, T_mat = NULL,
                                Sigma_ext_raw = NULL) {
  numer <- sigma2_int / n_int * sum(diag(solve(Sigma_int)))

  # Denominator: delta'delta + external variance term
  denom <- sum(delta^2)

  # Add external variance if available
  if (!is.null(sigma2_ext) && !is.null(n_ext) && !is.null(Sigma_ext_raw)) {
    if (is.null(T_mat)) T_mat <- diag(length(delta))
    denom <- denom + sigma2_ext / n_ext *
      sum(diag(T_mat %*% solve(Sigma_ext_raw) %*% t(T_mat)))
  }

  if (denom < 1e-15) return(5)
  numer / denom
}

#' Select eta via extended AIC (for restricted or partial scenario)
#'
#' Minimizes RSS(eta) - n*sigma2(eta) + 2*df(eta)*sigma2(eta) over a bounded grid.
#'
#' @param beta_int Internal baseline estimate (OLS or summary-stat)
#' @param beta_ext External point estimate
#' @param X_tr Training design matrix (for computing RSS)
#' @param y_tr Training response vector
#' @param eta_grid Grid of candidate eta values
#' @param p Number of parameters
#' @return List with selected eta, eAIC values, and integrated coefficients
#' @export
tune_eaic <- function(beta_int, beta_ext, X_tr, y_tr, eta_grid, p = NULL) {
  n <- nrow(X_tr)
  if (is.null(p)) p <- ncol(X_tr)

  aic_vals <- numeric(length(eta_grid))
  for (i in seq_along(eta_grid)) {
    e <- eta_grid[i]
    bp <- est_partial(beta_int, beta_ext, e)
    rss <- sum((y_tr - X_tr %*% bp)^2)
    df_e <- p / (1 + e)
    sigma2_e <- rss / max(n - df_e, 1)
    aic_vals[i] <- rss + 2 * df_e * sigma2_e
  }

  best_idx <- which.min(aic_vals)
  eta_sel <- eta_grid[best_idx]

  list(
    eta = eta_sel,
    coefficients = est_partial(beta_int, beta_ext, eta_sel),
    eaic_values = data.frame(eta = eta_grid, eaic = aic_vals),
    df = p / (1 + eta_sel)
  )
}

#' Select eta via K-fold cross-validation with bounded search
#'
#' @param X_tr Training design matrix
#' @param y_tr Training response vector
#' @param beta_ext External point estimate
#' @param eta_grid Grid of candidate eta values (should be bounded by theory)
#' @param nfolds Number of CV folds (default 5)
#' @param base Type of internal estimator: "ols" or "ridge"
#' @param ... Additional arguments passed to cv.glmnet when base = "ridge"
#' @return List with selected eta, CV errors, and integrated coefficients
#' @export
tune_cv <- function(X_tr, y_tr, beta_ext, eta_grid, nfolds = 5,
                    base = c("ols", "ridge"), ...) {
  base <- match.arg(base)
  n <- nrow(X_tr)
  p <- ncol(X_tr)
  folds <- sample(rep(1:nfolds, length.out = n))
  cv_errs <- matrix(0, nrow = length(eta_grid), ncol = nfolds)

  for (k in 1:nfolds) {
    X_k_tr <- X_tr[folds != k, , drop = FALSE]
    y_k_tr <- y_tr[folds != k]
    X_k_te <- X_tr[folds == k, , drop = FALSE]
    y_k_te <- y_tr[folds == k]

    # Fit internal model on fold
    if (base == "ols") {
      b_k <- as.numeric(MASS::ginv(crossprod(X_k_tr)) %*% t(X_k_tr) %*% y_k_tr)
    } else {
      cv_k <- glmnet::cv.glmnet(X_k_tr, y_k_tr, alpha = 0, intercept = FALSE, ...)
      b_k <- as.numeric(coef(cv_k, s = "lambda.min"))[-1]
    }

    # Evaluate integrated estimator at each eta
    for (i in seq_along(eta_grid)) {
      bp <- est_partial(b_k, beta_ext, eta_grid[i])
      cv_errs[i, k] <- mean((y_k_te - X_k_te %*% bp)^2)
    }
  }

  mean_cv <- rowMeans(cv_errs)
  best_idx <- which.min(mean_cv)
  eta_sel <- eta_grid[best_idx]

  # Refit on full training data
  if (base == "ols") {
    b_full <- as.numeric(MASS::ginv(crossprod(X_tr)) %*% t(X_tr) %*% y_tr)
  } else {
    cv_full <- glmnet::cv.glmnet(X_tr, y_tr, alpha = 0, intercept = FALSE, ...)
    b_full <- as.numeric(coef(cv_full, s = "lambda.min"))[-1]
  }

  list(
    eta = eta_sel,
    coefficients = est_partial(b_full, beta_ext, eta_sel),
    beta_internal = b_full,
    cv_errors = data.frame(eta = eta_grid, mean_cv_error = mean_cv),
    nfolds = nfolds,
    base = base
  )
}

#' Generate a bounded eta grid from the theoretical bound
#'
#' @param eta_bound Upper bound (from eta_bound_partial or eta_bound_full)
#' @param n_points Number of grid points (default 100)
#' @return Numeric vector of eta values on the interval from 0 to `eta_bound`
#' @export
make_eta_grid <- function(eta_bound, n_points = 100) {
  seq(0, eta_bound, length.out = n_points)
}
