#' @title Mahalanobis-Integrated Regression
#'
#' @description Unified interface for integrating external information into an
#'   internal regression model. Automatically detects the data-sharing scenario
#'   from the provided inputs and selects the appropriate estimator and tuning
#'   strategy.
#'
#' @param X_int Internal design matrix (n x p). If provided, summary statistics
#'   are computed automatically.
#' @param Y_int Internal response vector (n x 1).
#' @param beta_int Internal point estimate (p-vector). Computed from X_int, Y_int
#'   if not provided.
#' @param Sigma_int Internal Gram matrix X'X/n (p x p).
#' @param r_int Marginal correlations X'Y/n (p-vector). For restricted scenario.
#' @param Sigma_ref Reference LD/covariance matrix (p x p). For restricted scenario.
#' @param beta_ext External point estimate (p-vector). Required.
#' @param Sigma_ext External Gram matrix (p x p). If provided, uses full scenario.
#' @param sigma2_int Internal error variance. Estimated if X_int, Y_int provided.
#' @param sigma2_ext External error variance (for full scenario bounds).
#' @param n_int Internal sample size.
#' @param n_ext External sample size.
#' @param tuning Tuning method: "auto" (default), "cv", "eaic", or "fixed".
#' @param eta Fixed eta value (when tuning = "fixed").
#' @param eta_grid Custom grid for eta search. Default: bounded by theory.
#' @param nfolds Number of CV folds (default 5).
#' @param base Internal estimator type: "ols" (default) or "ridge".
#' @param lambda Ridge penalty. NULL = auto via cv.glmnet.
#' @param diagnostics Logical. Compute spectral diagnostics? (default TRUE)
#'
#' @return An object of class "mainger" containing coefficients, tuning results,
#'   and diagnostics.
#'
#' @examples
#' \dontrun{
#' # Partial sharing with individual data
#' fit <- mainger(X_int = X, Y_int = Y, beta_ext = theta_hat)
#' summary(fit)
#' plot(fit)
#'
#' # Restricted sharing with marginal correlations
#' fit_r <- mainger(r_int = XtY_over_n, Sigma_ref = LD_matrix,
#'                  beta_ext = gwas_betas, n_int = 2000)
#' }
#'
#' @export
mainger <- function(X_int = NULL, Y_int = NULL,
                    beta_int = NULL, Sigma_int = NULL,
                    r_int = NULL, Sigma_ref = NULL,
                    beta_ext,
                    Sigma_ext = NULL,
                    sigma2_int = NULL, sigma2_ext = NULL,
                    n_int = NULL, n_ext = NULL,
                    tuning = c("auto", "cv", "eaic", "fixed"),
                    eta = NULL, eta_grid = NULL, nfolds = 5,
                    base = c("ols", "ridge"), lambda = NULL,
                    diagnostics = TRUE) {

  tuning <- match.arg(tuning)
  base <- match.arg(base)
  cl <- match.call()

  # ===========================================================
  # Step 1: Compute summaries from individual data if provided
  # ===========================================================
  has_individual <- !is.null(X_int) && !is.null(Y_int)

  if (has_individual) {
    n <- nrow(X_int)
    p <- ncol(X_int)
    if (is.null(n_int)) n_int <- n
    if (is.null(Sigma_int)) Sigma_int <- crossprod(X_int) / n
    if (is.null(r_int)) r_int <- as.numeric(crossprod(X_int, Y_int) / n)

    if (is.null(beta_int)) {
      if (base == "ols" && p < n) {
        beta_int <- as.numeric(MASS::ginv(crossprod(X_int)) %*%
                                 crossprod(X_int, Y_int))
      } else {
        cv_fit <- glmnet::cv.glmnet(X_int, Y_int, alpha = 0,
                                     intercept = FALSE)
        beta_int <- as.numeric(coef(cv_fit, s = "lambda.min"))[-1]
        if (is.null(lambda)) lambda <- cv_fit$lambda.min
      }
    }

    if (is.null(sigma2_int)) {
      resid <- Y_int - X_int %*% beta_int
      df_resid <- max(n - p, 1)
      sigma2_int <- sum(resid^2) / df_resid
    }
  }

  p <- length(beta_ext)

  # ===========================================================
  # Step 2: Auto-detect scenario
  # ===========================================================
  if (!is.null(Sigma_ext) && !is.null(Sigma_int)) {
    scenario <- "full"
  } else if (!is.null(beta_int) && !is.null(Sigma_int)) {
    scenario <- "partial"
  } else if (!is.null(r_int) && !is.null(Sigma_ref)) {
    scenario <- "restricted"
  } else {
    stop("Cannot determine data-sharing scenario from provided inputs.\n",
         "Provide one of:\n",
         "  Full:       beta_int + Sigma_int + Sigma_ext\n",
         "  Partial:    beta_int + Sigma_int (or X_int + Y_int)\n",
         "  Restricted: r_int + Sigma_ref")
  }

  # Align dimensions (zero-pad if beta_ext is shorter)
  if (length(beta_ext) < length(beta_int)) {
    p_int <- length(beta_int)
    p_ext <- length(beta_ext)
    beta_ext_padded <- numeric(p_int)
    beta_ext_padded[1:p_ext] <- beta_ext
    beta_ext <- beta_ext_padded
    p <- p_int
    message(sprintf("Zero-padded external estimate from %d to %d dimensions.",
                    p_ext, p_int))
  }

  # ===========================================================
  # Step 3: Compute theoretical bound
  # ===========================================================
  if (scenario == "full" && !is.null(sigma2_int) && !is.null(sigma2_ext) &&
      !is.null(n_int) && !is.null(n_ext)) {
    delta <- beta_ext - beta_int
    eta_ub <- eta_bound_full(Sigma_int, Sigma_ext, delta,
                              sigma2_int, sigma2_ext, n_int, n_ext)
  } else if (scenario %in% c("partial", "restricted")) {
    Sig <- if (scenario == "partial") Sigma_int else Sigma_ref
    bi <- if (scenario == "partial") beta_int else {
      as.numeric(solve(Sig + ifelse(is.null(lambda), 0, lambda) * diag(p),
                       r_int))
    }
    eta_ub <- eta_bound_partial(bi, beta_ext, Sig,
                                 ifelse(is.null(sigma2_int), 1, sigma2_int))
  } else {
    eta_ub <- 5
  }

  # ===========================================================
  # Step 4: Select eta
  # ===========================================================
  if (is.null(eta_grid)) {
    eta_grid <- make_eta_grid(eta_ub, n_points = 100)
  }

  if (tuning == "auto") {
    tuning <- if (has_individual) "cv" else "eaic"
  }

  tuning_result <- NULL

  if (tuning == "fixed") {
    if (is.null(eta)) stop("Must provide eta when tuning = 'fixed'.")
    eta_sel <- eta

  } else if (tuning == "cv") {
    if (!has_individual) stop("CV tuning requires individual-level data (X_int, Y_int).")
    tuning_result <- tune_cv(X_int, Y_int, beta_ext, eta_grid,
                              nfolds = nfolds, base = base)
    eta_sel <- tuning_result$eta
    beta_int <- tuning_result$beta_internal

  } else if (tuning == "eaic") {
    if (!has_individual && scenario == "restricted") {
      # Need to construct X_tr-like quantity from summaries for RSS
      # Use Sigma_ref and r_int to evaluate the criterion analytically
      bi <- as.numeric(solve(Sigma_ref + ifelse(is.null(lambda), 0, lambda) * diag(p),
                             r_int))
      # Analytical eAIC using summary statistics
      aic_vals <- numeric(length(eta_grid))
      for (i in seq_along(eta_grid)) {
        e <- eta_grid[i]
        bp <- est_partial(bi, beta_ext, e)
        rss <- as.numeric(n_int * (t(bp) %*% Sigma_ref %*% bp -
                                    2 * t(r_int) %*% bp +
                                    t(beta_ext) %*% Sigma_ref %*% beta_ext))
        df_e <- p / (1 + e)
        sigma2_e <- rss / max(n_int - df_e, 1)
        aic_vals[i] <- rss - n_int * sigma2_e + 2 * df_e * sigma2_e
      }
      best_idx <- which.min(aic_vals)
      eta_sel <- eta_grid[best_idx]
      beta_int <- bi
      tuning_result <- list(eta = eta_sel,
                            eaic_values = data.frame(eta = eta_grid, eaic = aic_vals))
    } else if (has_individual) {
      tuning_result <- tune_eaic(beta_int, beta_ext, X_int, Y_int, eta_grid, p)
      eta_sel <- tuning_result$eta
    } else {
      stop("eAIC requires either individual data or restricted-scenario inputs.")
    }
  }

  # ===========================================================
  # Step 5: Compute final integrated estimate
  # ===========================================================
  if (scenario == "full") {
    coefficients <- est_full(beta_int, Sigma_int, beta_ext, Sigma_ext, eta_sel)
  } else if (scenario == "partial") {
    coefficients <- est_partial(beta_int, beta_ext, eta_sel)
  } else {
    coefficients <- est_restricted(r_int, Sigma_ref, beta_ext, eta_sel,
                                    lambda = ifelse(is.null(lambda), 0, lambda))
  }

  # ===========================================================
  # Step 6: Compute diagnostics
  # ===========================================================
  diag_out <- NULL
  if (diagnostics && scenario %in% c("full", "partial") && !is.null(Sigma_int)) {
    diag_out <- list()
    diag_out$df <- p / (1 + eta_sel)

    if (!is.null(Sigma_ext)) {
      sr <- spectral_condition_ratio(Sigma_int, Sigma_ext)
      diag_out$kappa <- sr$kappa
      diag_out$eigenvalues_k <- sr$eigenvalues

      delta <- beta_ext - beta_int
      decomp <- decompose_bound(Sigma_int, Sigma_ext, delta,
                                 sigma2_int, sigma2_ext, n_int, n_ext)
      diag_out$signal_to_bias <- decomp$signal_to_bias
      diag_out$geometric_penalty <- decomp$geometric_penalty
      diag_out$eta_star <- decomp$eta_star
      diag_out$eta_meta <- decomp$eta_meta

      # ---- NEW: direction-specific biases and concordance analysis ----
      db <- estimate_direction_biases(delta, Sigma_int, Sigma_ext,
                                       sigma2_ext, n_ext)
      diag_out$b_vals <- db$b_vals
      diag_out$k_vals <- db$k_vals

      # Concordance check at the selected eta
      diag_out$concordance <- check_concordance(db$k_vals, db$b_vals,
                                                 sigma2_int %||% 1,
                                                 n_int %||% 100, eta_sel)

      # Direct numerical MSE comparison at selected eta
      adv_tbl <- spectral_advantage(db$k_vals, db$b_vals,
                                     sigma2_int %||% 1, n_int %||% 100,
                                     eta_sel)
      diag_out$advantage_at_eta <- adv_tbl$advantage[1]
      diag_out$mse_partial_at_eta <- adv_tbl$mse_partial[1]
      diag_out$mse_full_at_eta <- adv_tbl$mse_full[1]

      # Full advantage curve across eta_grid for plotting
      diag_out$advantage_curve <- spectral_advantage(db$k_vals, db$b_vals,
                                                      sigma2_int %||% 1,
                                                      n_int %||% 100,
                                                      eta_grid)
    }
  }

  # ===========================================================
  # Assemble output
  # ===========================================================
  out <- list(
    coefficients = coefficients,
    eta = eta_sel,
    eta_bound = eta_ub,
    scenario = scenario,
    tuning_method = tuning,
    tuning_result = tuning_result,
    eta_grid = eta_grid,
    diagnostics = diag_out,
    internal = list(beta = beta_int, Sigma = Sigma_int,
                    sigma2 = sigma2_int, n = n_int),
    external = list(beta = beta_ext, Sigma = Sigma_ext,
                    sigma2 = sigma2_ext, n = n_ext),
    ridge = list(lambda = lambda, used = !is.null(lambda)),
    p = p,
    call = cl
  )

  class(out) <- "mainger"
  out
}
