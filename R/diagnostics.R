#' @title Diagnostic Functions for Privacy-Precision Trade-off
#' @name diagnostics
#' @description
#' Diagnostics are based on the corrected Proposition (see `prop2_corrected`):
#' the sign of the full-vs-partial MSE gap depends on whether the spectral
#' structure (eigenvalues k_j) is concordant or discordant with the direction-
#' specific heterogeneity (b_j). Large kappa alone does NOT imply full sharing
#' is better; it only means the spectral structure is informative.
NULL

#' Spectral condition ratio
#'
#' Computes the condition ratio of the cross-covariance eigenstructure:
#' kappa = k_max / k_min where k_j are eigenvalues of
#' Sigma_int^{1/2} Sigma_ext^{-1} Sigma_int^{1/2}.
#' When kappa = 1, the privacy cost of partial vs full sharing is exactly zero
#' (Proposition 1). When kappa > 1, the full-vs-partial comparison depends on
#' how the eigenvalues align with the direction-specific bias (use
#' `spectral_advantage` to evaluate).
#'
#' @param Sigma_int Internal Gram matrix (p x p)
#' @param Sigma_ext External Gram matrix (p x p)
#' @return List with:
#'   \itemize{
#'     \item `kappa`: the condition ratio k_max / k_min
#'     \item `eigenvalues`: sorted vector of k_j
#'     \item `eigenvectors`: orthonormal basis vectors (columns)
#'     \item `k_min`, `k_max`: extremes
#'     \item `Sigma_int_half`: the matrix square root of Sigma_int
#'   }
#' @export
spectral_condition_ratio <- function(Sigma_int, Sigma_ext) {
  eig_int <- eigen(Sigma_int, symmetric = TRUE)
  Sigma_int_half <- eig_int$vectors %*% diag(sqrt(pmax(eig_int$values, 0))) %*%
    t(eig_int$vectors)
  cross <- Sigma_int_half %*% solve(Sigma_ext) %*% Sigma_int_half
  cross_eig <- eigen(cross, symmetric = TRUE)
  k_vals <- pmax(cross_eig$values, 0)
  ord <- order(k_vals)

  list(
    kappa = max(k_vals) / max(min(k_vals), 1e-15),
    eigenvalues = k_vals[ord],
    eigenvectors = cross_eig$vectors[, ord, drop = FALSE],
    k_min = min(k_vals),
    k_max = max(k_vals),
    Sigma_int_half = Sigma_int_half
  )
}

#' Spectral advantage function A(eta)
#'
#' Computes A(eta) = MSE_partial(eta) - MSE_full(eta) using the per-direction
#' decomposition from Proposition 2 (corrected). When A > 0, full sharing
#' outperforms partial; when A < 0, partial outperforms full.
#'
#' Given estimates of the direction-specific biases b_j (expected squared
#' difference between theta_j and beta_1j in the diagonalized basis),
#' this function evaluates the exact MSE gap using
#' g_j(alpha) = (V + alpha^2 * b_j) / (1 + alpha)^2, where V = sigma2_int / n_int.
#'
#' @param k_vals Vector of eigenvalues k_j (from `spectral_condition_ratio`)
#' @param b_vals Vector of direction-specific external risks b_j (same length as k_vals)
#' @param sigma2_int Internal error variance
#' @param n_int Internal sample size
#' @param eta_grid Vector of eta values at which to evaluate A(eta)
#' @return A data.frame with columns eta, mse_partial, mse_full, advantage
#' @export
spectral_advantage <- function(k_vals, b_vals, sigma2_int, n_int, eta_grid) {
  V <- sigma2_int / n_int
  g_j <- function(alpha, b) (V + alpha^2 * b) / (1 + alpha)^2

  p <- length(k_vals)
  if (length(b_vals) != p) stop("k_vals and b_vals must have the same length.")

  out <- data.frame(eta = eta_grid, mse_partial = NA_real_,
                    mse_full = NA_real_, advantage = NA_real_)

  for (i in seq_along(eta_grid)) {
    eta <- eta_grid[i]
    mse_part <- sum(sapply(b_vals, function(b) g_j(eta, b)))
    mse_full <- sum(mapply(function(k, b) g_j(eta / k, b), k_vals, b_vals))
    out$mse_partial[i] <- mse_part
    out$mse_full[i]    <- mse_full
    out$advantage[i]   <- mse_part - mse_full
  }

  out
}

#' Check concordance between spectral structure and heterogeneity
#'
#' Evaluates the sufficient conditions from Proposition 2 (corrected):
#' full wins if (k_j - 1)(k_j * V - b_j * eta) <= 0 for all j, and
#' partial wins if (k_j - 1)(V - b_j * eta) >= 0 for all j.
#'
#' @param k_vals Eigenvalues from `spectral_condition_ratio`
#' @param b_vals Direction-specific biases in the diagonalized basis
#' @param sigma2_int Internal variance
#' @param n_int Internal sample size
#' @param eta The value of eta at which to check conditions
#' @return List with:
#'   \itemize{
#'     \item `full_wins_sufficient`: TRUE if full-wins condition holds for all j
#'     \item `partial_wins_sufficient`: TRUE if partial-wins condition holds for all j
#'     \item `concordance_per_direction`: character vector ("full", "partial", "indeterminate")
#'     \item `verdict`: overall assessment
#'     \item `conditions_full`, `conditions_partial`: the raw condition values per j
#'   }
#' @export
check_concordance <- function(k_vals, b_vals, sigma2_int, n_int, eta) {
  V <- sigma2_int / n_int
  cond_full <- (k_vals - 1) * (k_vals * V - b_vals * eta)
  cond_part <- (k_vals - 1) * (V - b_vals * eta)

  per_dir <- ifelse(cond_full <= 0 & cond_part <= 0, "full",
             ifelse(cond_full >= 0 & cond_part >= 0, "partial",
             "indeterminate"))

  full_all <- all(cond_full <= 1e-10)
  part_all <- all(cond_part >= -1e-10)

  verdict <- if (full_all && !part_all) {
    "CONCORDANT: full sharing is guaranteed to outperform partial at this eta."
  } else if (part_all && !full_all) {
    "DISCORDANT: partial sharing is guaranteed to outperform full at this eta."
  } else if (full_all && part_all) {
    "NEUTRAL: the eigenvalues are effectively all 1; full equals partial."
  } else {
    "INDETERMINATE: sufficient conditions for neither direction hold uniformly. Use spectral_advantage() for the exact numerical comparison."
  }

  list(
    full_wins_sufficient = full_all,
    partial_wins_sufficient = part_all,
    concordance_per_direction = per_dir,
    verdict = verdict,
    conditions_full = cond_full,
    conditions_partial = cond_part
  )
}

#' Decompose the beneficial range into Signal-to-Bias and Geometric Penalty
#'
#' Implements the decomposition of the full-sharing beneficial range bound:
#' eta_star = (Signal-to-Bias Adjustment) x (Geometric Mismatch Penalty) <= eta_meta.
#'
#' This is a sufficient condition on eta for the full estimator to improve over
#' the internal-only estimator. It does NOT resolve the full-vs-partial comparison
#' (use `spectral_advantage` or `check_concordance` for that).
#'
#' @param Sigma_int Internal Gram matrix
#' @param Sigma_ext External Gram matrix
#' @param delta Heterogeneity vector (beta_ext - beta_int)
#' @param sigma2_int Internal error variance
#' @param sigma2_ext External error variance
#' @param n_int Internal sample size
#' @param n_ext External sample size
#' @return List with decomposition components
#' @export
decompose_bound <- function(Sigma_int, Sigma_ext, delta,
                            sigma2_int, sigma2_ext, n_int, n_ext) {
  eig_int <- sort(eigen(Sigma_int, symmetric = TRUE, only.values = TRUE)$values)
  eig_ext <- sort(eigen(Sigma_ext, symmetric = TRUE, only.values = TRUE)$values)

  g_p <- max(eig_int)
  t_1 <- min(eig_ext)

  eig_int_full <- eigen(Sigma_int, symmetric = TRUE)
  Sigma_int_half <- eig_int_full$vectors %*%
    diag(sqrt(pmax(eig_int_full$values, 0))) %*% t(eig_int_full$vectors)
  cross <- Sigma_int_half %*% solve(Sigma_ext) %*% Sigma_int_half
  k_vals <- sort(eigen(cross, symmetric = TRUE, only.values = TRUE)$values)
  k_1 <- min(k_vals)

  delta_sq <- sum(delta^2)
  eta_meta <- (sigma2_int / n_int) / (sigma2_ext / n_ext)
  signal_to_bias <- (sigma2_int / n_int) / (delta_sq + sigma2_ext / n_ext)
  geometric_penalty <- k_1 * t_1 / g_p

  list(
    signal_to_bias = signal_to_bias,
    geometric_penalty = geometric_penalty,
    eta_star = signal_to_bias * geometric_penalty,
    eta_meta = eta_meta
  )
}

#' Estimate direction-specific biases from data
#'
#' Rotates the heterogeneity vector into the diagonalized basis to produce
#' b_j estimates. Specifically, b_j = <v_j, Sigma_int^{1/2} delta>^2 +
#' sigma2_ext / n_ext, where v_j are the eigenvectors of the cross matrix.
#'
#' @param delta Heterogeneity vector in the original basis
#' @param Sigma_int Internal Gram matrix
#' @param Sigma_ext External Gram matrix
#' @param sigma2_ext External error variance (NULL to use bias-only estimate)
#' @param n_ext External sample size
#' @return List with b_vals (direction-specific biases), k_vals (eigenvalues),
#'   and delta_rotated (delta projected onto the diagonalized basis)
#' @export
estimate_direction_biases <- function(delta, Sigma_int, Sigma_ext,
                                      sigma2_ext = NULL, n_ext = NULL) {
  sr <- spectral_condition_ratio(Sigma_int, Sigma_ext)
  delta_rot <- as.numeric(t(sr$eigenvectors) %*% (sr$Sigma_int_half %*% delta))
  b_bias <- delta_rot^2

  if (!is.null(sigma2_ext) && !is.null(n_ext)) {
    b_vals <- b_bias + sigma2_ext / n_ext
  } else {
    b_vals <- b_bias
  }

  list(b_vals = b_vals, k_vals = sr$eigenvalues, delta_rotated = delta_rot)
}

#' Run full diagnostics on a mainger fit
#'
#' Produces a diagnostic report including the spectral condition ratio,
#' concordance analysis, direct MSE comparison at the selected eta, and the
#' Signal-to-Bias decomposition. The concordance analysis replaces the old
#' interpretation that large kappa implies full is better.
#'
#' @param fit A mainger object
#' @return Invisibly returns the diagnostics list
#' @export
diagnose <- function(fit) {
  if (!inherits(fit, "mainger")) stop("Input must be a mainger object.")

  cat("Privacy-Precision Diagnostic Report\n")
  cat(strrep("-", 50), "\n")
  cat(sprintf("Scenario:           %s sharing\n", fit$scenario))
  cat(sprintf("Selected eta:       %.4f\n", fit$eta))
  cat(sprintf("Theoretical bound:  %.4f\n", fit$eta_bound))
  cat(sprintf("Effective df:       %.1f / %d\n", fit$p / (1 + fit$eta), fit$p))

  d <- fit$diagnostics
  if (is.null(d)) {
    cat("\n(No diagnostics available.)\n")
    return(invisible(NULL))
  }

  if (!is.null(d$kappa)) {
    cat("\nSPECTRAL STRUCTURE\n")
    cat(strrep("-", 50), "\n")
    cat(sprintf("Condition ratio kappa: %.3f\n", d$kappa))
    if (d$kappa < 1.1) {
      cat("  Sigma_2 is approximately proportional to Sigma_1.\n")
      cat("  By Proposition 1, partial sharing is nearly as good as full:\n")
      cat("  the privacy cost of not sharing Sigma_2 is negligible.\n")
    } else {
      cat("  The spectral structure is informative. Whether full sharing\n")
      cat("  outperforms partial depends on alignment between the eigenvalues\n")
      cat("  and the direction-specific bias (see concordance analysis).\n")
    }
  }

  if (!is.null(d$concordance)) {
    cat("\nCONCORDANCE ANALYSIS (at selected eta)\n")
    cat(strrep("-", 50), "\n")
    cat(d$concordance$verdict, "\n", sep = "")
    per_dir <- d$concordance$concordance_per_direction
    n_full <- sum(per_dir == "full")
    n_part <- sum(per_dir == "partial")
    n_indet <- length(per_dir) - n_full - n_part
    cat(sprintf("  Directions where full is preferred:    %d / %d\n",
                n_full, length(per_dir)))
    cat(sprintf("  Directions where partial is preferred: %d / %d\n",
                n_part, length(per_dir)))
    if (n_indet > 0)
      cat(sprintf("  Indeterminate directions:              %d / %d\n",
                  n_indet, length(per_dir)))
  }

  if (!is.null(d$advantage_at_eta)) {
    cat("\nSPECTRAL ADVANTAGE (numerical)\n")
    cat(strrep("-", 50), "\n")
    a <- d$advantage_at_eta
    cat(sprintf("MSE_partial - MSE_full at selected eta:  %+.5f\n", a))
    if (!is.null(d$mse_full_at_eta)) {
      if (a > 1e-8) {
        cat(sprintf("  Full sharing outperforms partial by %.2f%%.\n",
                    100 * a / d$mse_full_at_eta))
      } else if (a < -1e-8) {
        cat(sprintf("  Partial sharing outperforms full by %.2f%%.\n",
                    100 * abs(a) / d$mse_full_at_eta))
      } else {
        cat("  Full and partial give equal MSE at this eta.\n")
      }
    }
  }

  if (!is.null(d$signal_to_bias)) {
    cat("\nBENEFICIAL RANGE DECOMPOSITION\n")
    cat(strrep("-", 50), "\n")
    cat(sprintf("Signal-to-Bias Adjustment:    %.4f\n", d$signal_to_bias))
    cat(sprintf("Geometric Mismatch Penalty:   %.4f\n", d$geometric_penalty))
    cat(sprintf("Our bound (eta_star):         %.4f\n", d$eta_star))
    cat(sprintf("Meta-analysis weight:         %.4f\n", d$eta_meta))
  }

  cat("\n")
  invisible(d)
}
