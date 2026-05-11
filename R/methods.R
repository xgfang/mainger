#' @title Plot method for mainger objects
#' @param x A mainger object
#' @param type One of "eta_curve" (default), "diagnostics", or "advantage"
#' @param ... Additional arguments (not used)
#' @export
plot.mainger <- function(x, type = c("eta_curve", "diagnostics", "advantage"), ...) {
  type <- match.arg(type)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install it with install.packages('ggplot2').")
  }

  if (type == "eta_curve") {
    plot_eta_curve(x)
  } else if (type == "advantage") {
    plot_advantage(x)
  } else {
    plot_diagnostics(x)
  }
}

#' @keywords internal
plot_eta_curve <- function(fit) {
  # Get eta-MSE path from tuning result
  if (!is.null(fit$tuning_result$cv_errors)) {
    path_df <- fit$tuning_result$cv_errors
    colnames(path_df) <- c("eta", "criterion")
    ylab <- "CV error"
  } else if (!is.null(fit$tuning_result$eaic_values)) {
    path_df <- fit$tuning_result$eaic_values
    colnames(path_df) <- c("eta", "criterion")
    ylab <- "eAIC"
  } else {
    message("No tuning path available for plotting.")
    return(invisible(NULL))
  }

  p <- ggplot2::ggplot(path_df, ggplot2::aes(x = eta, y = criterion)) +
    ggplot2::geom_line(color = "#1F78B4", linewidth = 1) +
    ggplot2::geom_vline(xintercept = fit$eta_bound, linetype = "dotted",
                        color = "grey40", linewidth = 0.5) +
    ggplot2::geom_point(data = data.frame(eta = fit$eta,
                                           criterion = path_df$criterion[
                                             which.min(abs(path_df$eta - fit$eta))]),
                        ggplot2::aes(x = eta, y = criterion),
                        color = "#1F78B4", size = 3) +
    ggplot2::annotate("text", x = fit$eta_bound, y = max(path_df$criterion) * 0.98,
                      label = expression(eta[bound]^"*"),
                      hjust = -0.3, size = 3.5, color = "grey40") +
    ggplot2::labs(x = expression(eta), y = ylab,
                  title = sprintf("Tuning curve (%s sharing, %s)",
                                  fit$scenario, fit$tuning_method)) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  print(p)
  invisible(p)
}

#' @keywords internal
plot_diagnostics <- function(fit) {
  if (is.null(fit$diagnostics$eigenvalues_k)) {
    message("Spectral diagnostics require full-sharing scenario with Sigma_ext.")
    return(invisible(NULL))
  }

  k_vals <- fit$diagnostics$eigenvalues_k
  p <- length(k_vals)
  df <- data.frame(index = 1:p, k = k_vals)

  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = index, y = k)) +
    ggplot2::geom_point(color = "#7F77DD", size = 1.5) +
    ggplot2::geom_hline(yintercept = mean(k_vals), linetype = "dashed",
                        color = "grey50") +
    ggplot2::labs(x = "Eigenvalue index j",
                  y = expression(k[j]),
                  title = sprintf("Spectral structure (kappa = %.2f)",
                                  fit$diagnostics$kappa)) +
    ggplot2::theme_bw(base_size = 11)

  print(p1)
  invisible(p1)
}

#' @keywords internal
plot_advantage <- function(fit) {
  adv <- fit$diagnostics$advantage_curve
  if (is.null(adv)) {
    message("Advantage curve requires full-sharing scenario with Sigma_ext.")
    return(invisible(NULL))
  }

  # Highlight selected eta
  sel_idx <- which.min(abs(adv$eta - fit$eta))
  sel_pt <- adv[sel_idx, ]

  p <- ggplot2::ggplot(adv, ggplot2::aes(x = eta, y = advantage)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
    ggplot2::geom_line(color = "#378ADD", linewidth = 1.1) +
    ggplot2::geom_point(data = sel_pt,
                        ggplot2::aes(x = eta, y = advantage),
                        color = "#378ADD", size = 3) +
    ggplot2::geom_vline(xintercept = fit$eta_bound, linetype = "dotted",
                        color = "grey40") +
    ggplot2::annotate("text",
                      x = max(adv$eta) * 0.05,
                      y = max(adv$advantage),
                      label = "Full better",
                      color = "#1D9E75", hjust = 0, vjust = 1, size = 3.5) +
    ggplot2::annotate("text",
                      x = max(adv$eta) * 0.05,
                      y = min(adv$advantage),
                      label = "Partial better",
                      color = "#A32D2D", hjust = 0, vjust = 0, size = 3.5) +
    ggplot2::labs(x = expression(eta),
                  y = expression(A(eta) == MSE[partial] - MSE[full]),
                  title = "Spectral advantage curve") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  print(p)
  invisible(p)
}


#' @title Summary method for mainger objects
#' @param object A mainger object
#' @param n_top Number of top coefficients to show (default 10)
#' @param ... Additional arguments (not used)
#' @export
summary.mainger <- function(object, n_top = 10, ...) {
  cat("\nMahalanobis-Integrated Linear Regression\n")
  cat(strrep("-", 42), "\n")
  cat(sprintf("Scenario:           %s sharing\n", object$scenario))
  cat(sprintf("Tuning method:      %s\n", object$tuning_method))
  cat(sprintf("Selected eta:       %.4f\n", object$eta))
  cat(sprintf("Theoretical bound:  %.4f\n", object$eta_bound))
  cat(sprintf("Dimensions:         p = %d\n", object$p))
  cat(sprintf("Effective df:       %.1f\n", object$p / (1 + object$eta)))

  if (!is.null(object$internal$n)) {
    cat(sprintf("Internal n:         %d\n", object$internal$n))
  }
  if (!is.null(object$external$n)) {
    cat(sprintf("External n:         %d\n", object$external$n))
  }

  if (object$ridge$used) {
    cat(sprintf("Ridge lambda:       %.4f\n", object$ridge$lambda))
  }

  # Coefficient comparison
  if (!is.null(object$internal$beta)) {
    cat(sprintf("\nCoefficients (top %d by |change| from internal):\n", n_top))
    b_int <- object$internal$beta
    b_ext <- object$external$beta
    b_out <- object$coefficients
    change <- abs(b_out - b_int) / (abs(b_int) + 1e-8)
    top_idx <- head(order(change, decreasing = TRUE), n_top)

    coef_df <- data.frame(
      Internal = round(b_int[top_idx], 4),
      External = round(b_ext[top_idx], 4),
      Integrated = round(b_out[top_idx], 4),
      Pct_Change = sprintf("%+.1f%%", 100 * (b_out[top_idx] - b_int[top_idx]) /
                             (abs(b_int[top_idx]) + 1e-8))
    )
    if (!is.null(names(b_int))) {
      rownames(coef_df) <- names(b_int)[top_idx]
    }
    print(coef_df)
  }

  invisible(object)
}

#' @title Print method for mainger objects
#' @param x A mainger object
#' @param ... Additional arguments
#' @export
print.mainger <- function(x, ...) {
  cat(sprintf("mainger fit (%s sharing, eta = %.3f, p = %d)\n",
              x$scenario, x$eta, x$p))
  invisible(x)
}

#' @title Predict method for mainger objects
#' @param object A mainger object
#' @param newdata New design matrix (n_new x p)
#' @param ... Additional arguments
#' @return Predicted response vector (n_new x 1)
#' @export
predict.mainger <- function(object, newdata, ...) {
  if (ncol(newdata) != object$p) {
    stop(sprintf("newdata has %d columns, expected %d.", ncol(newdata), object$p))
  }
  as.numeric(newdata %*% object$coefficients)
}
