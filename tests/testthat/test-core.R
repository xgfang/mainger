test_that("partial estimator is convex combination", {
  set.seed(42)
  p <- 5
  b1 <- rnorm(p)
  b2 <- rnorm(p)

  for (eta in c(0, 0.5, 1, 5, 100)) {
    result <- est_partial(b1, b2, eta)
    expected <- b1 / (1 + eta) + eta * b2 / (1 + eta)
    expect_equal(result, expected, tolerance = 1e-12)
  }

  # eta = 0 recovers internal
  expect_equal(est_partial(b1, b2, 0), b1)
})

test_that("full estimator reduces to partial when Sigma_ext = c * Sigma_int", {
  set.seed(42)
  p <- 5
  A <- matrix(rnorm(p^2), p, p)
  Sigma_int <- crossprod(A) / p + 0.1 * diag(p)
  c_val <- 2.5
  Sigma_ext <- c_val * Sigma_int
  b1 <- rnorm(p)
  b2 <- rnorm(p)

  for (eta in c(0.5, 1, 3)) {
    full <- est_full(b1, Sigma_int, b2, Sigma_ext, eta)
    partial <- est_partial(b1, b2, c_val * eta)
    expect_equal(full, partial, tolerance = 1e-10)
  }
})

test_that("theoretical bound is positive", {
  set.seed(42)
  p <- 5
  b1 <- rnorm(p)
  b2 <- rnorm(p)
  A <- matrix(rnorm(p^2), p, p)
  Sigma <- crossprod(A) / p + 0.5 * diag(p)

  bound <- eta_bound_partial(b1, b2, Sigma, sigma2_int = 1)
  expect_true(bound > 0)
})

test_that("eAIC selects a reasonable eta", {
  set.seed(42)
  n <- 100; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(1, 0.5, 0, -0.5, 0.2)
  Y <- X %*% beta_true + rnorm(n, 0, 0.5)
  b_int <- as.numeric(solve(crossprod(X), crossprod(X, Y)))
  b_ext <- beta_true + rnorm(p, 0, 0.1)

  result <- tune_eaic(b_int, b_ext, X, Y, seq(0, 5, by = 0.1), p)
  expect_true(result$eta >= 0)
  expect_true(result$eta <= 5)
  expect_equal(length(result$coefficients), p)
})

test_that("zero_pad works correctly", {
  b_ext <- c(1, 2, 3)
  padded <- zero_pad(b_ext, 5)
  expect_equal(padded, c(1, 2, 3, 0, 0))

  padded2 <- zero_pad(b_ext, 5, shared_idx = c(2, 4, 5))
  expect_equal(padded2, c(0, 1, 0, 2, 3))
})

test_that("spectral condition ratio is 1 when Sigma_ext propto Sigma_int", {
  set.seed(42)
  p <- 5
  A <- matrix(rnorm(p^2), p, p)
  Sigma <- crossprod(A) / p + 0.5 * diag(p)

  sr <- spectral_condition_ratio(Sigma, 3 * Sigma)
  expect_equal(sr$kappa, 1, tolerance = 1e-8)
})

test_that("mainger auto-detects partial scenario", {
  set.seed(42)
  n <- 50; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- X %*% rnorm(p) + rnorm(n)
  b_ext <- rnorm(p)

  fit <- mainger(X_int = X, Y_int = Y, beta_ext = b_ext)
  expect_s3_class(fit, "mainger")
  expect_equal(fit$scenario, "partial")
  expect_equal(length(fit$coefficients), p)
  expect_true(fit$eta >= 0)
})


# ---------------------------------------------------------------------
# NEW TESTS: spectral advantage and concordance (from corrected Prop)
# ---------------------------------------------------------------------

test_that("spectral_advantage matches proof example (concordant case)", {
  # From the paper: V=1, b=(0.1, 10), k=(0.5, 2), eta=0.198
  # MSE_full = 0.521 + 0.909 = 1.430
  # MSE_partial at same eta is ~1.669 (at its optimum)
  adv <- spectral_advantage(k_vals = c(0.5, 2), b_vals = c(0.1, 10),
                            sigma2_int = 1, n_int = 1, eta_grid = 0.198)
  expect_equal(adv$mse_full, 1.430, tolerance = 0.01)
  expect_true(adv$advantage > 0)  # concordant: full wins
})

test_that("spectral_advantage matches proof example (discordant case)", {
  # Reversed eigenvalues: k=(2, 0.5) with same b=(0.1, 10)
  # MSE_full = 0.829 + 1.318 = 2.147, MSE_partial ~ 1.669
  adv <- spectral_advantage(k_vals = c(2, 0.5), b_vals = c(0.1, 10),
                            sigma2_int = 1, n_int = 1, eta_grid = 0.198)
  expect_equal(adv$mse_full, 2.147, tolerance = 0.01)
  expect_true(adv$advantage < 0)  # discordant: partial wins
})

test_that("check_concordance verdict matches discordant example", {
  cc <- check_concordance(k_vals = c(2, 0.5), b_vals = c(0.1, 10),
                          sigma2_int = 1, n_int = 1, eta = 0.198)
  expect_true(cc$partial_wins_sufficient)
  expect_false(cc$full_wins_sufficient)
})

test_that("check_concordance verdict is NEUTRAL when all k_j = 1", {
  cc <- check_concordance(k_vals = c(1, 1, 1), b_vals = c(0.1, 1, 10),
                          sigma2_int = 1, n_int = 1, eta = 0.5)
  # When k_j = 1, both conditions are trivially satisfied (product of zero)
  expect_true(cc$full_wins_sufficient)
  expect_true(cc$partial_wins_sufficient)
})

test_that("spectral_advantage is zero when all k_j = 1", {
  set.seed(42)
  adv <- spectral_advantage(k_vals = rep(1, 5), b_vals = runif(5, 0, 2),
                            sigma2_int = 1, n_int = 10,
                            eta_grid = seq(0, 5, length = 10))
  expect_true(all(abs(adv$advantage) < 1e-12))
})

test_that("estimate_direction_biases returns correct dimensions", {
  set.seed(42)
  p <- 5
  A <- matrix(rnorm(p^2), p, p)
  S1 <- crossprod(A) / p + 0.5 * diag(p)
  S2 <- 2 * S1 + 0.3 * diag(p)
  delta <- rnorm(p)

  db <- estimate_direction_biases(delta, S1, S2, sigma2_ext = 1, n_ext = 100)
  expect_equal(length(db$b_vals), p)
  expect_equal(length(db$k_vals), p)
  expect_true(all(db$b_vals >= 0))
})

test_that("diagnose runs without error on full-sharing fit", {
  set.seed(42)
  p <- 5
  A <- matrix(rnorm(p^2), p, p)
  S1 <- crossprod(A) / p + 0.5 * diag(p)
  S2 <- crossprod(matrix(rnorm(p^2), p, p)) / p + 0.5 * diag(p)
  b1 <- rnorm(p); b2 <- rnorm(p)

  fit <- mainger(beta_int = b1, Sigma_int = S1,
                 beta_ext = b2, Sigma_ext = S2,
                 sigma2_int = 1, sigma2_ext = 1,
                 n_int = 100, n_ext = 1000,
                 tuning = "fixed", eta = 1.0)

  expect_s3_class(fit, "mainger")
  expect_false(is.null(fit$diagnostics$concordance))
  expect_false(is.null(fit$diagnostics$advantage_at_eta))
  expect_output(diagnose(fit), "CONCORDANCE ANALYSIS")
})
