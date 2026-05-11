# mainger 0.2.0

## Theoretical update: corrected Proposition 2

This release updates the package to reflect the corrected Proposition 2
on the full-vs-partial MSE comparison. The earlier claim that full sharing
always weakly dominates partial sharing at every eta was incorrect: the
eigenvalues k_j of the cross-covariance matrix are fixed by the data and
are not free parameters that the estimator can optimize.

The corrected result characterizes the MSE gap through the alignment
(concordance) between the eigenvalues k_j and the direction-specific
external biases b_j. Full sharing outperforms partial when large k_j are
paired with large b_j (the spectral structure shields the unreliable
directions); partial outperforms full when this alignment is reversed.

## New functions

- `spectral_advantage(k_vals, b_vals, sigma2_int, n_int, eta_grid)` evaluates
  the exact MSE gap A(eta) = MSE_partial(eta) - MSE_full(eta) at any eta,
  using the per-direction decomposition g_j(alpha) = (V + alpha^2 b_j) / (1 + alpha)^2.
- `check_concordance(k_vals, b_vals, sigma2_int, n_int, eta)` evaluates the
  sufficient conditions (k_j - 1)(k_j V - b_j eta) <= 0 (full wins) and
  (k_j - 1)(V - b_j eta) >= 0 (partial wins) from the corrected Proposition.
- `estimate_direction_biases(delta, Sigma_int, Sigma_ext, sigma2_ext, n_ext)`
  rotates the heterogeneity vector into the diagonalized basis to produce
  b_j estimates suitable for input to the two functions above.

## Updated functions

- `diagnose()` now reports the concordance verdict at the selected eta and
  shows the numerical A(eta) value. The old interpretation ("large kappa
  means full sharing strongly recommended") has been removed; it was
  mathematically unsupported.
- `mainger()` populates new diagnostic fields `concordance`, `advantage_at_eta`,
  `mse_partial_at_eta`, `mse_full_at_eta`, `advantage_curve`, `b_vals`, `k_vals`
  when full-sharing data are provided.
- `plot.mainger()` accepts `type = "advantage"` to display the A(eta) curve
  with clearly labeled "Full better" and "Partial better" regions.

## Tests

New tests verify the numerical examples from the corrected proof:
- Concordant case (k = 0.5, 2; b = 0.1, 10): MSE_full approx 1.43 at eta = 0.198
- Discordant case (k = 2, 0.5; b = 0.1, 10): MSE_full approx 2.147 at eta = 0.198
- All k_j = 1 produces zero spectral advantage (regardless of b_j values)

## Vignette updates

The full-sharing example in the quickstart vignette now walks through
`plot(fit, type = "advantage")` and `check_concordance()`. The privacy cost
interpretation table has been rewritten to remove the incorrect claim that
large kappa uniformly favors full sharing.


# mainger 0.1.0

Initial release.
