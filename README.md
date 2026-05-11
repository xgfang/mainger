# mainger

An R package for **M**ahalanobis-integrated **R**egression under privacy
constraints. Provides a unified framework for integrating external summary
statistics into internal regression models under three data-sharing regimes
(full, partial, restricted), with closed-form estimators, theoretical bounds
on the beneficial tuning range, data-adaptive tuning via cross-validation
and extended AIC, and diagnostics for assessing the privacy-precision
trade-off.

This is an anonymous release accompanying a submission to STAIX 2026.

## Installation

From the local tarball:

```r
install.packages("mainger_0.2.0.tar.gz", repos = NULL, type = "source")
```

Or directly from a clone of the source:

```r
# install.packages("remotes")
remotes::install_local("path/to/mainger")
```

System requirements: R ≥ 4.2, plus `MASS`, `glmnet`, `Rcpp`, `RcppArmadillo`.
Building requires a C++ compiler (Rtools on Windows, Xcode CLT on macOS,
g++ on Linux).

## Quick start

```r
library(mainger)

# Partial sharing: have internal X / Y, external coefficient only
fit <- mainger(X_int = X, Y_int = y, beta_ext = theta)
summary(fit)

# Full sharing: also have external Gram matrix
fit <- mainger(X_int = X, Y_int = y,
               beta_ext = theta, Sigma_ext = S,
               sigma2_ext = 1, n_ext = 5000)

# Restricted sharing: only marginal correlations and a reference panel
fit <- mainger(r_int = r, Sigma_ref = S_ref,
               beta_ext = theta, n_int = 2000)
```

See `vignette("quickstart", package = "mainger")` for a complete walkthrough
of all three regimes.

## Three data-sharing regimes

| Regime       | Internal access            | External access                                |
|--------------|----------------------------|------------------------------------------------|
| Full         | $(X_1, Y_1)$               | $\hat\theta$ and $\hat\Sigma_2$                |
| Partial      | $(X_1, Y_1)$               | $\hat\theta$ only                              |
| Restricted   | $r_1 = X_1^\top Y_1 / n_1$ | $\hat\theta$ and a reference $\Sigma_{\mathrm{ref}}$ |

The top-level `mainger()` function auto-detects the regime from which inputs
you supply, computes a theoretical bound on the tuning parameter, selects
the integration weight via CV (full / partial) or eAIC (restricted), and
returns an integrated coefficient vector with diagnostics.

## Reproducibility

Pass `cv_seed = 548` (or any fixed integer) to make the CV fold assignment
deterministic across machines:

```r
fit <- mainger(X_int = X, Y_int = y, beta_ext = theta,
               tuning = "cv", cv_seed = 548)
```

The accompanying agent (separate repository) uses this same convention
to produce bit-identical results when re-running an analysis.

## License

MIT.
