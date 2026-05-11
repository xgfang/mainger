// fast_linalg.cpp
// Rcpp/RcppArmadillo routines for the restricted-vs-partial simulation.
// Compile with: Rcpp::sourceCpp("fast_linalg.cpp")

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------
// 1. Fast generation of N(0, Sigma) samples via Cholesky factor L
//    Returns n x p matrix  X = Z * L  where Z ~ iid N(0,1)
// ---------------------------------------------------------------
// [[Rcpp::export]]
arma::mat fast_rmvnorm(int n, const arma::mat& L) {
  int p = L.n_cols;
  arma::mat Z(n, p, fill::randn);   // iid N(0,1)
  return Z * L;                      // each row ~ N(0, L'L = Sigma)
}

// ---------------------------------------------------------------
// 2. Fast X'X / n  (Gram matrix, scaled)
// ---------------------------------------------------------------
// [[Rcpp::export]]
arma::mat fast_gram(const arma::mat& X) {
  int n = X.n_rows;
  return X.t() * X / n;
}

// ---------------------------------------------------------------
// 3. Fast X'Y / n  (marginal correlations, scaled)
// ---------------------------------------------------------------
// [[Rcpp::export]]
arma::vec fast_xty(const arma::mat& X, const arma::vec& Y) {
  int n = X.n_rows;
  return X.t() * Y / n;
}

// ---------------------------------------------------------------
// 4. Fast ridge solve:  (A + lambda * I)^{-1} b
//    Uses Cholesky since A + lambda*I is SPD for lambda > 0
// ---------------------------------------------------------------
// [[Rcpp::export]]
arma::vec fast_ridge_solve(const arma::mat& A, const arma::vec& b, double lambda) {
  int p = A.n_cols;
  arma::mat M = A + lambda * arma::eye(p, p);
  return arma::solve(M, b, arma::solve_opts::likely_sympd);
}

// ---------------------------------------------------------------
// 5. Fast matrix-vector product  X * beta
// ---------------------------------------------------------------
// [[Rcpp::export]]
arma::vec fast_mv(const arma::mat& X, const arma::vec& beta) {
  return X * beta;
}

// ---------------------------------------------------------------
// 6. Fast predictive R^2
// ---------------------------------------------------------------
// [[Rcpp::export]]
double fast_pred_r2(const arma::vec& Y_true, const arma::vec& Y_pred) {
  double ss_res = arma::accu(arma::square(Y_true - Y_pred));
  double mn = arma::mean(Y_true);
  double ss_tot = arma::accu(arma::square(Y_true - mn));
  if (ss_tot < 1e-15) return 0.0;
  return 1.0 - ss_res / ss_tot;
}

// ---------------------------------------------------------------
// 7. Fast MSPE
// ---------------------------------------------------------------
// [[Rcpp::export]]
double fast_mspe(const arma::vec& Y_true, const arma::vec& Y_pred) {
  return arma::mean(arma::square(Y_true - Y_pred));
}

// ---------------------------------------------------------------
// 8. Run the full eta grid search for one replication
//    Avoids R-level looping over eta values entirely.
//
//    beta_int    : p-vector, internal baseline (ridge or summary-based)
//    theta_ext   : p-vector, external point estimate
//    X_val       : n_val x p validation matrix
//    Y_val       : n_val-vector validation response
//    eta_grid    : vector of candidate eta values
//
//    Returns: list(best_eta, best_r2, all_r2)
// ---------------------------------------------------------------
// [[Rcpp::export]]
List fast_eta_search(const arma::vec& beta_int,
                     const arma::vec& theta_ext,
                     const arma::mat& X_val,
                     const arma::vec& Y_val,
                     const arma::vec& eta_grid) {

  int n_eta = eta_grid.n_elem;

  // Precompute the two prediction vectors once
  arma::vec pred_int = X_val * beta_int;    // n_val x 1
  arma::vec pred_ext = X_val * theta_ext;   // n_val x 1

  double mn_y = arma::mean(Y_val);
  double ss_tot = arma::accu(arma::square(Y_val - mn_y));

  arma::vec all_r2(n_eta);
  double best_r2 = -datum::inf;
  int best_idx = 0;

  for (int i = 0; i < n_eta; i++) {
    double eta = eta_grid(i);
    double w1 = 1.0 / (1.0 + eta);
    double w2 = eta  / (1.0 + eta);

    // pred = w1 * pred_int + w2 * pred_ext  (no new matrix multiply!)
    arma::vec pred = w1 * pred_int + w2 * pred_ext;
    double ss_res = arma::accu(arma::square(Y_val - pred));
    double r2 = (ss_tot > 1e-15) ? (1.0 - ss_res / ss_tot) : 0.0;

    all_r2(i) = r2;
    if (r2 > best_r2) {
      best_r2 = r2;
      best_idx = i;
    }
  }

  return List::create(
    Named("best_eta") = eta_grid(best_idx),
    Named("best_r2")  = best_r2,
    Named("all_r2")   = all_r2
  );
}

// ---------------------------------------------------------------
// 9. One full replication: generate data, fit all methods, return results.
//    This avoids repeated R<->C++ overhead for the inner loop.
//
//    L_pop        : p x p upper Cholesky of population covariance
//    beta_int     : p-vector true internal effects
//    beta_ext     : p-vector true external effects
//    sigma2       : noise variance
//    n_int, n_ext, n_ref, n_test : sample sizes
//    lambda_ridge : ridge penalty for internal
//    lambda_ext   : ridge penalty for external (small)
//    eta_grid     : vector of candidate eta values
//    n_val        : validation set size (carved from test)
// ---------------------------------------------------------------
// [[Rcpp::export]]
List run_one_rep(const arma::mat& L_pop,
                 const arma::vec& beta_int,
                 const arma::vec& beta_ext,
                 double sigma2,
                 int n_int, int n_ext, int n_ref, int n_test,
                 double lambda_ridge, double lambda_ext,
                 const arma::vec& eta_grid,
                 int n_val) {

  int p = L_pop.n_cols;

  // --- Generate all data matrices ---
  arma::mat X_int  = fast_rmvnorm(n_int, L_pop);
  arma::mat X_ext  = fast_rmvnorm(n_ext, L_pop);
  arma::mat X_ref  = fast_rmvnorm(n_ref, L_pop);
  arma::mat X_test = fast_rmvnorm(n_test, L_pop);

  double sd_noise = std::sqrt(sigma2);
  arma::vec Y_int  = X_int  * beta_int + sd_noise * arma::randn(n_int);
  arma::vec Y_ext  = X_ext  * beta_ext + sd_noise * arma::randn(n_ext);
  arma::vec Y_test = X_test * beta_int + sd_noise * arma::randn(n_test);

  // --- Sufficient statistics ---
  arma::mat Sigma_int = X_int.t() * X_int / n_int;
  arma::vec r_int     = X_int.t() * Y_int / n_int;

  arma::mat Sigma_ext = X_ext.t() * X_ext / n_ext;
  arma::vec r_ext     = X_ext.t() * Y_ext / n_ext;

  arma::mat Sigma_ref = X_ref.t() * X_ref / n_ref;

  // --- External point estimate (ridge) ---
  arma::vec theta_ext = arma::solve(Sigma_ext + lambda_ext * arma::eye(p, p),
                                     r_ext,
                                     arma::solve_opts::likely_sympd);

  // --- Internal baselines ---
  // Ridge on individual-level data
  arma::vec beta_ridge = arma::solve(Sigma_int + lambda_ridge * arma::eye(p, p),
                                      r_int,
                                      arma::solve_opts::likely_sympd);
  // Summary-stat based (uses reference LD)
  arma::vec beta_summary = arma::solve(Sigma_ref + lambda_ridge * arma::eye(p, p),
                                        r_int,
                                        arma::solve_opts::likely_sympd);

  // --- Split test into validation + evaluation ---
  arma::mat X_val  = X_test.rows(0, n_val - 1);
  arma::vec Y_val  = Y_test.subvec(0, n_val - 1);
  arma::mat X_eval = X_test.rows(n_val, n_test - 1);
  arma::vec Y_eval = Y_test.subvec(n_val, n_test - 1);

  // --- Eta search for Partial and Restricted ---
  List res_partial    = fast_eta_search(beta_ridge, theta_ext, X_val, Y_val, eta_grid);
  List res_restricted = fast_eta_search(beta_summary, theta_ext, X_val, Y_val, eta_grid);

  double best_eta_p = as<double>(res_partial["best_eta"]);
  double best_eta_r = as<double>(res_restricted["best_eta"]);

  // --- Evaluate on held-out data ---
  double w1_p = 1.0 / (1.0 + best_eta_p);
  double w2_p = best_eta_p / (1.0 + best_eta_p);
  arma::vec beta_p_opt = w1_p * beta_ridge + w2_p * theta_ext;

  double w1_r = 1.0 / (1.0 + best_eta_r);
  double w2_r = best_eta_r / (1.0 + best_eta_r);
  arma::vec beta_r_opt = w1_r * beta_summary + w2_r * theta_ext;

  arma::vec pred_ridge   = X_eval * beta_ridge;
  arma::vec pred_summary = X_eval * beta_summary;
  arma::vec pred_ext     = X_eval * theta_ext;
  arma::vec pred_partial = X_eval * beta_p_opt;
  arma::vec pred_restr   = X_eval * beta_r_opt;

  double r2_ridge    = fast_pred_r2(Y_eval, pred_ridge);
  double r2_summary  = fast_pred_r2(Y_eval, pred_summary);
  double r2_ext      = fast_pred_r2(Y_eval, pred_ext);
  double r2_partial  = fast_pred_r2(Y_eval, pred_partial);
  double r2_restr    = fast_pred_r2(Y_eval, pred_restr);

  double mspe_ridge   = fast_mspe(Y_eval, pred_ridge);
  double mspe_summary = fast_mspe(Y_eval, pred_summary);
  double mspe_ext     = fast_mspe(Y_eval, pred_ext);
  double mspe_partial = fast_mspe(Y_eval, pred_partial);
  double mspe_restr   = fast_mspe(Y_eval, pred_restr);

  return List::create(
    Named("r2_ridge")      = r2_ridge,
    Named("r2_summary")    = r2_summary,
    Named("r2_ext")        = r2_ext,
    Named("r2_partial")    = r2_partial,
    Named("r2_restricted") = r2_restr,
    Named("mspe_ridge")    = mspe_ridge,
    Named("mspe_summary")  = mspe_summary,
    Named("mspe_ext")      = mspe_ext,
    Named("mspe_partial")  = mspe_partial,
    Named("mspe_restricted") = mspe_restr,
    Named("eta_partial")   = best_eta_p,
    Named("eta_restricted") = best_eta_r
  );
}
