#include <RcppArmadillo.h>
arma::mat j_mvrnorm(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

arma::mat j_mvrnormFast(int n, arma::vec mu, arma::mat chol_sigma) {
  int ncols = chol_sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * chol_sigma;
}