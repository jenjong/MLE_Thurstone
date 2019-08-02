# ifndef JBASE
# define JBASE

# include <RcppArmadillo.h>
arma::mat j_mvrnorm(int n, arma::vec mu, arma::mat sigma);
arma::mat j_mvrnormFast(int n, arma::vec mu, arma::mat chol_sigma);
# endif
