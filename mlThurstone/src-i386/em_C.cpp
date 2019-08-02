#include <RcppArmadillo.h>
#include "rtn.h"
#include "j_base.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List E_fun_C (int i,
           arma::umat pi_mat, 
           Rcpp::List rankIndex_list, 
           arma::colvec mu_e,
           arma::mat Sig_e,
           arma::mat Omg_e,
           int burn_num,
           int restore_num)
{

  int p = pi_mat.n_cols;
  int n = pi_mat.n_rows;
  colvec rank_index = rankIndex_list["rank_index"];
  colvec freq_v = rankIndex_list["freq_v"];
  colvec Ez = arma::zeros(n);
  mat Covz = arma::zeros(p, p);
  int rank_index_length = rank_index.n_rows;
  int x_idx =0;
  
  for (int j = 0; j<rank_index_length; j++)
  {
    if (rank_index(j) ==  i)
    {
     x_idx = j; 
     break;
    }
  }

  uvec x = arma::vectorise(pi_mat.row(x_idx));
  mat chol_sigma = arma::chol(Sig_e);
  mat chol_sigma_tmp = 5*arma::eye(p,p);
  mat z = arma::zeros(1,p);
  uvec zo(p), zo_tmp(p);
  zo = arma::sort_index(z, "descend");
  int init_iter = 1;  
  while(true)
  {
    if (init_iter>1000)
    {
      z = j_mvrnormFast(1, mu_e, chol_sigma_tmp);
    } else {
      z = j_mvrnormFast(1, mu_e, chol_sigma);
    }
    
    zo = sort_index(z, "descend") + 1;
    zo_tmp = zo - x;
    if (!any(zo_tmp)) break;
    init_iter++;
  }
  
  vec zv(p);
  for (int j=0; j<p;j++)
  {
    zv(j) = z(j);
  }

  int j;
  int j_c = 0;
  int j_k = 0;
  double m;
  double m1;
  double v;
  double a,b;
  // j_c : from 0 to (p-1)
  // j denotes location index in x
  //
  for (int iter=0; iter<(p*burn_num); iter++)
  {
    m1 = 0;
    j = x(iter%p)-1;
    for (int k=0; k<p; k++)
    {
      if (k==j) continue;
      m1 = m1 + Omg_e(j,k)*(zv(k) - mu_e(k));
    }
    m = mu_e(j) - m1/Omg_e(j,j);
    v = pow(1/Omg_e(j,j), 0.5);
     if (j_c==0)
    {
      b = R_PosInf;
    } else {
      b = zv(x(j_c-1)-1);
    }

    if (j_c==(p-1))
    {
      a = R_NegInf;
    }   else {
      a = zv(x(j_c+1)-1);
    }
    zv(j) = rtn(m, v,  a, b);
    j_c ++;
    if (j_c > (p-1)) j_c = 0;
  }
  
  
  vec zm = zeros<vec>(p);
  mat zc = zeros<mat>(p,p);
  
  
  for (int iter=0; iter<(p*restore_num); iter++)
  {
    m1 = 0;
    j = x(iter%p)-1;
    for (int k=0; k<p; k++)
    {
      if (k==j) continue;
      m1 = m1 + Omg_e(j,k)*(zv(k) - mu_e(k));
    }
    m = mu_e(j) - m1/Omg_e(j,j);
    v = pow(1/Omg_e(j,j), 0.5);
    if (j_c==0)
    {
      b = R_PosInf;
    } else {
      b = zv(x(j_c-1)-1);
    }
    
    if (j_c==(p-1))
    {
      a = R_NegInf;
    }   else {
      a = zv(x(j_c+1)-1);
    }
    zv(j) = rtn(m, v,  a, b);
    j_c ++;
    if (j_c > (p-1))
    {
      zm = zm + zv/(restore_num);
      zc = zc + zv*zv.t()/(restore_num);
      j_c = 0;
      j_k++;
    }
  }
  
  Ez = zm*freq_v(i-1);
  Covz = zc*freq_v(i-1);

  return Rcpp::List::create(Rcpp::Named("Ez")=Ez,
                            Rcpp::Named("Covz")=Covz);

  // return Rcpp::List::create(Rcpp::Named("Ez")=Ez,
  //                           Rcpp::Named("Covz")=Covz,
  //                           Rcpp::Named("Sig_e")=Sig_e,
  //                           Rcpp::Named("Omg_e")=Omg_e,
  //                           Rcpp::Named("aa")=rank_index,
  //                           Rcpp::Named("x_idx")=x_idx,
  //                           Rcpp::Named("x")=x, 
  //                           Rcpp::Named("chol_sigma")=chol_sigma,
  //                           Rcpp::Named("z")=z,
  //                           Rcpp::Named("zo")=zo,
  //                           Rcpp::Named("zv")=zv,
  //                           Rcpp::Named("a")=a,
  //                           Rcpp::Named("b")=b,
  //                           Rcpp::Named("m")=m,
  //                           Rcpp::Named("v")=v,
  //                           Rcpp::Named("zero")=zm,
  //                           Rcpp::Named("zmat")=zc
  //                           );
}
