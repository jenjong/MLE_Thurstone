#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]

SEXP rtn(double a, double b, double m, double s){
  
  double x;
  SEXP v;
  
  Function tr("rtruncnorm");

  v = tr(Named("n") = 1, Named("a") = a, Named("b") = b, 
     Named("mean")=m,Named("sd") = s);
  return v;
}
      
                
  


  
  