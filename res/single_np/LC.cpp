// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec LC(arma::mat xc, arma::vec yc, arma::mat eval,  arma::vec h)
{
  int n_y = yc.n_elem;
  int n_eval = eval.n_rows;
  int J = eval.n_cols;
  int i,j,i1;

  arma::vec est;
  est.zeros(n_eval);
  arma::vec Kvec;
  arma::vec Kvec1;
  arma::vec idx;
  idx.zeros(n_y);
  arma::vec z1;
  arma::vec z2;
  arma::vec z;
  double t1;

  for(i = 0; i < n_eval; i++){
    Kvec.ones(n_y);
    
    for(j = 0; j < J; j++){
      z1 = xc.col(j)/h(j);
      z2.ones(n_y);
      z2 = z2 * (eval(i,j)/h(j));
      z = z1-z2;
      for(i1 = 0; i1 < n_y; i1++){
        if(z(i1)<1 && z(i1)>-1){
          idx(i1) = 1;
        }else{
          idx(i1) = 0;
        }
      }
      arma::vec t;
      t = 1.0 - z%z;
      Kvec1 = 0.75*t%idx/h(j);
      Kvec = Kvec%Kvec1;
    }

    arma::vec multi;
    multi = Kvec%yc;
    t1 = sum(multi)/sum(Kvec);
    // t1 = sum(multi);
    est(i) = t1;
  }

  
  return est;
}







