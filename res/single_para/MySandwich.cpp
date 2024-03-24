// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat sandwich(arma::mat X, arma::mat A, arma::mat A_bar, arma::mat Sigma_hat, arma::vec Nc_group, const int R, const int N, const int spill, const int design)
{
  int i1;
  int i2;

  arma::mat V_hat;
  V_hat.zeros(R,R);

  arma::mat xc_left;
  arma::mat xc_right;
  arma::mat v;
  arma::vec e;
  e.ones(N);

  xc_left.zeros(N,4);
  xc_right.zeros(N,4);
  xc_left.col(0) = e;
  xc_right.col(0) = e;

  arma::mat t;

  for(i1 = 0; i1 < R; i1++)
    for(i2 = 0; i2 < R; i2++){

          t = X.row(i1);
          xc_left.col(1) = trans(t);
          t = A.row(i1);
          xc_left.col(2) = trans(t);
          t = X.row(i2);
          xc_right.col(1) = trans(t);
          t = A.row(i2);
          xc_right.col(2) = trans(t);

          if(design == 0 || spill == 0)
          {
            v = (Sigma_hat(i1,i2)*inv(trans(xc_left.submat(0,0,N-1,2))*xc_left.submat(0,0,N-1,2))
                   *trans(xc_left.submat(0,0,N-1,2))*xc_right.submat(0,0,N-1,2)
                   *inv(trans(xc_right.submat(0,0,N-1,2))*xc_right.submat(0,0,N-1,2)));
            V_hat(i1,i2) = v(2,2);
          }
          if(spill == 1 && design == 1)
          {
            t = A_bar.row(i1);
            xc_left.col(3) = trans(t);
            t = A_bar.row(i2);
            xc_right.col(3) = trans(t);
            v = (Sigma_hat(i1,i2)*inv(trans(xc_left)*xc_left)*trans(xc_left)*xc_right)*inv(trans(xc_right)*xc_right);
            V_hat(i1,i2) = sum(sum(v.submat(2,2,3,3),1));
          }
          if(spill == 1 && design == 2)
          {
            t = A_bar.row(i1);
            xc_left.col(3) = trans(t);
            t = A_bar.row(i2);
            xc_right.col(3) = trans(t);
            if(Nc_group(i1)==1 && Nc_group(i2)==1)
            {
              v = (Sigma_hat(i1,i2)*inv(trans(xc_left.submat(0,0,N-1,2))*xc_left.submat(0,0,N-1,2))
                     *trans(xc_left.submat(0,0,N-1,2))*xc_right.submat(0,0,N-1,2)
                     *inv(trans(xc_right.submat(0,0,N-1,2))*xc_right.submat(0,0,N-1,2)));
              V_hat(i1,i2) = v(2,2);
            }
            if(Nc_group(i1)>1 && Nc_group(i2)==1)
            {
              v = (Sigma_hat(i1,i2)*inv(trans(xc_left)*xc_left)
                     *trans(xc_left)*xc_right.submat(0,0,N-1,2)
                     *inv(trans(xc_right.submat(0,0,N-1,2))*xc_right.submat(0,0,N-1,2)));
              V_hat(i1,i2) = sum(sum(v.submat(2,2,3,2),1));
            }
            if(Nc_group(i1)==1 && Nc_group(i2)>1)
            {
              v = (Sigma_hat(i1,i2)*inv(trans(xc_left.submat(0,0,N-1,2))*xc_left.submat(0,0,N-1,2))
                     *trans(xc_left.submat(0,0,N-1,2))*xc_right
                     *inv(trans(xc_right)*xc_right));
              V_hat(i1,i2) = sum(sum(v.submat(2,2,2,3),1));
            }
            if(Nc_group(i1)>1 && Nc_group(i2)>1)
            {
              v = (Sigma_hat(i1,i2)*inv(trans(xc_left)*xc_left)*trans(xc_left)*xc_right)*inv(trans(xc_right)*xc_right);
              V_hat(i1,i2) = sum(sum(v.submat(2,2,3,3),1));
            }
              
          }

        }

  return V_hat;
}
