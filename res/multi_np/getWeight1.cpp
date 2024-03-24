// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::cube getWeight(arma::cube X, arma::cube A, arma::vec idx, arma::vec ww, const int target_policy, const int d_w)
{
  int R = X.n_rows;
  int M = X.n_cols;
  int N = X.n_slices;
  int r, i1, i, t1, t2, j1, j2, j, t;

  arma::mat Stilde0, Stilde, beta, beta1;
  beta.zeros(R,d_w);
  arma::mat G, G1;
  double kappa, x_c1, x_c2, x_c1_next, x_c2_next;
  double pi = 3.141592654;
  arma::vec bases_eval0, bases_eval1, idx1;
  arma::vec delta1, delta2;
  arma::cube w;
  w.zeros(R,M,N/2);
  arma::vec eigval;
  arma::mat eigvec;

  
  // estimate on the first subsample (idx)
  for(r=0; r < R; r++){

    G.zeros(d_w,d_w);
    G1.zeros(d_w,d_w);
    Stilde0 = X.row(r);

    for(i1=0; i1 < N/2; i1++){

      i = idx(i1)-1;

      // Stilde = (X(r,:,i)+X_bar(r,:,i))/2;
      Stilde = Stilde0.col(i);

      for(t1=0; t1 < M-1; t1++)
        for(t2=0; t2 < M-1; t2++){

          x_c1 = Stilde(t1);
          x_c2 = Stilde(t2);
          x_c1_next = Stilde(t1+1);
          x_c2_next = Stilde(t2+1);

          bases_eval0.ones(d_w);
          bases_eval1.ones(d_w);
          for(j1=0; j1<(d_w-1)/2; j1++){
            bases_eval0(2*j1+1) = sin(2*(j1+1)*pi*x_c1);
            bases_eval0(2*j1+2) = cos(2*(j1+1)*pi*x_c1);
            bases_eval1(2*j1+1) = sin(2*(j1+1)*pi*x_c1_next);
            bases_eval1(2*j1+2) = cos(2*(j1+1)*pi*x_c1_next);
          }
          delta1 = bases_eval0*ww(r)*(A(r,t1,i)==target_policy)-bases_eval1;
          bases_eval0.ones(d_w);
          bases_eval1.ones(d_w);
          for(j1=0; j1<(d_w-1)/2; j1++){
            bases_eval0(2*j1+1) = sin(2*(j1+1)*pi*x_c2);
            bases_eval0(2*j1+2) = cos(2*(j1+1)*pi*x_c2);
            bases_eval1(2*j1+1) = sin(2*(j1+1)*pi*x_c2_next);
            bases_eval1(2*j1+2) = cos(2*(j1+1)*pi*x_c2_next);
          }
          delta2 = bases_eval0*ww(r)*(A(r,t2,i)==target_policy)-bases_eval1;

          kappa = exp(-(x_c1_next-x_c2_next)*(x_c1_next-x_c2_next)/2);
          for(j1=0; j1<d_w; j1++)
            for(j2=0; j2<d_w; j2++){
              G1(j1,j2) = delta1(j1) * delta2(j2);
            }
          G = G + kappa * G1;
        }
    }

    eig_sym(eigval, eigvec, G);
    beta.row(r) = trans(eigvec.col(0));
  }
    
  //evaluate on the second subsample (idx1)
  if(idx(0)==1){idx1 = idx+14;}
  if(idx(0)==16){idx1 = idx-16;}

  for(r=0; r < R; r++){

    beta1 = - beta.row(r);
    Stilde0 = X.row(r);

    for(i1=0; i1<N/2; i1++){

      i = idx1(i1);
      Stilde = Stilde0.col(i);
      for(t=1; t<M; t++){
        bases_eval0.ones(d_w);
        for(j1=0; j1<(d_w-1)/2; j1++){
          bases_eval0(2*j1+1) = sin(2*(j1+1)*pi*Stilde(t));
          bases_eval0(2*j1+2) = cos(2*(j1+1)*pi*Stilde(t));
        }
        w(r,t-1,i1) = sum(bases_eval0%trans(beta1));
      }
    }
  }

  return w;
}







