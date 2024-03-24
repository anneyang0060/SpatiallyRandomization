getWeight1 <- function(X,X_bar, A, A_bar, idx, target_policy, d_w){
  
  # estimate on the first subsample (idx)
  beta <- matrix(0,R,d_w)
  
  for(r in 1:R){
    
    G <- 0
    
    for(i in idx){
      
      # Stilde <- scale_fn(X[r,,i]+X_bar[r,,i])
      Stilde <- (X[r,,i]+X_bar[r,,i])/2
      
      for(t1 in 1:(M-1))
        for(t2 in 1:(M-1)){
          
          x_c1 <- Stilde[t1]
          x_c2 <- Stilde[t2]
          x_c1_next <- Stilde[t1+1]
          x_c2_next <- Stilde[t2+1]
          
          bases_eval0 <- 1; bases_eval1 <- 1 
          for(i1 in 1:((d_w-1)/2)){
            bases_eval0 <- c(bases_eval0,c(sin(2*i1*pi*x_c1),cos(2*i1*pi*x_c1)))
            bases_eval1 <- c(bases_eval1,c(sin(2*i1*pi*x_c1_next),cos(2*i1*pi*x_c1_next)))
          }          
          delta1 <- bases_eval0*(A[r,t1,i]==target_policy & A_bar[r,t1,i]==target_policy)/max(0.5^(Nc_group[r]),0.05)-bases_eval1
          bases_eval0 <- 1; bases_eval1 <- 1 
          for(i1 in 1:((d_w-1)/2)){
            bases_eval0 <- c(bases_eval0,c(sin(2*i1*pi*x_c2),cos(2*i1*pi*x_c2)))
            bases_eval1 <- c(bases_eval1,c(sin(2*i1*pi*x_c2_next),cos(2*i1*pi*x_c2_next)))
          }          
          delta2 <- bases_eval0*(A[r,t2,i]==target_policy & A_bar[r,t2,i]==target_policy)/max(0.5^(Nc_group[r]),0.05)-bases_eval1
          kappa <- exp(-(x_c1_next-x_c2_next)^2/2)
          
          G <- G + matrix(delta1,ncol=1) %*% matrix(delta2,nrow=1) * kappa
        }
    }
    beta[r,] <- eigen(G)$vector[,d_w]
  }
  
  # evaluate on the second subsample (idx1)
  idx1 <- setdiff(1:N,idx)
  w <- array(0, dim = c(R,M,N/2))

  for(r in 1:R){

    beta1 <- beta[r,]

    for(i in idx1){

      Stilde <- (X[r,,i]+X_bar[r,,i])/2
      for(t in 2:M){
        bases_eval <- 1
        for(i1 in 1:((d_w-1)/2)){
          bases_eval <- c(bases_eval,c(sin(2*i1*pi*Stilde[t]),cos(2*i1*pi*Stilde[t])))
        }
        w[r,t-1,i-15] <- sum(bases_eval*beta1)
      }
    }
  }
  return(w)
}

