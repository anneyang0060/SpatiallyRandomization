########### statistics of interest
stat_est <- function(y,design,spill,X,A){
  
  A_bar <- array(0,dim=c(R,N))
  for(i in 1:R){
    nb1 <- length(nb[[i]])
    if(nb1>1){A_bar[i,] <- colMeans(A[nb[[i]],])}
    if(nb1==1){A_bar[i,] <- A[nb[[i]],]}
  }
  
  #### OLS ATE est
  if(spill==0|design==0)
  {
    ## region-wise and timepiece-wise estimation
    coeff <- sapply(1:R,function(i){sapply(1:M,function(j){
      xc <- cbind(rep(1,N), X[i,], A[i,])
      return((ginv(t(xc)%*%xc)%*%t(xc)%*%matrix(y[i,],N,1))[,1])
    })})
    coeff <- matrix(coeff, 3,R)
    gamma_hat <- coeff[3,]
    ATE_hat <- sum(gamma_hat)
  }
  if(design==1&spill==1)
  {
    coeff <- sapply(1:R,function(i){sapply(1:M,function(j){
      xc <- cbind(rep(1,N), X[i,], A[i,],A_bar[i,])
      return((ginv(t(xc)%*%xc)%*%t(xc)%*%matrix(y[i,],N,1))[,1])
    })})
    coeff <- matrix(coeff,4,R)
    gamma_hat <- coeff[3,]; theta_hat <- coeff[4,]
    ATE_hat <- sum(gamma_hat+theta_hat)
  }
  if(design==2&spill==1)
  {
    coeff1 <- sapply(which(Nc_group==1),function(i){sapply(1:M,function(j){
      xc <- cbind(rep(1,N), X[i,], A[i,])
      return((ginv(t(xc)%*%xc)%*%t(xc)%*%matrix(y[i,],N,1))[,1])
    })})
    coeff1 <- matrix(coeff1, nrow=3)
    
    coeff2 <- sapply(which(Nc_group>1),function(i){sapply(1:M,function(j){
      xc <- cbind(rep(1,N), X[i,], A[i,],A_bar[i,])
      return((ginv(t(xc)%*%xc)%*%t(xc)%*%matrix(y[i,],N,1))[,1])
    })})
    coeff2 <- matrix(coeff2,nrow=4)
    
    coeff <- matrix(0,4,R)
    coeff[1:3, which(Nc_group==1)] <- coeff1
    coeff[, which(Nc_group>1)] <- coeff2
    gamma_hat <- coeff[3,]; theta_hat <- coeff[4,]
    ATE_hat <- sum(gamma_hat+theta_hat)
    
  }
  
  # #### OLS var est (sandwich)
  if(gamma_strength == 0){
    ### noise est
    e_hat <- matrix(0, R, N)
    for(i in 1:R){
        if(spill == 0|design==0){
          e_hat[i,] <- y[i,] - (coeff[1,i] + X[i,] * coeff[2,i] + coeff[3,i] * A[i,])
        }else{
          e_hat[i,] <- y[i,] - (coeff[1,i] + X[i,] * coeff[2,i] + coeff[3,i] * A[i,]
                                          + coeff[4,i] * A_bar[i,])
        }
      }
    if(design==0|spill==0){A_bar <- matrix(0,R,N)}
    Sigma_hat0 <- 1/(N-3-spill) * (e_hat %*% t(e_hat))
    V_hat0 <- sandwich(X, A, A_bar, Sigma_hat0, Nc_group, R, N, spill, design)
    var_hat <- sum(V_hat0)
  }

  if(gamma_strength == 0){
    res <- list(ATE_hat,var_hat); names(res) <- c('ATE_est','var_est')
  }else{
      res <- list(ATE_hat); names(res) <- c('ATE_est')
  }

  return(res)
  # return(ATE_hat)
}

# stat_est <- function(y,design,spill,X,A){
#   ### estimate
#   ATE_hat <- ATE_est(y,design,spill,X,A)
#   # smo_res <- smo_ATE_est(R,M,lat,lon,ols_res,h1,h2)
#   
#   ### bootstrap
#   if(stren==0)
#   {
#     ATE_b <- c()
#     for(i_b in 1:Nb){
#       set.seed(i_b)
#       idx_b <- sample(1:N, N, replace = TRUE)
#       X_b <- X[,idx_b]; y_b <- y[,idx_b]; A_b <- A[,idx_b]
#       ols_res_b <- ATE_est(y_b,design,spill,X_b,A_b)
#       ATE_b <- c(ATE_b, ols_res_b)
#       # ATE_smo_b <- c(ATE_b, smo_ATE_est(R,M,lat,lon,ols_res_b,h1,h2)$ATE_smo)
#     }
#     var_hat <- var(ATE_b)
#   }
#   else
#   {
#     var_hat <- 0
#   }
#   # var_smo <- var(ATE_smo_b)
#   # res <- list(c(ols_res$ATE_hat,smo_res$ATE_smo),c(var_hat,var_smo)); names(res) <- c('ATE_est','var_est')
#   res <- list(ATE_hat, var_hat); names(res) <- c('ATE_est','var_est')
#   
#   return(res)
#   
# }

model_aggre <- function(R,design,spill){
  
  ### treatments generation
  {
    set.seed(sd)
    A <- matrix(0,R,N)
    for(i in 1:N)
    {
      a <- rbinom(R,1,0.5);
      if(design==0){A[,i] <- a[1]}
      if(design==1){A[,i] <- a}
      if(design==2)
      {
        for(c in 1:N_cluster)
        {
          A[((c-1)*9+1):(c*9),i] <- a[c]
        }
      }
    }
  }  
  
  ### responses generation
  y <- outcome_gene(design, spill, A)
  
  ### estimate
  res <- stat_est(y, design, spill,X, A)
  
  return(res)
  
}











