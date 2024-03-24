########### smoothing fns
# Epan <- function(z){
#   return( 3/4 * (1-z^2) * (abs(z)<1) )
# } 
# 
# LC <- function(x, y, eval, h, N, d){
#   
#   est <- c()
#   
#   if(d==1){
#     
#     EV <- length(eval)
#     
#     for(i in 1:EV){
#       
#       Kvec <- Epan((x - eval[i])/h)/h
#       est <- c(est, sum(Kvec*y)/sum(Kvec))
#       
#     }
#     
# 
#   }else{
#     
#     EV <- nrow(eval)
# 
#     for(i in 1:EV){
#       
#       Kvec <- (Epan((x[,1] - eval[i,1])/h)/h
#                 *Epan((x[,2] - eval[i,2])/h)/h)
#       est <- c(est, sum(Kvec*y)/sum(Kvec))
#       
#     }
#     
#   }
#   
#   return(est)
# }
# 
# smo_mat <- function(R,M,lat,lon,y,h1,h2){
# 
#   ## temporal smo
#   for(i in 1:R){
#     y[i,] <- LC((1:M)/M, y[i,], (1:M)/M, h2, M, 1) 
#   }
# 
#   ## spatial smo
#   for(i in 1:M){
#     y[,i] <- LC(cbind(lat,lon), y[,i], cbind(lat,lon), h1, R, 2)
#   }
# 
#   return(y)
#   
# }

########### statistics of interest
OLS_ATE_est <- function(y,design,spill,A,X,nb){
  
  A_bar <- array(0,dim=c(R,M,N))
  for(i in 1:R)
    for(j in 1:M){
      nb1 <- length(nb[i][[1]])
      if(nb1>1){A_bar[i,j,] <- colMeans(as.matrix(A[nb[i][[1]],j,]))}
      if(nb1==1){A_bar[i,j,] <- A[nb[i][[1]],j,]}
    }
  
  if(spill==0|design==0){
    
    ## outcome coefficient estimation
    coeff <- sapply(1:R,function(i){sapply(1:M,function(j){
      xc <- cbind(rep(1,N), X[i,j,], A[i,j,])
      return((ginv(t(xc)%*%xc)%*%t(xc)%*%matrix(y[i,j,],N,1))[,1])
    })})
    coeff <- array(coeff, dim=c(3,M,R))
    alpha_hat <- t(matrix(coeff[1,,],M,R))
    beta_hat <- t(matrix(coeff[2,,],M,R))
    gamma_hat <- t(matrix(coeff[3,,],M,R))
    theta_hat <- t(matrix(0,M,R))
    ## state coefficient estimation
    coeff <- sapply(1:R,function(i){sapply(1:(M-1),function(j){
      xc <- cbind(rep(1,N), X[i,j,], A[i,j,])
      return((ginv(t(xc)%*%xc)%*%t(xc)%*%matrix(X[i,j+1,],N,1))[,1])
    })})
    coeff <- array(coeff, dim=c(3,M-1,R))
    alpha_X_hat <- t(matrix(coeff[1,,],M-1,R))
    beta_X_hat <- t(matrix(coeff[2,,],M-1,R))
    gamma_X_hat <- t(matrix(coeff[3,,],M-1,R))
    theta_X_hat <- t(matrix(0,M-1,R))
                 
  }else{
    
    ## outcome coefficient estimation
    coeff <- sapply(1:R,function(i){sapply(1:M,function(j){
      xc <- cbind(rep(1,N), X[i,j,], A[i,j,],A_bar[i,j,])
      return((ginv(t(xc)%*%xc)%*%t(xc)%*%matrix(y[i,j,],N,1))[,1])
    })})
    coeff <- array(coeff, dim=c(4,M,R))
    alpha_hat <- t(matrix(coeff[1,,],M,R))
    beta_hat <- t(matrix(coeff[2,,],M,R))
    gamma_hat <- t(matrix(coeff[3,,],M,R))
    theta_hat <- t(matrix(coeff[4,,],M,R))
    ## state coefficient estimation
    coeff <- sapply(1:R,function(i){sapply(1:(M-1),function(j){
      xc <- cbind(rep(1,N), X[i,j,], A[i,j,],A_bar[i,j,])
      return((ginv(t(xc)%*%xc)%*%t(xc)%*%matrix(X[i,j+1,],N,1))[,1])
    })})
    coeff <- array(coeff, dim=c(4,M-1,R))
    alpha_X_hat <- t(matrix(coeff[1,,],M-1,R))
    beta_X_hat <- t(matrix(coeff[2,,],M-1,R))
    gamma_X_hat <- t(matrix(coeff[3,,],M-1,R))
    theta_X_hat <- t(matrix(coeff[4,,],M-1,R))
    
  }
  
  ## compute ATE
  {
    c_hat <- c()
    for(k in 1:(M-1)){
      for(tau in (k+1):M){
        if(tau==k+1){b <- beta_hat[,tau]}
        if(tau>k+1){
          p <- 1
          for(j in (k+1):(tau-1)){p <- p*beta_X_hat[,j]}
          b <- b + beta_hat[,tau] * p
        }
      }
      c_hat <- cbind(c_hat,b)
    }
    ATE_hat <- sum(gamma_hat+theta_hat) + sum(c_hat * (gamma_X_hat+theta_X_hat))
  }
  
  res <- list(ATE_hat, alpha_X_hat, beta_X_hat, gamma_X_hat, theta_X_hat,
              alpha_hat, beta_hat, gamma_hat, theta_hat)
  names(res) <- c('ATE_hat','alpha_X_hat', 'beta_X_hat', 'gamma_X_hat', 'theta_X_hat',
                  'alpha_hat', 'beta_hat', 'gamma_hat', 'theta_hat')
  
  return(res)
}

# smo_ATE_est <- function(R,M,lat,lon,ols_res,h1,h2){
#   
#   ## smoothing
#   res1 <- ols_res
#   for(i in 2:5){
#     res1[[i]] <- smo_mat(R,M-1,lat,lon,res1[[i]],h1,h2)
#   }
#   for(i in 6:9){
#     res1[[i]] <- smo_mat(R,M,lat,lon,res1[[i]],h1,h2)
#   }
#   ## compute ATE
#   attach(res1, warn.conflicts = FALSE)
#   {
#     c_smo <- c()
#     for(k in 1:(M-1)){
#       for(tau in (k+1):M){
#         if(tau==k+1){b <- beta_hat[,tau]}
#         if(tau>k+1){
#           p <- 1
#           for(j in (k+1):(tau-1)){p <- p*beta_X_hat[,j]}
#           b <- b + beta_hat[,tau] * p
#         }
#       }
#       c_smo <- cbind(c_smo,b)
#     }
#     ATE_smo <- sum(gamma_hat+theta_hat) + sum(c_smo * (gamma_X_hat+theta_X_hat))
#   }
#   res <- list(ATE_smo, alpha_X_hat, beta_X_hat, gamma_X_hat, theta_X_hat,
#               alpha_hat, beta_hat, gamma_hat, theta_hat)
#   names(res) <- c('ATE_smo','alpha_X_smo', 'beta_X_smo', 'gamma_X_smo', 'theta_X_smo',
#                   'alpha_smo', 'beta_smo', 'gamma_smo', 'theta_smo')
#   detach(res1)
#   
#   return(res)
# }

stat_est <- function(y,design,spill,A,nb,X,e,lat,lon){

  gap_lat <- sort(unique(lat))[2]-sort(unique(lat))[1]
  gap_lon <- sort(unique(lon))[2]-sort(unique(lon))[1]
  h1 <- 1.5*max(gap_lat, gap_lon)*(M==1)+0.8*min(gap_lat,gap_lon)*(M>1)
  h2 <- 1.5/M
  
  ### estimate
  ols_res <- OLS_ATE_est(y,design,spill,A,X,nb)
  # smo_res <- smo_ATE_est(R,M,lat,lon,ols_res,h1,h2)
  
  ### bootstrap
  if(stren==0){
    ATE_b <- c()
    for(i_b in 1:Nb){
      set.seed(i_b)
      idx_b <- sample(1:N, N, replace = TRUE)
      X_b <- X[,,idx_b]; y_b <- y[,,idx_b]; A_b <- A[,,idx_b]
      ols_res_b <- OLS_ATE_est(y_b,design,spill,A_b,X_b,nb)
      ATE_b <- c(ATE_b, ols_res_b$ATE_hat)
      # ATE_smo_b <- c(ATE_b, smo_ATE_est(R,M,lat,lon,ols_res_b,h1,h2)$ATE_smo)
    }
    var_hat <- var(ATE_b)
  }else{
    var_hat <- 0
  }
  # var_smo <- var(ATE_smo_b)
  # res <- list(c(ols_res$ATE_hat,smo_res$ATE_smo),c(var_hat,var_smo)); names(res) <- c('ATE_est','var_est')
  res <- list(ols_res$ATE_hat, var_hat); names(res) <- c('ATE_est','var_est')
  
  return(res)
}
    
model_aggre <- function(R,M,design,spill){

  ### treatments generation
  {
    A <- array(0,dim=c(R,M,N))
    if(design==0)
    {
      for(i in 1:N){a <- rbinom(1,1,0.5); A[,,i] <- a}
    }
    if(design==1)
    {
      for(i in 1:N){
        a <- rbinom(R,1,0.5); A[,,i] <- a
        # if(M>1){A[,seq(2,M,2),i] <- 1-a}
      }
    }
    if(design==2)
    {
      for(i in 1:N)
        for(c in 1:N_cluster)
        {
          a <- rbinom(1,1,0.5)
          A[((c-1)*9+1):(c*9),,i] <- a
        }
    }
    
  }
  
  ### responses generation
  X <- state_gene(design, spill, A, nb, X, e_X, alpha_X, beta_X, gamma_X, theta_X)
  y <- outcome_gene(design, spill, A, nb, X, e, alpha, beta, gamma, theta)
  
  ### estimate
  res <- stat_est(y, design, spill, A, nb, X, e, lat,lon)
  
  return(res)
  
}


