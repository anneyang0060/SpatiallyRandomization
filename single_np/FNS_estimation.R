########### statistics of interest
ATE_est_fn <- function(X,y,design,spill,A,nb){
  
  A_bar <- matrix(0,R,N)
  X_bar <- matrix(0,R,N)
  for(i in 1:R){
    
    if(length(nb[[i]])>1){
      A_bar[i,] <- colMeans(A[nb[[i]],])
      X_bar[i,] <- colMeans(X[nb[[i]],])
    }else{
      A_bar[i,] <- (A[nb[[i]],])
      X_bar[i,] <- (X[nb[[i]],])
    }
    
  }
  
  ATE_est <- 0
  
  for(fold in 1:2){
    
    if(fold==1){idx <- 1:(N/2)}else{idx <- (N/2+1):N}
    
    idx1 <- setdiff(1:N,idx)
    
    xc <- cbind((as.vector(X[,idx])+as.vector(X_bar[,idx])-6)/4, 
                (rep(lat,N/2)+rep(lon,N/2))/2, as.vector(A[,idx]))
    
    eval0 <- cbind((as.vector(X[,idx1])+as.vector(X_bar[,idx1])-6)/4, 
                   (rep(lat,N/2)+rep(lon,N/2))/2, rep(0,R*N/2))
    eval1 <- cbind((as.vector(X[,idx1])+as.vector(X_bar[,idx1])-6)/4, 
                   (rep(lat,N/2)+rep(lon,N/2))/2, rep(1,R*N/2))
    
    h <- c(0.3, 1.2/2*(sort(unique(as.vector(lat)))[2]+sort(unique(as.vector(lon)))[2]),0.3)
    
    if(spill==1&design >= 1){
      
      xc <- cbind(xc, as.vector(A_bar[,idx]))
      
      eval0 <- cbind(eval0, rep(0,R*N/2))
      eval1 <- cbind(eval1, rep(1,R*N/2))
      
      # h <- c(h, 1.2*sort(unique(as.vector(A_bar)))[2])
      h <- c(h, 0.8)
    }
    yc <- as.vector(y[,idx])
    
    y_est1 <- LC(xc,yc,eval1,h); y_est1[is.na(y_est1)] <- mean(y_est1[!is.na(y_est1)])
    y_est0 <- LC(xc,yc,eval0,h); y_est0[is.na(y_est0)] <- mean(y_est0[!is.na(y_est0)])
    
    if(design==0|spill==0){
      
      ww <- rep(1/0.5, R*N/2)
      # ww <- rep(0,R*N)
      r1 <- (as.vector(y[,idx1]) - y_est1)*(as.vector(A[,idx1])==1)*ww + y_est1
      r0 <- (as.vector(y[,idx1]) - y_est0)*(as.vector(A[,idx1])==0)*ww + y_est0
      
    }
    if(design==1&spill==1){
      
      ww <- sapply(1:R, function(k){max(0.5^length(nb[[k]]), 0.025)})
      # ww <- sapply(1:R, function(k){0.5^length(nb[[k]])})
      ww <- rep(1/ww, N/2)
      # ww <- rep(0,R*N)
      r1 <- (as.vector(y[,idx1]) - y_est1)*(as.vector(A[,idx1]+A_bar[,idx1])==2)*ww + y_est1
      r0 <- (as.vector(y[,idx1]) - y_est0)*(as.vector(A[,idx1]+A_bar[,idx1])==0)*ww + y_est0
      
    }
    if(design==2&spill==1){
      
      ww <- sapply(1:R, function(k){max(0.5^length(Nc_group[k]), 0.025)})
      # ww <- sapply(1:R, function(k){0.5^length(Nc_group[k])})
      ww <- rep(1/ww, N/2)
      # ww <- rep(0,R*N)
      r1 <- (as.vector(y[,idx1]) - y_est1)*(as.vector(A[,idx1]+A_bar[,idx1])==2)*ww + y_est1
      r0 <- (as.vector(y[,idx1]) - y_est0)*(as.vector(A[,idx1]+A_bar[,idx1])==0)*ww + y_est0
      
    }
    
    
    ATE_est <- ATE_est + sum(r1 - r0)
  }
  
  ATE_est <- ATE_est/N
  
  return(ATE_est)
}

stat_est <- function(X,y,design,spill,A,nb){
  
  #### DR ATE est
  ATE_est <- ATE_est_fn(X,y,design,spill,A,nb)
  
  #### var est
  if(N_b>0){
    ATE_est_b <- c()
    for(i_b in 1:N_b){
      set.seed(i_b)
      idx_b <- sample(1:N, N, replace = TRUE)
      X_b <- X[,idx_b]; y_b <- y[,idx_b]; A_b <- A[,idx_b]
      ATE_est_b <- c(ATE_est_b, ATE_est_fn(X_b,y_b,design,spill,A_b,nb))
    }
    var_est <- var(ATE_est_b[which(!is.na(ATE_est_b))])
    res <- list(ATE_est,var_est); names(res) <- c('ATE_est','var_est')
  }else{
    res <- list(ATE_est); names(res) <- c('ATE_est')
  }
  
  
  return(res)
}

model_aggre <- function(X,R,design,spill,sd, nb){
  
  set.seed(sd)
  ### treatments generation
  {
    A <- matrix(0,R,N)
    if(design==0)
    {
      for(i in 1:N){a <- rbinom(1,1,0.5); A[,i] <- a}
    }
    if(design==1)
    {
      for(i in 1:N){
        a <- rbinom(R,1,0.5); A[,i] <- a
      }
    }
    if(design==2)
    {
      for(i in 1:N)
        for(c in 1:N_cluster)
        {
          a <- rbinom(1,1,0.5)
          A[((c-1)*9+1):(c*9),i] <- a
        }
    }
  }  
  
  ### responses generation
  y <- outcome_gene(spill, A)
  
  ### estimate
  res <- stat_est(X,y,design,spill,A,nb)
  
  return(res)
  
}











