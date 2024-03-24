scale_fn <- function(x){(x-min(x))/(max(x)-min(x))}

getQ <- function(X, A,  A_bar, y, idx){
  
  idx1 <- setdiff(1:N,idx)
  r1_hat <- array(0,dim = c(R,M,N/2)); q1_hat <- r1_hat
  r0_hat <- r1_hat; q0_hat <- r0_hat
  r_hat <- r1_hat; q_hat1 <- r_hat; q_hat0 <- q_hat1
  
  for(t in 1:M){
    
    x <- (as.vector(X[,t,idx]))
    loc <- (rep(lat,N/2)+rep(lon,N/2))/2
    a <- as.vector(A[,t,idx])
    a_bar <- as.vector(A_bar[,t,idx])
    xc <- cbind(x/2,loc, a, a_bar)
    yc <- as.vector(y[,t,idx])
    
    x <- (as.vector(X[,t,idx1]))
    a <- as.vector(A[,t,idx1])
    a_bar <- as.vector(A_bar[,t,idx1])
    eval1 <- cbind(x/2, loc, rep(1,N/2),rep(1,N/2))
    eval0 <- cbind(x/2, loc, rep(0,N/2),rep(0,N/2))
    
    h <- c(0.2,0.4, 0.6, 0.6)
    
    y_est1 <- LC(xc,yc,eval1,h)
    y_est1[is.na(y_est1)] <- mean(y_est1[!is.na(y_est1)])
    y_est0 <- LC(xc,yc,eval0,h)
    y_est0[is.na(y_est0)] <- mean(y_est0[!is.na(y_est0)])
    y_est <- LC(xc,yc,cbind(x/2, loc, a, a_bar),h)
    y_est[is.na(y_est)] <- mean(y_est[!is.na(y_est)])
    r1_hat[,t,] <- matrix(y_est1,R,N/2)
    r0_hat[,t,] <- matrix(y_est0,R,N/2)
    r_hat[,t,] <- matrix(y_est,R,N/2)
    
  }
  
  q1_hat[,M,] <- r1_hat[,M,]; q0_hat[,M,] <- r0_hat[,M,]
  q_hat1[,M,] <- r_hat[,M,]; q_hat0[,M,] <- r_hat[,M,]
  for (t in (M-1):1){
    q1_hat[,t,] <- apply(r1_hat[,t:M,],c(1,3),sum)
    q_hat1[,t,] <- q1_hat[,t,] - r1_hat[,t,] + r_hat[,t,]
    q0_hat[,t,] <- apply(r0_hat[,t:M,],c(1,3),sum)
    q_hat0[,t,] <- q0_hat[,t,] - r0_hat[,t,] + r_hat[,t,]
  }

  q_res <- list(q1_hat,q0_hat,q_hat1,q_hat0)
  names(v_res) <- c('q1_hat','q0_hat','q_hat1','q_hat0')
  
  return(v_res)
  
}

########### statistics of interest
ATE_est_fn <- function(X,y,design,spill,A,nb, mu1_, mu0_, bp){
  
  A_bar <- array(0,dim=c(R,M,N))
  X_bar <- array(0,dim=c(R,M,N))
  for(i in 1:R)
    for(t in 1:M){
      
      if(length(nb[[i]])>1){
        A_bar[i,t,] <- colMeans(A[nb[[i]],t,])
        X_bar[i,t,] <- colMeans(X[nb[[i]],t,])
      }else{
        A_bar[i,t,] <- A[nb[[i]],t,]
        X_bar[i,t,] <- X[nb[[i]],t,]
      }
    }
  
  ATE_est <- 0
  mu1 <- array(0, dim = c(R,M,N)); mu0 <- mu1
  
  for(fold in 1:2){
    
    if(fold==1){idx <- 1:(N/2)}else{idx <- (N/2+1):N}# for omega and q estimation
    idx1 <- setdiff(1:N,idx)
    
    if(bp==FALSE)
    {
      # w1 <- array(1, dim = c(R,M,N/2)); w0 <- w1
      if(design==0)
      {ww <- array(1/0.5,dim=c(R,M,N/2))}
      if(design==1)
      {
        ww <- sapply(1:R, function(k){max(0.5^(length(nb[[k]])+1),0.1)})
        ww <- array(rep(1/ww,M*N/2),dim=c(R,M,N/2))
      }
      if(design==2)
      {
        ww <- sapply(1:R, function(k){max(0.5^(Nc_group[k]),0.05)})
        ww <- array(rep(1/ww,M*N/2),dim=c(R,M,N/2))
      }
      w0 <- getWeight(X,  A,  idx, ww, target_policy = 0, d_w)
      w1 <- getWeight(X,  A,  idx, ww, target_policy = 1, d_w)
      w0 <- scale_fn(abs(w0)); w1 <- scale_fn(abs(w1))
      # w0 <- w0-min(w0); w1 <- w1 - min(w1)
      mu1[,,idx1] <- w1*(A[,,idx1]==1)*ww
      mu0[,,idx1] <- w0*(A[,,idx1]==0)*ww
    }
    else
    {
      mu1 <- mu1_; mu0 <- mu0_
    }
    q_res <- getQ(X,  A, A_bar, y, idx)
    attach(q_res)
    rm(v_res)
    V1 <- sum(mu1[,,idx1]*(y[,,idx1]-q_hat1))+sum(q1_hat)
    V0 <- sum(mu0[,,idx1]*(y[,,idx1]-q_hat0))+sum(q0_hat)
    detach(q_res)
    ATE_est <- ATE_est + (V1-V0)/N
  }
  
  res <- list(ATE_est, mu1, mu0)
  names(res) <- c('ATE_est','mu1','mu0')
  return(res)
}

stat_est <- function(X,y,design,spill,A,nb, mu1_, mu0_){
  
  # X <- (X-min(X)) / (max(X)-min(X))
  #### DR ATE est
  ATE_est_res <- ATE_est_fn(X,y,design,spill,A,nb, mu1_, mu0_, FALSE)
  ATE_est <- ATE_est_res$ATE_est
  mu1 <- ATE_est_res$mu1
  mu0 <- ATE_est_res$mu0
  
  #### var est
  if(stren==0){
    ATE_est_b <- c()
    for(i_b in 1:N_b){
      ATE_est_b1 <- c()
      set.seed(i_b)
      idx_b <- sample(1:N, N, replace = TRUE)
      X_b <- X[,,idx_b]; y_b <- y[,,idx_b]; A_b <- A[,,idx_b]
      ATE_est_b <- c(ATE_est_b, ATE_est_fn(X_b,y_b,design,spill,A_b,nb, mu1[,,idx_b], mu0[,,idx_b], TRUE)$ATE_est)
      # ATE_est_b <- c(ATE_est_b, ATE_est_fn(X_b,y_b,design,spill,A_b,nb, 0, 0, TRUE)$ATE_est)
    }
    var_est <- var(ATE_est_b[which(!is.na(ATE_est_b))])
  }else{
    var_est <- 0
  }
  
  res <- list(ATE_est,var_est,mu1,mu0)
  names(res) <- c('ATE_est','var_est','mu1','mu0')
  
  return(res)
}

model_aggre <- function(R,design,spill,sd, nb, mu1_, mu0_){
  
  set.seed(sd)
  
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
  X <- state_gene(A)
  y <- outcome_gene(X, spill, A)
  
  ### estimate
  res <- stat_est(X,y,design,spill,A,nb, mu1_, mu0_)
  
  return(res)
  
}











