########### smoothing fns
Epan <- function(z){
  return( 3/4 * (1-z^2) * (abs(z)<1) )
} 

LC <- function(x, y, eval, h, N, d){
  
  est <- c()
  
  if(d==1){
    
    EV <- length(eval)
    
    for(i in 1:EV){
      
      Kvec <- Epan((x - eval[i])/h)/h
      est <- c(est, sum(Kvec*y)/sum(Kvec))
      
    }
    
    
  }else{
    
    EV <- nrow(eval)
    
    for(i in 1:EV){
      
      Kvec <- (Epan((x[,1] - eval[i,1])/h)/h
               *Epan((x[,2] - eval[i,2])/h)/h)
      est <- c(est, sum(Kvec*y)/sum(Kvec))
      
    }
    
  }
  
  return(est)
}

smo_mat <- function(y,h){
  
  ## spatial smo
  for(i in 1:ncol(y)){
    y[,i] <- LC(cbind(lat,lon), y[,i], cbind(lat,lon), h, R, 2)
  }
  

  return(y)
  
}

reformat <- function(data,idx){
  
  TI <- 144/M
  if(idx<7)
  {
    data1 <- data
    data1[which(data1$time==0),idx] <- 0
    data1[which(data1$time>0),idx] <- data[which(data$time<143),idx]
    data[,idx] <- data[,idx] - data1[,idx]
    data[which(data$time==143),idx] <- 0
  }
  
  y <- array(0, dim = c(R, M, N))
  for(i in 1:N)
    for(r in 1:R){
      z <- data[which(data$region==(r-1)&data$date==(stat_date+i-1)),idx]
      for(j in 1:M){
        y[r,j,i] <- sum(z[((j-1)*TI+1):(j*TI)])
      }
    }
  
  if(idx==7){
    y<-round(y/TI)
  }
  
  return(y)
}

########### statistics of interest
OLS_ATE_est <- function(y,design,spill,A,X,nb){
  
  A_bar <- array(0,dim=c(R,M,N))
  for(i in 1:R)
    for(j in 1:M){
      nb1 <- length(nb[i][[1]])
      if(nb1>1){A_bar[i,j,] <- colSums(as.matrix(A[nb[i][[1]],j,]))}
      if(nb1==1){A_bar[i,j,] <- A[nb[i][[1]],j,]}
      if(j>1){A_bar[i,j,] <- (A[i,j-1,] + A_bar[i,j,])/(nb1+1)}else{A_bar[i,j,] <- A_bar[i,j,]/nb1}
    }
  # A_bar <- 1*(A_bar>0)
  
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
      xc <- cbind(rep(1,N), X[i,j,], A[i,j+1,])
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
      xc <- cbind(rep(1,N), X[i,j,], A[i,j,], A_bar[i,j,])
      return((ginv(t(xc)%*%xc)%*%t(xc)%*%matrix(y[i,j,],N,1))[,1])
    })})
    coeff <- array(coeff, dim=c(4,M,R))
    alpha_hat <- t(matrix(coeff[1,,],M,R))
    beta_hat <- t(matrix(coeff[2,,],M,R))
    gamma_hat <- t(matrix(coeff[3,,],M,R))
    theta_hat <- t(matrix(coeff[4,,],M,R))
    
    ## state coefficient estimation
    coeff <- sapply(1:R,function(i){sapply(1:(M-1),function(j){
      xc <- cbind(rep(1,N), X[i,j,], A[i,j+1,], A_bar[i,j+1,])
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
    DE_hat <- sum(gamma_hat[Rs,1:M0]+theta_hat[Rs,1:M0])
    IE_hat <- sum(c_hat[Rs,1:(M0-1)] * (gamma_X_hat[Rs,1:(M0-1)]+theta_X_hat[Rs,1:(M0-1)]))
    ATE_hat <- DE_hat+IE_hat
  }
  

  return(c(ATE_hat, DE_hat, IE_hat))
 
  return(res)
}





