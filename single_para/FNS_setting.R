########### fns to generate coefficients
coeff_fn0 <- function(x,a){
  n <- length(a)
  res <- c()
  for(i in 1:length(x)){
    res <- c(res, a[1] + sum(a[seq(2,n,2)]*cos(pi*seq(2,n,2)*x[i])) +
               sum(a[seq(3,n,2)]*sin(pi*seq(3,n,2)*x[i])))
  }
  return(res)
}

coeff_fn <- function(k,lat,lon){
  cs <- matrix(runif(3*k,0,1),3,k)
  y <- coeff_fn0(lat,cs[1,]) + coeff_fn0(lon, cs[2,])
  return(y)
  # return(cs)
}

########### fn to generate regions of given shape and number
symmetry <- function(coor0,A,B,C){
  x0 <- coor0[1]
  y0 <- coor0[2]
  x <- x0 - 2*A * (A*x0+B*y0+C) / (A^2+B^2)
  y <- y0 - 2*B * (A*x0+B*y0+C) / (A^2+B^2)
  coor1 <- c(x,y)
  return(coor1)
}

cluster_translation <- function(coors0, dx, dy){
  coors <- coors0
  coors[, 1] <- coors0[,1] + dx
  coors[, 2] <- coors0[,2] + dy
  return(coors)
}

region_gene <- function(side, R){
  
  coor <- matrix(0, 144, 2)
  cluster <- rep(1:16, each=9)
  #coordinates
  if(side==3)
  {
    # generate the coordinates of regions in the 1st cluster
    {
      coor[1, ] <- c(1, sqrt(3)/3)
      coor[2, ] <- c(2, 2*sqrt(3)/3)
      coor[c(3,5),] <- rbind(coor[1, ]+c(2,0), coor[1, ]+c(4,0))
      coor[4,] <- coor[2,]+c(2,0)
      coor[6:8,] <- coor[1:3,]+matrix(rep(c(1,sqrt(3)),each=3),3,2)
      coor[9,] <- coor[1,] + c(2,2*sqrt(3))
    }
    
    # region coordinates in the 2nd cluster (use symmetry)
    {
      A <- sqrt(3); B <- 1; C <- -6*sqrt(3)
      for(i in 10:18){coor[i,] <- symmetry(coor[i-9,],A,B,C)}
    }
    
    # generate other clusters according to cluster 1 and 2
    {
      # cluster 3, 5, 7
      coor[(2*9+1):(3*9),] <- cluster_translation(coor[1:9,], 6, 0)
      coor[(4*9+1):(5*9),] <- cluster_translation(coor[1:9,], 12, 0)
      coor[(6*9+1):(7*9),] <- cluster_translation(coor[1:9,], 18, 0)
      # cluster 4, 6, 8
      coor[(3*9+1):(4*9),] <- cluster_translation(coor[10:18,], 6, 0)
      coor[(5*9+1):(6*9),] <- cluster_translation(coor[10:18,], 12, 0)
      coor[(7*9+1):(8*9),] <- cluster_translation(coor[10:18,], 18, 0)
      # cluster 9 -- 18
      coor[73:144,] <- cluster_translation(coor[1:72,], 3, 3*sqrt(3))
    }
    
    if(R==36){coor <- coor[1:36, ]; cluster <- cluster[1:36]}
    if(R==81){coor <- coor[c(1:45,73:108), ]; cluster <- cluster[c(1:45,73:108)]}
  }
  if(side == 4)
  {
    coor[1:9, 1] <- rep(c(1,3,5),3)
    coor[1:9, 2] <- rep(c(1,3,5),each = 3)
    coor[10:18,] <- cluster_translation(coor[1:9,],6,0)
    coor[19:36,] <- cluster_translation(coor[1:18,],12,0)
    coor[37:72,] <- cluster_translation(coor[1:36,],0,6)
    coor[73:144,] <- cluster_translation(coor[1:72,],0,12)
    if(R==36){coor <- coor[c(1:18,37:54), ]; cluster <- cluster[c(1:18,37:54)]}
    if(R==81){coor <- coor[c(1:27,37:63,73:99), ]; cluster <- cluster[c(1:27,37:63,73:99)]}
  }
  if(side == 6)
  {
    # cluster 1
    coor[1,] <- c(2, sqrt(3))
    coor[2,] <- c(5, 2*sqrt(3))
    coor[3,] <- c(8, sqrt(3))
    coor[4:6,] <- cluster_translation(coor[1:3,],0,2*sqrt(3))
    coor[7:9,] <- cluster_translation(coor[1:3,],0,4*sqrt(3))
    # other clusters
    coor[10:18,] <- cluster_translation(coor[1:9,], 9, -sqrt(3))
    coor[19:36,] <- cluster_translation(coor[1:18,], 18, 0)
    coor[37:72,] <- cluster_translation(coor[1:36,], 0, 6*sqrt(3))
    coor[73:144,] <- cluster_translation(coor[1:72,], 0, 12*sqrt(3))
    if(R==36){coor <- coor[c(1:18,37:54), ]; cluster <- cluster[c(1:18,37:54)]}
    if(R==81){coor <- coor[c(1:27,37:63,73:99), ]; cluster <- cluster[c(1:27,37:63,73:99)]}
  }
  #neighbors
  nb <- list() 
  thre <- 2/3*sqrt(3)*(side==3) + 2*(side==4) + 2*sqrt(3)*(side==6)
  for(r in 1:R){
    dis <- sapply(1:R, function(i){sqrt((coor[i,1]-coor[r,1])^2+(coor[i,2]-coor[r,2])^2)})
    nb1 <- which(dis<=thre+0.1)
    nb1 <- nb1[nb1!=r]
    nb[[r]] <- nb1
  }
  # No. clusters of the group (the neighbors and itself)
  Nc_group <- c()
  for(i in 1:R){Nc_group[i] <- length(unique(cluster[c(i, nb[[i]])]))}
  
  mylist <- list(coor, nb, cluster, Nc_group)
  names(mylist) <- c('coor', 'nb', 'cluster', 'Nc_group')
  return(mylist)
}

########### noise generation under different convarianece structures
noise_gene <- function(mode, sd, R, M, rho, sigma, dis_mat){
  
  e <- matrix(0,R,N)
  
  # exchangeable
  if(mode==1){
    cov_eta <- matrix(rho,R,R); diag(cov_eta) <- 1
    set.seed(sd)
    e <- t(rmvn(N,rep(0,R),cov_eta))
  }
  
  # exponential
  if(mode==2){
    cov_eta <- (rho)^(dis_mat)
    set.seed(sd)
    e <- t(rmvn(N,rep(0,R),cov_eta))
  }
  
  # m dependent
  if(mode==3){
    set.seed(123)
    v1 <- 0.25*runif(R)+0.75; v1[ceiling(R*rho):R] <- 0
    cov_eta <- v1%*%t(v1); diag(cov_eta) <- 1
    set.seed(sd)
    e <- t(rmvn(N,rep(0,R),cov_eta))
  }
  
  # if(mode==3){
  #   m <- ceiling(rho*R)
  #   e0 <- matrix(rnorm(N*rep(R+m)),nrow=N)
  #   for(j in 1:R){
  #     e[j,1,] <- rowSums(e0[,j:(j+M)])
  #   }
  #   e <- e/sum(1:(m+1))
  # }


  return(e)
}

########### outcomes
outcome_gene <- function(design, spill, A){
  
  y <- matrix(0,R,N)
  A_bar <- matrix(0,R,N)
  
  for(i in 1:R){
      if(spill==1){
        
        nb1 <- length(nb[[i]])
        if(nb1>1){A_bar[i,] <- colMeans(as.matrix(A[nb[i][[1]],]))}
        if(nb1==1){A_bar[i,] <- A[nb[[i]],]}
        y[i,] <- alpha[i] + X[i,] * beta[i] + gamma[i] * A[i,] + theta[i] * A_bar[i,] + e[i,]

      }
      if(spill==0){
        y[i,] <- alpha[i] + X[i,] * beta[i] + gamma[i] * A[i,] + e[i,]
      }
    }
  
  return(y)
}

