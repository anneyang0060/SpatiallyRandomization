############# table 5
setwd('~/spatio blessing/res/multi_para')
Nsim <- 500
N_struct <- 3

mse <- function(spill, design, side, rho){
  
  filename <- paste0('side',side,'_rho',rho,'_eX0.3.Rdata')
  
  load(filename)
  
  stat <- matrix(0, N_struct, Nsim)
  
  for(j in 1:N_struct)
    
    for(i in 1:Nsim){
      
      res1 <- res_c[[i]][[1]][[design]][[j]][[2]]
      stat[j,i] <- res1$ATE_est
      
    }
  
  return(rowMeans(stat^2))
}

emp1 <- c()
for(s in c(3,4,6))
  for(k in c(0.9,0.6,0.3)){
    
    res0 <- mse(spill, design=1, side=s, rho=k)
    res1 <- mse(spill, design=2, side=s, rho=k)
    res2 <- mse(spill, design=3, side=s, rho=k)
    emp1 <- rbind(emp1, res1/res0)
    emp1 <- rbind(emp1, res2/res0)
    
  }
emp1 <- round(emp1,3)

setwd('~/spatio blessing/res/multi_np')
mse <- function(design, side, rho){
  
  filename <- paste0('side',side,'_rho',rho,'_thre0.075.Rdata')
  
  load(filename)
  
  stat <- matrix(0, N_struct, Nsim)
  
  for(j in 1:N_struct)
    
    for(i in 1:Nsim){
      
      res1 <- res_c[[i]][[j]][[design]][[1]]
      stat[j,i] <- res1$ATE_est
      if(is.na(stat[j,i])){stat[j,i] <- 0}
    }
  
  return(rowMeans(stat^2))
  # return(rowMeans(stat))
}

emp2 <- c()
for(s in c(3,4,6))
  for(k in c(0.9, 0.6,0.3)){
    
    res0 <- mse(design=1, side=s, rho=k)
    res1 <- mse(design=2, side=s, rho=k)
    res2 <- mse(design=3, side=s, rho=k)
    emp2 <- rbind(emp1, res1/res0)
    emp2 <- rbind(emp1, res2/res0)
  }
emp2 <- round(emp2,3)
 

cbind(emp1, emp2)


############ figure 6 in the main

filename <- 'D:/ImportantFiles/projects/3. Finished/SpatioBless/figs/'
filename <- paste0(filename,'multi_spill1_thre0.075.pdf')

pdf(filename, height = 7 , width= 12)

par(mai=c(0.6,0.6,0.3,0.3),omi=c(0.1,0.1,0,0),mfcol=c(3,6))

setwd('D:/ImportantFiles/projects/3. Finished/SpatioBless/code/res/multi_para')
strengths <- (seq(0,0.01,0.001)*100)[-2]
N_struct <- 3
teststat <- function(stren, design,side, rho){
  
  filename <- paste0('side',side,'_rho',rho,'_eX0.3.Rdata')
  
  load(filename)
  
  stat <- matrix(0, N_struct, Nsim)
  
  for(j in 1:N_struct)
    
    for(i in 1:Nsim){
      
      ate <- res_c[[i]][[stren]][[design]][[j]][[2]]$ATE_est
      v <- res_c[[i]][[1]][[design]][[j]][[2]]$var_est
      tt <- ate / sqrt(v)
      if(!is.na(tt)){stat[j,i] <- tt}
    }
  
  return(stat)
}
for(k in c(0.9,0.6,0.3))
  for(s in c(6,4,3))
  {
    a <- c(); b <- c(); c <- c()
    for(i in (1:(length(strengths)+1))[-2]){
      res <- teststat(stren=i, design=1, side=s, rho=k)
      a <- rbind(a, rowMeans(res > qnorm(0.95)))
      res <- teststat( stren=i, design=2, side=s, rho=k)
      b <- rbind(b, rowMeans(res > qnorm(0.95)))
      res <- teststat( stren=i, design=3, side=s, rho=k)
      c <- rbind(c, rowMeans(res > qnorm(0.95)))
    }
    if(k==0.9&s==3)
    {plot(strengths, a[,3],type='l',lty=1,ylim = c(0,1),
          ylab = 'P_reject',xlab = 'effect(%)',cex.axis=1.3,cex.lab=1.3)}
    if(k==0.9&s>3)
    {plot(strengths, a[,3],type='l',lty=1,ylim = c(0,1),
          ylab = 'P_reject',xlab = '',xaxt="n",cex.axis=1.3,cex.lab=1.3)}
    if(k<0.9&s==3)
    {plot(strengths, a[,3],type='l',lty=1,ylim = c(0,1),
          ylab = '',xlab = 'effect(%)',yaxt="n",cex.axis=1.3,cex.lab=1.3)}
    if(k<0.9&s>3)
    {plot(strengths, a[,3],type='l',lty=1,ylim = c(0,1),
          ylab = '',xlab = '',yaxt="n",xaxt="n",cex.axis=1.3,cex.lab=1.3)}
    lines(strengths,b[,3],lty=1,col='red')
    lines(strengths,c[,3],lty=1,col='blue')
    lines(strengths,a[,2], lty=2)
    lines(strengths,b[,2],lty=2, col='red')
    lines(strengths,c[,2],lty=2,col='blue')
    lines(strengths,a[,1], lty=3)
    lines(strengths,b[,1],lty=3, col='red')
    lines(strengths,c[,1],lty=3,col='blue')
  }

# strengths <- c(0, 0.005, 0.01, 0.015, 0.02)*100
strengths <- seq(0,0.02,0.0025)*100
setwd('D:/ImportantFiles/projects/3. Finished/SpatioBless/code/res/multi_np')
N_struct <- 3 #3
teststat <- function( stren, design,side, rho){
  
  filename <- paste0('side',side,'_rho',rho,'_thre0.075.Rdata')
  
  load(filename)
  
  stat <- matrix(0, N_struct, Nsim)
  
  for(j in 1:N_struct)
    
    for(i in 1:Nsim){
      
      ate <- res_c[[i]][[j]][[design]][[stren]]$ATE_est
      v <- res_c[[i]][[j]][[design]][[1]]$var_est
      tt <- ate / sqrt(v)
      if(!is.na(tt)){stat[j,i] <- tt}
    }
  
  return(stat)
}
for(k in c(0.9,0.6,0.3))
  for(s in c(6,4,3))
  {
    a <- c(); b <- c(); c <- c()
    for(i in 1:length(strengths)){
      res <- teststat(stren=i, design=1, side=s, rho=k)
      a <- rbind(a, rowMeans(res > qnorm(0.95)))
      res <- teststat( stren=i, design=2, side=s, rho=k)
      b <- rbind(b, rowMeans(res > qnorm(0.95)))
      res <- teststat( stren=i, design=3, side=s, rho=k)
      c <- rbind(c, rowMeans(res > qnorm(0.95)))
    }
    if(k==0.9&s==3)
    {plot(strengths, a[,3],type='l',lty=1,ylim = c(0,1),
          ylab = 'P_reject',xlab = 'effect(%)',cex.axis=1.3,cex.lab=1.3)}
    if(k==0.9&s>3)
    {plot(strengths, a[,3],type='l',lty=1,ylim = c(0,1),
          ylab = 'P_reject',xlab = '',xaxt="n",cex.axis=1.3,cex.lab=1.3)}
    if(k<0.9&s==3)
    {plot(strengths, a[,3],type='l',lty=1,ylim = c(0,1),
          ylab = '',xlab = 'effect(%)',yaxt="n",cex.axis=1.3,cex.lab=1.3)}
    if(k<0.9&s>3)
    {plot(strengths, a[,3],type='l',lty=1,ylim = c(0,1),
          ylab = '',xlab = '',yaxt="n",xaxt="n",cex.axis=1.3,cex.lab=1.3)}
    lines(strengths,b[,3],lty=1,col='red')
    lines(strengths,c[,3],lty=1,col='blue')
    lines(strengths,a[,2], lty=2)
    lines(strengths,b[,2],lty=2, col='red')
    lines(strengths,c[,2],lty=2,col='blue')
    lines(strengths,a[,1], lty=3)
    lines(strengths,b[,1],lty=3, col='red')
    lines(strengths,c[,1],lty=3,col='blue')
  }
dev.off()

