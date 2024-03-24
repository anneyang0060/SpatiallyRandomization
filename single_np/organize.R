# setwd('D:/ImportantFiles/projects/3. Finished/SpatioBless/code/res/single_np')
setwd('~/spatio blessing/res/single_np')

############ parameters and functions
Nsim <- 500

Rs <- c(36, 81, 144)

N_struct <- 3

mse <- function(spill,  design, mode, side, rho){
  
  filename <- paste0('mode',mode,'/side',side,'_rho',rho,'_thre0.025_e1.Rdata')
  
  load(filename)
  
  stat <- matrix(0, N_struct, Nsim)
  mse1 <- rep(0, N_struct)
  
  for(j in 1:N_struct){
    
    for(i in 1:Nsim){
      
      res1 <- res_c[[i]][[1]][[design]][[j]][[spill]]
      stat[j,i] <- res1$ATE_est
      
    }
    
    mse1[j] <- mean(stat[j, ]^2)
  }
  return(mse1)
}


############# table 3 and 4
for(spill in c(2,1)){
  emp1 <- c()
  
  for(j in 1:3)
    for(s in c(3,4,6))
      for(k in c(0.9,0.6,0.3)){
        
        res0 <- mse(spill, design=1, mode=j, side=s, rho=k)
        res1 <- mse(spill, design=2, mode=j, side=s, rho=k)
        res2 <- mse(spill, design=3, mode=j, side=s, rho=k)
        emp1 <- rbind(emp1, res1/res0)
        emp1 <- rbind(emp1, res2/res0)
        
      }
  emp1 <- round(emp1,3)
  emp <- matrix(0,18,9)
  emp[,1:3] <- emp1[1:18,]
  emp[,4:6] <- emp1[19:36,]
  emp[,7:9] <- emp1[37:54,]
  print(emp)
}



########## figure 5 in the main; Fig 3 and 4 in the supplement  
strengths <- 100*c(0, 0.0025,0.005, 0.0075, 0.01, 0.0125,0.015, 0.02, 0.025, 0.03)[1:7]

teststat <- function(spill, stren, design, mode, side, rho){
  
  filename <- paste0('mode',mode,'/side',side,'_rho',rho,'_thre0.025_e1.Rdata')
  
  load(filename)
  
  stat <- matrix(0, N_struct, Nsim)
  
  for(j in 1:N_struct)
    
    for(i in 1:Nsim){
      
      res0 <- res_c[[i]][[1]][[design]][[j]][[spill]]
      res1 <- res_c[[i]][[stren]][[design]][[j]][[spill]]
      stat[j,i] <- res1$ATE_est / sqrt(res0$var_est)
      
    }
  
  return(stat)
}

for(j in 1:3){
filename <- 'D:/ImportantFiles/projects/3. Finished/SpatioBless/figs/'
filename <- paste0(filename,'single_np_mode',j,'_thre0.025.pdf')
pdf(filename, height = 7 , width= 12)
par(mai=c(0.6,0.6,0.3,0.3),omi=c(0.1,0.1,0,0),mfcol=c(3,6))
for(k in c(0.9,0.6,0.3))
  for(s in c(6,4,3))
  {
    a <- c(); b <- c(); c<- c()
    for(i in 1:7){
      res <- teststat(spill=2, stren=i, design=1, mode=j, side=s, rho=k)
      a <- rbind(a, rowMeans(res > qnorm(0.95)))
      res <- teststat(spill=2, stren=i, design=2, mode=j, side=s, rho=k)
      b <- rbind(b, rowMeans(res > qnorm(0.95)))
      res <- teststat(spill=2, stren=i, design=3, mode=j, side=s, rho=k)
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
    lines(strengths,c[,1],lty=3, col='blue')
  }
for(k in c(0.9,0.6,0.3))
  for(s in c(6,4,3))
  {
    a <- c(); b <- c(); c <- c()
    for(i in 1:7){
      res <- teststat(spill=1, stren=i, design=1, mode=j, side=s, rho=k)
      a <- rbind(a, rowMeans(res > qnorm(0.95)))
      res <- teststat(spill=1, stren=i, design=2, mode=j, side=s, rho=k)
      b <- rbind(b, rowMeans(res > qnorm(0.95)))
      res <- teststat(spill=1, stren=i, design=3, mode=j, side=s, rho=k)
      c <- rbind(c, rowMeans(res > qnorm(0.95)))
    }
    if(k==0.9&s==3)
    {plot(strengths, a[,3],type='l',lty=1,ylim = c(0,1),#col='blue',
          ylab = 'P_reject',xlab = 'effect(%)',cex.axis=1.3,cex.lab=1.3)}
    if(k==0.9&s>3)
    {plot(strengths, a[,3],type='l',lty=1,ylim = c(0,1),#col='blue',
          ylab = 'P_reject',xlab = '',xaxt="n",cex.axis=1.3,cex.lab=1.3)}
    if(k<0.9&s==3)
    {plot(strengths, a[,3],type='l',lty=1,ylim = c(0,1),#col='blue',
          ylab = '',xlab = 'effect(%)',yaxt="n",cex.axis=1.3,cex.lab=1.3)}
    if(k<0.9&s>3)
    {plot(strengths, a[,3],type='l',lty=1,ylim = c(0,1),#col='blue',
          ylab = '',xlab = '',yaxt="n",xaxt="n",cex.axis=1.3,cex.lab=1.3)}
    lines(strengths,b[,3],lty=1,col='red')
    lines(strengths,c[,3],lty=1,col='blue')
    lines(strengths,a[,2], lty=2)
    lines(strengths,b[,2],lty=2, col='red')
    lines(strengths,c[,2],lty=2,col='blue')
    lines(strengths,a[,1], lty=3)
    lines(strengths,b[,1],lty=3, col='red')
    lines(strengths,c[,1],lty=3, col='blue')
  }
dev.off()
}

