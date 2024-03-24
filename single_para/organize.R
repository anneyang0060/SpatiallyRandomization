# setwd('D:/ImportantFiles/projects/3. Finished/SpatioBless/code/res/single_para')
setwd('~/spatio blessing/res/single_para')

############ parameters and functions

Nsim <- 500

Rs <- c(36, 81, 144)

N_struct <- 3

mse <- function(spill, design, mode, side, rho){
  
  filename <- paste0('mode',mode,'/side',side,'_rho',rho,'_.Rdata')
  
  load(filename)
  
  stat <- matrix(0, N_struct, Nsim)
  mse1 <- rep(0, N_struct)
  
  for(j in 1:N_struct){
    
    for(i in 1:Nsim){
      
      res1 <- res_c[[i]][[1]][[design]][[j]][[spill]]
      stat[j,i] <- res1$ATE_est

    }
    
  }

  return(rowMeans(stat^2))
}

############# table 1
for(j in 1:3)
  for(s in c(3,4,6))
    for(k in c(0.9,0.6,0.3)){
      
      res0 <- mse(spill=1, design=1, mode=j, side=s, rho=k)
      res1 <- mse(spill=1, design=2, mode=j, side=s, rho=k)
      res2 <- mse(spill=1, design=3, mode=j, side=s, rho=k)
      emp2 <- rbind(emp2, res1/res0)
      emp2 <- rbind(emp2, res2/res0)
      
    }
emp2 <- round(emp2,3)

############# table 2
emp1 <- c(); emp2 <- c()

for(j in 1:3)
  for(s in c(3,4,6))
    for(k in c(0.9,0.6,0.3)){
      
      res0 <- mse(spill=2, design=1, mode=j, side=s, rho=k)
      res1 <- mse(spill=2, design=2, mode=j, side=s, rho=k)
      res2 <- mse(spill=2, design=3, mode=j, side=s, rho=k)
      emp1 <- rbind(emp1, res1/res0)
      emp1 <- rbind(emp1, res2/res0)
    }
emp1 <- round(emp1,3)

# emp <- matrix(0, 36, 9)
# emp[1:18,1:3] <- emp1[1:18, ]
# emp[1:18,4:6] <- emp1[19:36, ]
# emp[1:18,7:9] <- emp1[37:54, ]
# emp[19:36,1:3] <- emp2[1:18, ]
# emp[19:36,4:6] <- emp2[19:36, ]
# emp[19:36,7:9] <- emp2[37:54, ]
# emp

########## figure 4 in the main; figure 1 and 2 in the supplement
teststat <- function(spill, stren, design, mode, side, rho){
  
  filename <- paste0('mode',mode,'/side',side,'_rho',rho,'_.Rdata')
  
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

strengths <- seq(0, 0.02, 0.0025)
strengths <- strengths*100

for(j in 1:3){
filename <- 'D:/ImportantFiles/projects/3. Finished/SpatioBless/figs/'
# filename <- '~/spatio blessing/figs/'
filename <- paste0(filename,'single_para_mode',j,'.pdf')
pdf(filename, height = 7 , width= 12)
par(mai=c(0.6,0.6,0.3,0.3),omi=c(0.1,0.1,0,0),mfcol=c(3,6))
for(k in c(0.9,0.6,0.3))
  for(s in c(6,4,3))
  {
    a <- c(); b <- c(); c <- c()
    for(i in 1:9){
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
    lines(strengths,c[,1],lty=3,col='blue')
  }
for(k in c(0.9,0.6,0.3))
  for(s in c(6,4,3))
  {
    a <- c(); b <- c(); c <- c()
    for(i in 1:9){
      res <- teststat(spill=1, stren=i, design=1, mode=j, side=s, rho=k)
      a <- rbind(a, rowMeans(res > qnorm(0.95)))
      res <- teststat(spill=1, stren=i, design=2, mode=j, side=s, rho=k)
      b <- rbind(b, rowMeans(res > qnorm(0.95)))
      res <- teststat(spill=1, stren=i, design=3, mode=j, side=s, rho=k)
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
}



######### generate fig 2
delta_i <- 0.2
delta_c <- 0.1
alpha_i <- 0.9
alpha_c <- 0.6

ratio_i <- function(R,rho,r,mode){
  if(mode==1){return(r/R/rho+(r+1)^2/R)}
  if(mode==2){return((1+r)^2/R/rho^delta_i)}
  if(mode==3){return((1+r)^2/R/rho*alpha_i)}
}
ratio_c <- function(R,rho,r,c,mode){
  if(mode==1){return(r/sqrt(c)/R/rho+(sqrt(c)+r)^2/R)}
  if(mode==2){return((sqrt(c)+r)^2/R/rho^delta_c)}
  if(mode==3){return((sqrt(c)+r)^2/R/rho*alpha_c)}
}

c <- 9
r <- 4
Rs <- 50:300
filename <- 'D:/ImportantFiles/projects/3. Finished/SpatioBless/figs/illustration.pdf'
# filename <- '~/spatio blessing/figs/illustration.pdf'
pdf(filename, height = 3 , width= 9)
par(mai=c(0.8,0.8,0.2,0.2),omi=c(0.1,0.1,0,0),mfrow=c(1,3))
for(mode in c(1,2,3)){
  res_i <- c(); res_c <- c()
  for(rho in c(0.3,0.6,0.9)){
    res_i <- rbind(res_i, sapply(Rs, function(R){ratio_i(R,rho,r,mode)}))
    res_c <- rbind(res_c, sapply(Rs, function(R){ratio_c(R,rho,r,c,mode)}))
  }
  plot(Rs, res_i[1,],type = 'l',ylim = range(rbind(res_i,res_c)),
       ylab = 'mse ratios',xlab = 'R')
  lines(Rs, res_c[1,],lty = 2)
  lines(Rs, res_i[2,],col = 'blue')
  lines(Rs, res_c[2,],col = 'blue', lty = 2)
  lines(Rs, res_i[3,],col = 'red')
  lines(Rs, res_c[3,],col = 'red', lty = 2)
}
dev.off()

# ######### generate fig
# tau <- seq(0,10,0.01)
# d1 <- 1
# d2 <- 3*d1
# p1 <- 1-pnorm(qnorm(0.95)-tau/d1)
# p2 <- 1-pnorm(qnorm(0.95)-tau/d2)
# 
# 
# # filename <- 'D:/ImportantFiles/projects/3. Finished/SpatioBless/figs/illustration_power_compare.pdf'
# filename <- '~/spatio blessing/figs/illustration_power_compare.pdf'
# pdf(filename, height = 4 , width= 6)
# par(mai=c(0.8,0.8,0.2,0.2),omi=c(0.1,0.1,0,0),mfrow=c(1,1))
# plot(tau,p1,type='l',xlab = 'signal strength',ylab = 'power',
#      col='red')
# lines(tau,p2,col='blue')
# dev.off()

