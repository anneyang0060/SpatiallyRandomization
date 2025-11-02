library(MASS)

setwd('C:/Users/annie/OneDrive/ImportantFiles/projects/3. Finished/SpatioBless/code and data/RealDataAnalysis')
source('FNS.R')

data <- read.csv('V1_CityE_expand_AB.csv')
M <- length(unique(data$time))
R <- length(unique(data$region))
N <- length(unique(data$date))
adj_mat <- matrix(c(0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
          1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0,
          1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
          0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0,
          0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0,
          0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1,
          0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0),R,R)


A <- reformat(data$is_exp)
X <- reformat(data$cnt_call)
y <- reformat(data$gmv)

nb <- list() 
for(r in 1:R){
  nb1 <- which(adj_mat[r,]!=0)
  nb[[r]] <- nb1
}

eps_hat <- fitted_resid(y, design=1, spill=1, A, X, nb)
eps_hat_city <- apply(eps_hat, MARGIN=c(2,3), sum)
y_hat <- y[,2:M,]-eps_hat
y_hat_city <- apply(y_hat, MARGIN=c(2,3), sum)
y_city <- apply(y[,2:M,], MARGIN=c(2,3), sum)

color <- c('lightblue','lightgray','lightgreen','lightpink','lightsalmon','thistle','lightsteelblue')

pdf('realdata.pdf', height = 2.2 , width= 10)
par(mai=c(0.8,0.8,0.2,0.2),omi=c(0.1,0.1,0,0),mfrow=c(1,4))
for(r in c(1,11)){
  
  y1 <- as.vector(y[r,2:M,])
  y1_hat <- as.vector(y_hat[r,,])
  std <- sd(y1)
  plot(eps_hat[r,,1]/std,type='l', ylim = range(eps_hat[r,,]/std),
       xlab = 'time', ylab = 'fitted residuals',lwd=1.3,col=color[1])
  for(i in 2:N){lines(eps_hat[r,,i]/std,col=color[i%%7+1],lwd=1.3)}  
  
  plot(y1,y1_hat, pch=16, col='lightsteelblue',
       xlab = 'y', ylab = 'y_hat')
  lines(0:max(y1),0:max(y1),lwd=2)
  
  hist(y1, main=r)
}
dev.off()

Nb <- 200
res <- stat_est(y, design=1, spill=1, A, nb, X)
p_val <- 1-pnorm(res$ATE_est/sqrt(res$var_est))  


