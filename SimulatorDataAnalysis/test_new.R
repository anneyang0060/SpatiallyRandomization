library(MASS)
library(foreach)
library(doParallel)

setwd('~/SimulatorDataAnalysis')
source('FNSnew.R')


N <- 20
R <- 64
Nb <- 100
stat_date <- 20161109



#################### parametric
Nsim <- 100
Mcl <- 50
Rs <- (1:R)
M <- 24
M0 <- M
### compute lat, lon and nb
r <- sqrt(R)
lat <- rep((2*(1:r)-1) * (1/(2*r)), r)
lon <- rep((2*(1:r)-1) * (1/(2*r)), each = r)
nb <- list() #neighbors of each region
dis_mat <- c()
for(r in 1:R){
  dis <- sapply(1:R, function(i){sqrt((lat[i] - lat[r])^2 + (lon[i] - lon[r])^2)})
  dis_mat <- rbind(dis_mat,dis)
  nb1 <- which(dis < 1.01*(1/sqrt(R)))
  nb1 <- nb1[nb1!=r]
  nb[[r]] <- nb1
}

######## power
for(DD in c(0.7,0.71,0.72,0.75))

  for(design in c(0,1))

    for(M in c(24))

      for(spill in c(0,1)){

        TI <- 6

        cl <- makeCluster(Mcl)
        registerDoParallel(cl)
        res_c <- foreach(sd=0:(Nsim-1), .packages = c('MASS')) %dopar% {

          ### read in data
          if(DD>0.7){
            load(paste0('./data/design',design,'_R',R,'_M',M,'_DD',DD,'_spill',spill,'/sd',sd,'.Rdata'))
          }else{
            load(paste0('./data/design',0,'_R',R,'_M',M,'_DD',DD,'_spill',spill,'/sd',sd,'.Rdata'))
          }

          ### estimate
          ATE_hat <- OLS_ATE_est(y,design,spill,A,X,nb)

          ### bootstrap
          ATE_b <- c()
          DE_b <- c();IE_b <- c()
          for(i_b in 1:Nb){
            set.seed(i_b)
            idx_b <- sample(1:N, N, replace = TRUE)
            X_b <- X[,,idx_b] ; y_b <- y[,,idx_b]; A_b <- A[,,idx_b]
            # ATE_b <- c(ATE_b, OLS_ATE_est(y_b,design,spill,A_b,X_b,nb))
            res_est <- OLS_ATE_est(y_b,design,spill,A_b,X_b,nb)
            ATE_b <- c(ATE_b,res_est[1]); DE_b <- c(DE_b,res_est[2]); IE_b <- c(IE_b,res_est[3])
          }
          # var_hat <- var(ATE_b)
          var_hat <- c(var(ATE_b),var(DE_b),var(IE_b))
          res <- c(ATE_hat,var_hat)

          return(res)

        }
        stopImplicitCluster()
        stopCluster(cl)

        save(res_c, file = paste0('res/design',design,'_R',R,'_M',M,'_DD',DD,'_spill',spill,'_5.Rdata'))

      }


for(DD in c(0.7, 0.71, 0.72, 0.75))
  
  for(design in c(0,1))
    
    for(M in c(24))
      
      for(spill in c(0,1)){
        
        filename <- paste0('res/design',design,'_R',R,'_M',M,'_DD',DD,'_spill',spill,'_5.Rdata')
        
        load(filename)
        
        # test_stat <- sapply(1:100, function(i){res_c[[i]]})
        # print(mean(test_stat>0.95))
        
        test_stat1 <- sapply(1:100, function(i){res_c[[i]][1]/sqrt(res_c[[i]][4])})
        test_stat2 <- sapply(1:100, function(i){res_c[[i]][2]/sqrt(res_c[[i]][5])})
        test_stat3 <- sapply(1:100, function(i){res_c[[i]][3]/sqrt(res_c[[i]][6])})

        print(mean(test_stat1>qnorm(0.95)))
        # print(mean(test_stat1>qnorm(0.95)))
        # print(c(mean(test_stat1>qnorm(0.95)),mean(test_stat2>qnorm(0.95)),mean(test_stat3>qnorm(0.95))))
      }


####### mse
load('data/res_mc.Rdata')
ate_mc <- colMeans(effect_mc)
M <- 24; gmv <- mean(gmv)

DD <- 0.71
idx <- which(c(0.71,0.72,0.75)==DD)
mse <- c()

for(design in c(0,1))
  
  for(spill in c(0,1)){
    
    filename <- paste0('res/design',design,'_R',R,'_M',M,'_DD',DD,'_spill',spill,'_5.Rdata')
    
    load(filename)
    
    ate_hat <- sapply(1:100, function(i){res_c[[i]][1]})
    var_hat <- sapply(1:100, function(i){res_c[[i]][4]})
    mse <- c(mse, mean(var_hat + (ate_hat-ate_mc[idx])^2))
    # mse <- c(mse, var(ate_hat)+(mean(ate_hat)-ate_mc[idx])^2)
    
  }

mse[3:4]/mse[1:2]

DD <- 0.7
mse <- c()

for(design in c(0,1))
  
  for(spill in c(0,1)){
    
    filename <- paste0('res/design',design,'_R',R,'_M',M,'_DD',DD,'_spill',spill,'_5.Rdata')
    
    load(filename)
    
    ate_hat <- sapply(1:100, function(i){res_c[[i]][1]})
    var_hat <- sapply(1:100, function(i){res_c[[i]][4]})
    mse <- c(mse, mean(var_hat + (ate_hat)^2))
    # mse <- c(mse, var(ate_hat)+(mean(ate_hat))^2)

  }

mse[3:4]/mse[1:2]

