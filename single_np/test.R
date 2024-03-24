setwd('~/spatio blessing/single_np')
source('FNS_setting.R')
source('FNS_estimation.R')
setwd('~/spatio blessing/res/single_np')

### main process
# mode <- 1; side <- 4

for(mode in c(1,2,3))

  for(rho in c(0.9, 0.6, 0.3))

    for(side in c(3,4,6))
      # mode= 1; rho = 0.8; side = 3
    {
      
      t0 <- Sys.time()
      
      filename <- paste0('mode',mode,'/side',side,'_rho',rho,'_thre0.025_e1.Rdata')
      
      cl <- makeCluster(Mcl)
      registerDoParallel(cl)
      res_c <- foreach(sd=sds[1:Nsim], .packages = c('MASS' ,'mvnfast')) %dopar% {
        
        ll <- which(sds==sd) 
        res <- list(); ns <- c()
        
        for(stren in strengths){
          
          res0 <- list(); res1 <- list(); res2 <- list(); st_ns <- c()
          if(stren == 0){N_b <- N_b0}else{N_b <- 0 }

          for(R in Rs){
            # stren = strengths[4] ; R=50
            st_ns <- c(st_ns, paste0('R',R))
            
            ### generate regions
            {
              mylist <- region_gene(4, R)
              coor <- mylist$coor
              lat <- (coor[,1]-min(coor[,1]))/(max(coor[,1])-min(coor[,1]))
              lon <- (coor[,2]-min(coor[,2]))/(max(coor[,2])-min(coor[,2]))
              dis_mat <- c()
              for(r in 1:R){
                dis <- sapply(1:R, function(i){0.5*sqrt((lat[i]-lat[r])^2+(lon[i]-lon[r])^2)})
                dis_mat <- rbind(dis_mat,dis)
              }
              mylist <- region_gene(side,R)
              coor <- mylist$coor
              lat <- (coor[,1]-min(coor[,1]))/(max(coor[,1])-min(coor[,1]))
              lon <- (coor[,2]-min(coor[,2]))/(max(coor[,2])-min(coor[,2]))
              nb <- mylist$nb # neighbors of each region
              cluster <- mylist$cluster
              Nc_group <- mylist$Nc_group
              N_cluster <- R/9
              rm(mylist)
            }
            
            ### generate noises and X
            {
              set.seed(sd)
              e <- noise_gene(mode, sd, R, M, rho, sigma, dis_mat)
              X <- rnorm(5*R*N,EX,DX)
              X <- X[which(X>3&X<5)]
              X <- matrix(X[1:(R*N)],R,N)
            }
            
            {
              
              ## global
              res_temp <- list(model_aggre(X, R, design = 0,spill = 0, sd, nb),
                               model_aggre(X, R, design = 0,spill = 1, sd, nb))
              names(res_temp) <- c('nospill', 'spill')
              res0 <- c(res0, list(res_temp))

              ## individual randomization
              res_temp <- list(model_aggre(X, R, design = 1,spill = 0, sd, nb),
                               model_aggre(X, R, design = 1,spill = 1, sd, nb))
              names(res_temp) <- c('nospill','spill')
              res1 <- c(res1, list(res_temp))
              
              ## cluster
              res_temp <- list(model_aggre(X, R, design = 2,spill = 0, sd, nb),
                               model_aggre(X, R, design = 2,spill = 1, sd, nb))
              names(res_temp) <- c('nospill','spill')
              res2 <- c(res2, list(res_temp))

            }
            
          }
          
          names(res0) <- st_ns; names(res1) <- st_ns; names(res2) <- st_ns
          
          res_givenstrength <- list(res0, res1, res2)
          names(res_givenstrength) <- c('global', 'individual','cluster')
          
          res <- c(res, list(res_givenstrength))
          n_stren <- which(strengths == stren)
          ns <- c(ns, paste0('strength',n_stren-1))
        }
        
        names(res) <- ns
        
        return(res)
      }  
      stopCluster(cl)
      save(res_c,  file = filename)
      print(paste0('mode',mode,'_side',side,'_rho',rho))
      
      t1 <- Sys.time()
      runtime <- difftime(t1,t0,units = 'mins')
      print(paste0('runtime = ',round(runtime,2),' mins'))
    }
  
