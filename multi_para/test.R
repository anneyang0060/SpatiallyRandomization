setwd('~/spatio blessing/multi_para')
source('FNS_setting.R')
source('FNS_estimation_ATE.R')
setwd('~/spatio blessing/res/multi_para')

### main process
for(rho in c(0.3, 0.6, 0.9))
  
  for(side in c(3,4,6)){
    
    t0 <- Sys.time()
    
    filename <- paste0('side',side,'_rho',rho,'_eX0.3.Rdata')
    
    cl <- makeCluster(Mcl)
    registerDoParallel(cl)
    res_c <- foreach(sd=sds[1:Nsim], .packages = c('MASS','mvnfast')) %dopar% {
      
      ll <- which(sds==sd) 
      res <- list(); ns <- c()
      
      for(stren in strengths){
        res0  <- list(); res1 <- list(); res2 <- list(); st_ns <- c()
        
        for(nr_struc in 1:nrow(st_structure)){
          
          R <- st_structure[nr_struc,1]; M <- st_structure[nr_struc,2]
          st_ns <- c(st_ns, paste0('R',R,'M',M))
          
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
          
          ### set regression coefficients
          {
            set.seed(123)
            alpha <- 8+2*coeff_fn(3,lat,lon) 
            # alpha <- alpha/sum(alpha)*500; 
            beta <- coeff_fn(3,lat,lon)
            gamma_strength <- sum(alpha + EX * beta) * stren
            gamma <- alpha/sum(alpha)*gamma_strength; theta <- beta/sum(beta)*gamma_strength*0.6
            alpha <- matrix(alpha,R,M); beta <- matrix(beta,R,M)
            gamma <- matrix(gamma,R,M); theta <- matrix(theta,R,M)
            set.seed(321)
            alpha_X <- 4+2*coeff_fn(3,lat,lon) #; alpha_X <- alpha_X/sum(alpha_X)*200; 
            beta_X <- coeff_fn(3,lat,lon); beta_X <- (beta_X-min(beta_X))/2/(max(beta_X)-min(beta_X)) + 0.3
            gamma_X_strength <- sum(alpha_X + EX * beta_X) * stren
            gamma_X <- alpha_X/sum(alpha_X)*gamma_X_strength; theta_X <- beta_X/sum(beta_X)*gamma_X_strength*0.6
            alpha_X <- matrix(alpha_X,R,M); beta_X <- matrix(beta_X,R,M)
            gamma_X <- matrix(gamma_X,R,M); theta_X <- matrix(theta_X,R,M)
          } 
          
          ### generate noises and X
          {
            set.seed(sd)
            e <- noise_gene(sd, R, M, rho, sigma, dis_mat)
            e_X <- noise_gene(sd+1e4, R, M, rho, sigma, dis_mat)*0.3
            X <- array(rnorm(R*M*N,EX,DX),dim=c(R,M,N)); X[,2:M,] <- 0
          }
          
          
          {
            # ## global
            # res_temp <- list(model_aggre(R, M, design = 0,spill = 0), 
            #                  model_aggre(R, M, design = 0,spill = 1))
            # names(res_temp) <- c('nospill','spill')
            # res0<- c(res0, list(res_temp))
            # 
            # ## individual
            # res_temp <- list(model_aggre(R, M, design = 1,spill = 0), 
            #                  model_aggre(R, M, design = 1,spill = 1))
            # names(res_temp) <- c('nospill','spill')
            # res1 <- c(res1, list(res_temp))
            # 
            # ## cluster
            # res_temp <- list(model_aggre(R, M, design = 2,spill = 0), 
            #                  model_aggre(R, M, design = 2,spill = 1))
            # names(res_temp) <- c('nospill','spill')
            # res2 <- c(res2, list(res_temp))
            
            ## global
            res_temp <- list(0, model_aggre(R, M, design = 0,spill = 1))
            names(res_temp) <- c('nospill','spill')
            res0<- c(res0, list(res_temp))
            
            ## individual
            res_temp <- list(0, model_aggre(R, M, design = 1,spill = 1))
            names(res_temp) <- c('nospill','spill')
            res1 <- c(res1, list(res_temp))
            
            ## cluster
            res_temp <- list(0, model_aggre(R, M, design = 2,spill = 1))
            names(res_temp) <- c('nospill','spill')
            res2 <- c(res2, list(res_temp))
            
          }
        }
        
        names(res0) <- st_ns; names(res1) <- st_ns; names(res2) <- st_ns
        
        res_givenstrength <- list(res0, res1, res2)
        names(res_givenstrength) <- c('global', 'individual', 'cluster')
        
        res <- c(res, list(res_givenstrength))
        n_stren <- which(strengths == stren)
        ns <- c(ns, paste0('strength',n_stren-1))
      }
      
      names(res) <- ns
      return(res)
    }  
    stopCluster(cl)
    save(res_c,  file = filename)
    print(paste0('side',side,'_rho',rho))
    
    t1 <- Sys.time()
    runtime <- difftime(t1,t0,units = 'mins')
    print(paste0('runtime = ',round(runtime,2),' mins'))
  }

