setwd('~/spatio blessing/multi_np')
source('FNS_setting.R')
source('FNS_estimation.R')
setwd('~/spatio blessing/res/multi_np')

### main process
# Rs <- c(36)
# M <- 12
# v_mis <- 0.5
# strengths <- c(0, 0.01, 0.02, 0.03)
# mode <- 2
# rho = 0.9
# Mcl = 50; Nsim = 50
# for(side in c(3,4,6))
mode <- 2
side <- 6
for(rho in c(0.3,0.6,0.9))
  
  # for(side in c(3,4,6))
    # R <- 36; stren <- 0; sd=sds[1]
    # rho <- 0.9; side <- 3; Mcl <-50; Nsim <-50
  {
    t0 <- Sys.time()
    
    filename <- paste0('side',side,'_rho',rho,'_thre0.075.Rdata')
    
    cl <- makeCluster(Mcl)
    registerDoParallel(cl)
    res_c <- foreach(sd=sds[1:Nsim], .packages = c('MASS' ,'mvnfast','Rcpp')) %dopar% {
      
      sourceCpp('LC.cpp')
      sourceCpp('getWeight1.cpp')
      ll <- which(sds==sd) 
      res <- list(); st_ns <- c()
      mu0_0 <- 0; mu1_0 <- 0; mu0_1 <- 0; mu1_1 <- 0
      
      for(R in Rs){
        
        res0 <- list(); res1 <- list(); res2 <- list(); ns <- c()
        st_ns <- c(st_ns, paste0('R',R,'M',M))
        
        for(stren in strengths){
          
          # t0 <- Sys.time()
          
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
          
          ### generate noises
          {
            set.seed(sd)
            e <- 0.5 * noise_gene(sd, R, M, rho, sigma, dis_mat)
            E <- 0.1 * noise_gene(sd, R, M, rho, sigma, dis_mat)
          }
          
          ### estimate
          {
            ## global
            res_g <- model_aggre(R, design = 0,spill = 1, sd, nb, mu1_0, mu0_0)
            res0 <- c(res0, list(res_g[1:2]))
            
            ## individual
            res_i <- model_aggre(R, design = 1,spill = 1, sd, nb, mu1_1, mu0_1)
            res1 <- c(res1, list(res_i[1:2]))
            
            ## cluster
            res_c <- model_aggre(R, design = 2,spill = 1, sd, nb, mu1_1, mu0_1)
            res2 <- c(res2, list(res_c[1:2]))
            
            # if(stren==0){
            #   mu1_0 <- res_g$mu1; mu0_0 <- res_g$mu0
            #   mu1_1 <- res_i$mu1; mu0_1 <- res_i$mu0
            #   mu1_2 <- res_c$mu1; mu0_2 <- res_c$mu0
            #   rm(res_g, res_i, res_c)
            # }
          }
          
          # t1 <- Sys.time()
          # runtime <- difftime(t1,t0,units = 'mins')
          # print(paste0('runtime = ',round(runtime,2),' mins'))
          
          n_stren <- which(strengths == stren)
          ns <- c(ns, paste0('strength',n_stren-1))
        }
        
        names(res0) <- ns; names(res1) <- ns; names(res2) <- ns
        
        res_givenR <- list(res0, res1, res2)
        names(res_givenR) <- c('global', 'individual', 'cluster')
        
        res <- c(res, list(res_givenR))
      }
      
      names(res) <- st_ns
      
      return(res)
    }  
    stopCluster(cl)
    
    save(res_c,  file = filename)
    print(paste0('mode',mode,'_side',side,'_rho',rho))
    
    t1 <- Sys.time()
    runtime <- difftime(t1,t0,units = 'mins')
    print(paste0('runtime = ',round(runtime,2),' mins'))
  }

# est0_st<-c();for(i in 1:Nsim){est0_st<-c(est0_st,res_c[[i]][[1]][[2]][[1]][[1]])}
# est0_city<-c();for(i in 1:Nsim){est0_city<-c(est0_city,res_c[[i]][[1]][[1]][[1]][[1]])}
# est1_st<-c();for(i in 1:Nsim){est1_st<-c(est1_st,res_c[[i]][[1]][[2]][[2]][[1]])}
# est1_city<-c();for(i in 1:Nsim){est1_city<-c(est1_city,res_c[[i]][[1]][[1]][[2]][[1]])}
# est2_st<-c();for(i in 1:Nsim){est2_st<-c(est2_st,res_c[[i]][[1]][[2]][[3]][[1]])}
# est2_city<-c();for(i in 1:Nsim){est2_city<-c(est2_city,res_c[[i]][[1]][[1]][[3]][[1]])}
# v_city<-c();for(i in 1:Nsim){v_city<-c(v_city,res_c[[i]][[1]][[1]][[1]][[2]])}
# v_st<-c();for(i in 1:Nsim){v_st<-c(v_st,res_c[[i]][[1]][[2]][[1]][[2]])}
# mean(est0_city/sqrt(v_city)>qnorm(0.95))
# mean(est1_city/sqrt(v_city)>qnorm(0.95))
# mean(est2_city/sqrt(v_city)>qnorm(0.95))
# mean(est0_st/sqrt(v_st)>qnorm(0.95))
# mean(est1_st/sqrt(v_st)>qnorm(0.95))
# mean(est2_st/sqrt(v_st)>qnorm(0.95))
