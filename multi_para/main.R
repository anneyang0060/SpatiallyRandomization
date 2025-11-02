rm(list = ls())
library(foreach)
library(doParallel)

Nsim <- 500; Mcl <- 100

N <- 30
# noise_scale <- 0.5
set.seed(123)
Nsds <- 500
sds <- unique(round(1e5*runif(1000)))[1:Nsds]
sigma <- 1
Rs <- c(36, 81, 144)
Ms <- c(12)
st_structure <- cbind(rep(Rs,each=length(Ms)), rep(Ms,length(Rs)))
EX <- 4; DX <- 1
strengths <- seq(0,0.02,0.001)
mis_specify = 0
Nb <- 100
source('~/multi_para/test.R')





