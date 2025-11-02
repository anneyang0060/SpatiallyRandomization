library(foreach)
library(doParallel)

Nsim <- 500; Mcl <- 100

EX <- 4; DX <- 1
N <- 30
# noise_scale <- 0.5
set.seed(123)
Nsds = Nsim
sds <- unique(round(1e5*runif(2*Nsds)))[1:Nsds]
sigma <- 1
Rs <- c(36, 81, 144)
M <- 1
v_mis <- 0.5
strengths <- seq(0, 0.02, 0.0025)

source('~/single_para/test.R')



