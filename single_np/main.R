library(foreach)
library(doParallel)

Nsim <- 500; Mcl <- 50

EX <- 4; DX <- 1
# N <- 100
noise_scale <- 0.5
set.seed(123)
Nsds = Nsim
sds <- unique(round(1e5*runif(2*Nsds)))[1:Nsds]
sigma <- 1
Rs <- c(36,81,144)
M <- 1
strengths <- c(0, 0.0025,0.005, 0.0075, 0.01, 0.0125,0.015, 0.02, 0.025, 0.03)

N_b0 <- 100

N <- 30
source('~/spatio blessing/single_np/test.R')

