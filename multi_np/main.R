rm(list = ls())
library(foreach)
library(doParallel)

Nsim <- 50; Mcl <- 50

EX <- 1; DX <- 1
N <- 30
set.seed(123)
Nsds = Nsim
sds <- unique(round(1e5*runif(2*Nsds)))[1:Nsds]
sigma <- 1
Rs <- c(36,81,144)
M <- 12
v_mis <- 0.5
strengths <- seq(0,0.03,0.0025)
d_w <- 11
# 0, 0.005,0.01,0.02,0.05

N_b <- 100
source('~/multi_np/test.R')

