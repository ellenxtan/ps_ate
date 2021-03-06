# setwd("C:/Users/xiaoq/Box/lab/proj_causal_infer/code/simulation_code")
rm(list=ls())
set.seed(20190712)
source("sim_main.R")
source("sim_fn.R")

beta_true = c(-3,-5,0.2)
Main(beta_true)

beta_true = c(-5,-1,-2)
Main(beta_true)

beta_true = c(-7,3,0.2)
Main(beta_true)
