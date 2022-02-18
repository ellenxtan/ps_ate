## Method for estimating a principal stratum causal effect
## Main function wrapping up the helper functions

Main <- function(beta_true) {
  #> Parameter Settings ---------- 
  pars_init = list()
  pars_init$beta_true = beta_true
  pars_init$r = 1000  # number of dataset created
  pars_init$b = 500  # number of bootstrapping samples
  pars_init$n = 4000 # number of simulated subjects
  pars_init$m = 4  # level of Xi, choose 3 or 4 for now
  pars_init$beta_init = pars_init$beta_true
  pars_init$optimConst = 1#pars_init$n
  pars_init$optimMethod = "BFGS"
  pars_init$lambda = 0
  
  cat("n=", pars_init$n, "b=", pars_init$b, "r=", pars_init$r, 
      "m=", pars_init$m, "beta_true=", pars_init$beta_true, "\n")
  
  #> Data Settings ----------
  pars_init$pz = c(.5, .5) # true:(.5,.5)
  stopifnot(sum(pars_init$pz)==1)
  pars_init$px = rep(0.25, pars_init$m) #c(.25, .2, .17, .38) # true:(.25,.20,.17,.38)
  
  pars_init$pS01_x = c(0.3, 0.25, 0.25, 0.2) # true:(.28,.23,.16,.22)
  pars_init$pS00_x = 1 - pars_init$pS01_x
  stopifnot(sum(pars_init$pS01_x)==1)
  
  betaFn = function(x, beta, y){
    return(exp(sum(beta*c(1,y,x)))/(1+exp(sum(beta*c(1,y,x)))))
  }
  X = seq(pars_init$m)-1
  pars_init$pS11_S00Y00x = sapply(X, betaFn, beta=pars_init$beta_true, y=0)  # true 5yr:(0.00619 0.16281 0.32567 0.08709)
  pars_init$pS11_S00Y01x = sapply(X, betaFn, beta=pars_init$beta_true, y=1)  # true 5yr:(0.00031 0.00974 0.02384 0.0048)
  
  pars_init$pY11_S00S10Y00x = 0.50
  pars_init$pY11_S00S10Y01x = 0.60
  pars_init$pY11_S00S11Y00x = 0.85  #rep(0.85, m)
  pars_init$pY11_S00S11Y01x = 0.90  #rep(0.90, m)
  pars_init$pY11_S01S11Y00x = 0.85  #rep(0.85, m)
  pars_init$pY11_S01S11Y01x = 0.90  #rep(0.90, m)
  
  pars_init$pY01_S00x = c(0.7, 0.65, 0.60, 0.55)  # true:(.77,.78,.66,.6)
  pars_init$pY01_S01x = pars_init$pY01_S00x * 1.2#c(.93,.81,.94,.82)# true:(.93,.81,.94,.82)
  # increase = 1.2
  # pars_init$pY01_S01x = increase * pars_init$pY01_S00x
  pars_init$pY00_S01x = 1 - pars_init$pY01_S01x
  pars_init$pY00_S00x = 1 - pars_init$pY01_S00x
  
  #> ThetaTrue ----------
  GetThetaTrue(pars_init)
  
  #> Parallel Computing ----------
  library(foreach)
  library(doParallel)
  # choose cluster numbers
  cl <- parallel::makeCluster(detectCores(logical = TRUE))
  print(cl)
  registerDoParallel(cl)
  
  resBoot = Estimation(pars_init, GetThetaTrue, GetSimulateData, EstimateBeta, EstimateEffect, GetOptimBeta_sum)
  bootTheta = resBoot$bootTheta
  thetaTrue = resBoot$thetaTrue
  dataTheta = resBoot$dataTheta
  resAnalysis = Analysis(dataTheta, bootTheta, thetaTrue, pars_init$n, pars_init$b, pars_init$r)
  
  filename = paste0("sim_", pars_init$beta_true[1], "_", 
                    pars_init$beta_true[2], "_",
                    pars_init$beta_true[3], "_",
                    "n",pars_init$n,"b",pars_init$b,"r",pars_init$r, "_",
                    Sys.Date(), ".RData")
  print(filename)
  save.image(filename)
  
  cat("Finished!\n")
}


# end parallel
#stopCluster(cl)

