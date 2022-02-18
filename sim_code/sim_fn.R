## Method for estimating a principal stratum causal effect
## Helper functions

#> Calculate ground truth theta -----------------------------------------------
GetThetaTrue = function(pars_init) {  # m: number of Xi category
  pz = pars_init$pz
  px = pars_init$px
  
  pS01_x = pars_init$pS01_x
  pS00_x = pars_init$pS00_x

  pS11_S00Y00x = pars_init$pS11_S00Y00x
  pS11_S00Y01x = pars_init$pS11_S00Y01x
  
  pY11_S00S11Y00x = pars_init$pY11_S00S11Y00x
  pY11_S00S11Y01x = pars_init$pY11_S00S11Y01x
  pY11_S01S11Y00x = pars_init$pY11_S01S11Y00x
  pY11_S01S11Y01x = pars_init$pY11_S01S11Y01x
  
  pY01_S00x = pars_init$pY01_S00x
  pY01_S01x = pars_init$pY01_S01x
  pY00_S01x = pars_init$pY00_S01x
  pY00_S00x = pars_init$pY00_S00x
  
  p1 = sum(px*pS01_x + px*pS00_x*pY00_S00x*pS11_S00Y00x +
             px*pS00_x*pY01_S00x*pS11_S00Y01x)
  p2 = sum( px*pS01_x*pY01_S01x + px*pS00_x*pY01_S00x*pS11_S00Y01x)
  p3 = sum(px*pS01_x*pY00_S01x*pY11_S01S11Y00x + px*pS00_x*pY00_S00x*pS11_S00Y00x*pY11_S00S11Y00x +
             px*pS01_x*pY01_S01x*pY11_S01S11Y01x + px*pS00_x*pY01_S00x*pS11_S00Y01x*pY11_S00S11Y01x)
  thetaTrue = (p3-p2)/p1
  
  return(thetaTrue)
}



#> Date Setup -----------------------------------------------
GetSimulateData = function (pars_init, pars_out) {
  n = pars_init$n
  m = pars_init$m
  index = seq(n)
  level = seq(m)
  
  # simulate Xi from a multinomial distribution
  px = pars_init$px
  X = rep(0:(m-1), n*px)
  #X = sample(seq(0, m-1), n, replace=T, prob=pars_init$px)
  pxFn = function(x, X){sum(X==x)/length(X)}
  pars_out$px = sapply(0:(m-1), pxFn, X)
  
  # simulate Si_0 given Xi from a Bernoulli distribution
  S0 = rep(NA, n)
  pS01_x = pars_init$pS01_x
  for (i in index) {
    p = pS01_x[X[i]+1]  # X[i]=0:3, starts from 0, needs to +1 to start from 1
    S0[i] = sample(c(0, 1), 1, prob=c(1-p, p))  # Pr[Si_0 = 1] = p
  }
  p_S00Fn = function(x, X, S0){sum(X==x & S0==1)/sum(X==x)}
  pars_out$pS01_x = sapply(0:(m-1), p_S00Fn, X, S0)
  
  # simulate Yi_0 given Si_0 and Xi
  pY01_S00x = pars_init$pY01_S00x
  pY01_S01x = pars_init$pY01_S01x
  Y0 = rep(NA, n)
  for (x in level) {
    for (i in index) {
      if (S0[i]==0 & X[i]==x-1) {
        p=pY01_S00x[x]
        Y0[i] = sample(c(0, 1), 1, prob=c(1-p, p))
      }
      if (S0[i]==1 & X[i]==x-1) {
        p=pY01_S01x[x]
        Y0[i] = sample(c(0, 1), 1, prob=c(1-p, p))
      }
    }
  }
  p_s0Fn = function(x, s0, X, S0, Y0){sum(X==x & S0==s0 & Y0==1)/sum(X==x & S0==s0)}
  pars_out$pY01_S00x = sapply(0:(m-1), p_s0Fn, s0=0, X, S0, Y0)
  pars_out$pY01_S01x = sapply(0:(m-1), p_s0Fn, s0=1, X, S0, Y0)
  
  # simulate Si_1 given Si_0, Yi_0, Xi (use logit)
  pS11_S00Y00x = pars_init$pS11_S00Y00x
  pS11_S00Y01x = pars_init$pS11_S00Y01x
  S1 = rep(NA, n)
  for (x in level) {
    for (i in index) {
      if (S0[i]==1) {S1[i]=1}  # monotonicity
      if (S0[i]==0 & Y0[i]==0 & X[i]==x-1) {
        p = pS11_S00Y00x[x]
        S1[i] = sample(c(0, 1), 1, prob=c(1-p, p))
      }
      if (S0[i]==0 & Y0[i]==1 & X[i]==x-1) {
        p = pS11_S00Y01x[x]
        S1[i] = sample(c(0, 1), 1, prob=c(1-p, p))
      }
    }
  }
  ps1Fn_2 = function(x, y, X, Y0, S0, S1){sum(S0==0 & Y0==y & X==x & S1==1)/sum(S0==0 & Y0==y & X==x)}
  pars_out$pS11_S00Y00x = sapply(0:(m-1), ps1Fn_2, y=0, X, Y0, S0, S1)
  pars_out$pS11_S00Y01x = sapply(0:(m-1), ps1Fn_2, y=1, X, Y0, S0, S1)
  
  # simulate Yi_1 given Si_0, Si_1, Yi_0
  Y1 = rep(NA, n)
  for (i in index) {
    if (S0[i]==0 & S1[i]==0 & Y0[i]==0) {p = pars_init$pY11_S00S10Y00x}
    if (S0[i]==0 & S1[i]==0 & Y0[i]==1) {p = pars_init$pY11_S00S10Y01x}
    if (S0[i]==0 & S1[i]==1 & Y0[i]==0) {p = pars_init$pY11_S00S11Y00x}
    if (S0[i]==0 & S1[i]==1 & Y0[i]==1) {p = pars_init$pY11_S00S11Y01x}
    if (S0[i]==1 & S1[i]==1 & Y0[i]==0) {p = pars_init$pY11_S01S11Y00x}
    if (S0[i]==1 & S1[i]==1 & Y0[i]==1) {p = pars_init$pY11_S01S11Y01x}
    Y1[i] = sample(c(0, 1), 1, prob=c(1-p, p))
  }
  pars_out$pY11_S00S10Y00x = sum(S0==0 & S1==0 & Y0==0 & Y1==1) / sum(S0==0 & S1==0 & Y0==0)
  pars_out$pY11_S00S10Y01x = sum(S0==0 & S1==0 & Y0==1 & Y1==1) / sum(S0==0 & S1==0 & Y0==1)
  pars_out$pY11_S00S11Y00x = sum(S0==0 & S1==1 & Y0==0 & Y1==1) / sum(S0==0 & S1==1 & Y0==0)
  pars_out$pY11_S00S11Y01x = sum(S0==0 & S1==1 & Y0==1 & Y1==1) / sum(S0==0 & S1==1 & Y0==1)
  pars_out$pY11_S01S11Y00x = sum(S0==1 & S1==1 & Y0==0 & Y1==1) / sum(S0==1 & S1==1 & Y0==0)
  pars_out$pY11_S01S11Y01x = sum(S0==1 & S1==1 & Y0==1 & Y1==1) / sum(S0==1 & S1==1 & Y0==1)
  
  # simulate Zi (lastly)
  pz = pars_init$pz[1]
  Z = c()
  for (i in 1:(m*2)) {
    if (i %% 2 == 0) {
      Z = c(Z, rep(1, n/(m*2)))
    } else {
      Z = c(Z, rep(0, n/(m*2)))
    }
  }
  pars_out$pz = c(sum(Z==0)/length(Z), sum(Z==1)/length(Z))
  
  pars_out$m = pars_init$m
  pars_out$n = pars_init$n
  pars_out$b = pars_init$b
  pars_out$r = pars_init$r
  pars_out$beta_true = pars_init$beta_true
  pars_out$beta_init = pars_init$beta_init
  pars_out$optimMethod = pars_init$optimMethod
  pars_out$lambda = pars_init$lambda
  pars_out$pY00_S01x = 1 - pars_out$pY01_S01x
  pars_out$pY00_S00x = 1 - pars_out$pY01_S00x
  pars_out$pS00_x = 1 - pars_out$pS01_x
  pars_out$optimConst = pars_init$optimConst
  
  simulateData = as.data.frame(cbind(index, Z, X, S0, S1, Y0, Y1))
  
  full = list(simulateData=simulateData, pars_out=pars_out)
  return(full)
}

#> Obj fn to est betas -----------------------------------------------
GetOptimBeta_sum = function(beta, GL, GR, lambda, level, optimConst) {
  m = length(level)
  GM = array(NA, dim=c(2, m))
  beta0 = beta[1]
  beta1 = beta[2]
  beta2 = beta[3]
  for (y in 1:2) {
    for (x in level) {
      GM[y,x] = exp(beta0 + beta1*(y-1) + beta2*(x-1)) / (1 + exp(beta0 + beta1*(y-1) + beta2*(x-1)))
    }
  }
  #sum((log((1+GL)/(1-GL)) - log((1+colSums(GM * GR))/(1-colSums(GM * GR))))^2)
  sum((GL - colSums(GM * GR))^2) * optimConst
}


EstimateBeta = function(mydata, beta_init, optimMethod, lambda, optimConst, GetOptimBetaFn) {
  # setup
  Z = mydata$Z
  X = mydata$X
  S0 = mydata$S0
  S1 = mydata$S1
  Y0 = mydata$Y0
  Y1 = mydata$Y1
  n = nrow(mydata)  # number of subjects
  m = length(unique(X))  # number of Xi levels
  level = seq(m)

  # compute N_zsx: N00x, N01x, N10x, N11x. sum=n
  N00x = rep(NA, m)
  N01x = rep(NA, m)
  N10x = rep(NA, m)
  N11x = rep(NA, m)
  for (x in level) {
    N00x[x] = sum(S0==0 & X==x-1)
    N01x[x] = sum(S0==1 & X==x-1)
    N10x[x] = sum(S1==0 & X==x-1)
    N11x[x] = sum(S1==1 & X==x-1)
  }

  # estimate p_zsx: p00x, p01x, p11x. sum=1
  p00x = rep(NA, m)
  p01x = rep(NA, m)
  p11x = rep(NA, m)
  for (i in level) {
    if (N11x[i]/N10x[i] >= N01x[i]/N00x[i]) {
      p00x[i] = N10x[i] / (N10x[i]+N11x[i])
      p11x[i] = N01x[i] / (N00x[i]+N01x[i])
      p01x[i] = 1 - p00x[i] - p11x[i]
    }
    else {
      p00x[i] = (N00x[i]+N10x[i]) / (sum(N00x[i],N01x[i],N10x[i],N11x[i]))
      p11x[i] = (N01x[i]+N11x[i]) / (sum(N00x[i],N01x[i],N10x[i],N11x[i]))
      p01x[i] = 0
    }
  }

  # GL(x): Pr{Si(1) = 1|Si(0) = 0,Xi = x}
  GL = p01x / (1 - p11x)
  
  # GR(x,y): Pr{Yi(0) = y|Si(0) = 0,Xi = x}
  pY01_S00x = rep(NA, m)
  pY01_S01x = rep(NA, m)
  for (x in level) {
    pY01_S00x[x] = sum(Z==0 & S0==0 & Y0==1 & X==x-1) / sum(Z==0 & S0==0 & X==x-1)
    pY01_S01x[x] = sum(Z==0 & S0==1 & Y0==1 & X==x-1) / sum(Z==0 & S0==1 & X==x-1)
  }
  GR = rbind(1-pY01_S00x, pY01_S00x)
  # GR[y, x] = sum(S0==0 & Y0==y-1 & X==x-1) / sum(S0==0 & X==x-1)
  
  res = optim(par=beta_init, fn=GetOptimBetaFn, GL=GL, GR=GR, lambda=lambda, method=optimMethod,
              level=level, optimConst=optimConst, control=list(maxit=1e7))

  # estimate p1=Pr[Y(1)=1 | S(1)=1]
  p1 = sum(Z==1 & S1==1 & Y1==1) / sum(Z==1 & S1==1)
  
  full = list(beta_hat=res$par, obj_loss=res$value, converge=res$converge, 
              GL=GL, GR=GR, p00x=p00x, p11x=p11x, p01x=p01x, p1=p1,
              pY01_S00x=pY01_S00x, pY01_S01x=pY01_S01x)
  return(full)
}


#> Estimate Causal Effect -----------------------------------------------
EstimateEffect = function(mydata, beta_hat, p00x, p11x, p01x, p1, pY01_S00x, pY01_S01x) {
  # setup
  Z = mydata$Z
  X = mydata$X
  S0 = mydata$S0
  S1 = mydata$S1
  Y0 = mydata$Y0
  Y1 = mydata$Y1
  n = nrow(mydata)  # number of subjects
  m = length(unique(X))  # number of Xi levels
  level = seq(m)
  
  # estimate Pr[S_i(0)=0|X_i=x] & Pr[S_i(0)=1|X_i=x] & Pr[S_i(1)=1|X_i=x]
  pS00_x = 1 - p11x
  pS01_x = p11x
  pS11_x = 1 - p00x
  
  # estimate Pr[X_i=x]
  px = rep(NA, m)
  for (x in level) {
    px[x] = sum(X==x-1) / n
  }
  
  # estimate Pr[S_i(1)=1|S_i(0)=0, Y_i(0)=1, X_i=x]
  pS11_S00Y01x = rep(NA, m)
  for (x in level) {
    pS11_S00Y01x[x] = (exp( sum( beta_hat * c(1, 1, x-1 ))) / (1 + exp( sum( beta_hat * c(1, 1, x-1)))))
  }
  
  # estimate Pr[Y_i(0)=1, S_i(1)=1|X_i=x]
  pY01S11_x = pS11_S00Y01x * pY01_S00x * pS00_x + 1 * pY01_S01x * pS01_x
  
  # estimate p0
  p0 = sum(pY01S11_x * px) / sum(pS11_x * px)
  
  # estimate causal effect
  theta = p1-p0
  #cat("Causal effect is:", theta, "\n")
  #print(theta)
  full = c(theta, beta_hat)
  return(full)
}


##### Main Estimation Function #####
Estimation = function(pars_init, GetThetaTrue, GetSimulateData, EstimateBeta, EstimateEffect, GetOptimBeta){
  n = pars_init$n
  b = pars_init$b
  r = pars_init$r
  m = pars_init$m
  beta_true = pars_init$beta_true
  beta_init = pars_init$beta_init
  optimMethod = pars_init$optimMethod
  lambda = pars_init$lambda
  optimConst = pars_init$optimConst
  thetaTrue = GetThetaTrue(pars_init)

  dataTheta = dataBeta0s = dataBeta1s = dataBeta2s = rep(NA, r)
  dataP00xs = dataP11xs = dataP01xs = array(NA, c(r, m))
  imputeTheta = rep(NA, r)
  bootTheta = bootBeta0s = bootBeta1s = bootBeta2s = array(NA, c(b, r))
  bootP00xs = bootP11xs = bootP01xs = vector()
  t0 = proc.time()
  cvr90 = 0
  cvr95 = 0
  for (i in 1:r) {   
    t1 = proc.time()

    # create dataset
    getsimulateData = GetSimulateData(pars_init, list())
    simulateData = getsimulateData$simulateData
    pars_out = getsimulateData$pars_out

    # get estimates based on the dataset
    tmp = EstimateBeta(simulateData, beta_init, optimMethod, lambda, optimConst, GetOptimBeta)
    dataP00xs[i,] = tmp$p00x
    dataP11xs[i,] = tmp$p11x
    dataP01xs[i,] = tmp$p01x
    resData = EstimateEffect(simulateData, tmp$beta_hat, tmp$p00x, tmp$p11x, tmp$p01x, 
                             tmp$p1, tmp$pY01_S00x, tmp$pY01_S01x)
    dataTheta[i] = resData[1]
    dataBeta0s[i] = resData[2]
    dataBeta1s[i] = resData[3]
    dataBeta2s[i] = resData[4]
    
    # impute beta with beta_true
    imputeTheta[i] = EstimateEffect(simulateData, beta_true, tmp$p00x, tmp$p11x, tmp$p01x,
                             tmp$p1, tmp$pY01_S00x, tmp$pY01_S01x)[1]

    # get bootstrap samples - parallel computing
    resBoot = foreach(j=1:b, .combine="rbind") %dopar% {
      #getBootData = GetSimulateData(pars_out, list())
      #bootdata = getBootData$simulateData
      bootdata = simulateData[sample(nrow(simulateData), nrow(simulateData), replace=TRUE), ]
      tmp = EstimateBeta(bootdata, beta_init, optimMethod, lambda, optimConst, GetOptimBeta)
      resBoot = EstimateEffect(bootdata, tmp$beta_hat, tmp$p00x, tmp$p11x, tmp$p01x, 
                             tmp$p1, tmp$pY01_S00x, tmp$pY01_S01x)
      return(resBoot)
    }
    bootTheta[,i] = resBoot[,1]
    bootBeta0s[,i] = resBoot[,2]
    bootBeta1s[,i] = resBoot[,3]
    bootBeta2s[,i] = resBoot[,4]

    sorted = sort(bootTheta[,i])
    lo90 = sorted[b*(0.1/2)]
    up90 = sorted[b*(1 - 0.1/2)]
    lo95 = sorted[b*(0.05/2)]
    up95 = sorted[b*(1 - 0.05/2)]
    
    if (lo90 < thetaTrue & thetaTrue < up90){
      cvr90 = cvr90 + 1  # coverage
    }
    if (lo95 < thetaTrue & thetaTrue < up95){
      cvr95 = cvr95 + 1  # coverage
    }
    
    cat("i=", i, "cvr90=", round(cvr90/i,3), "cvr95=", round(cvr95/i,3), 
        "n=", n, "b=", b, "r=", r, 
        "betaTrue=", beta_true, "thetaTrue=",thetaTrue,
        "meanThetaHat=", round(mean(dataTheta, na.rm=T),3), 
        "bias=", round(mean(dataTheta, na.rm=T)-thetaTrue,3),
        "betaHat=", round(mean(dataBeta0s, na.rm=T),3), 
        round(mean(dataBeta1s, na.rm=T),3), round(mean(dataBeta2s, na.rm=T),3), 
        "time=", round((t1-t0)[3]/60,2), "mins \n")
  }
  full = list(dataTheta=dataTheta, thetaTrue=thetaTrue,
              dataBeta0s=dataBeta0s, dataBeta1s=dataBeta1s, dataBeta2s=dataBeta2s,
              dataP00xs=dataP00xs, dataP11xs=dataP11xs, dataP01xs=dataP01xs,
              imputeTheta=imputeTheta, bootTheta=bootTheta,
              bootBeta0s=bootBeta0s, bootBeta1s=bootBeta1s, bootBeta2s=bootBeta2s,
              bootP00xs=bootP00xs, bootP01xs=bootP01xs, bootP11xs=bootP11xs)
  return(full)
}


#> Analyze Results -----------------------------------------------
Analysis = function(dataTheta, bootTheta, thetaTrue, n, b, r){
  alpha90 = 0.1
  alpha95 = 0.05
  thetaHat = mean(dataTheta)
  thetaBias = mean(dataTheta) - thetaTrue
  thetaVar = var(dataTheta)
  MSEthetaHat = thetaBias^2 + thetaVar
  compare_bias_var = (4*thetaBias)^2 < thetaVar
  cat("thetaTrue:", thetaTrue, "thetaHat", thetaHat, "bootBias:", mean(bootTheta)-thetaTrue, "dataBias:", thetaBias, "Var:", thetaVar, "Compare bias & var:", compare_bias_var, "MSE:", MSEthetaHat, "\n")
  # bootBias and dataBias needs to be equal
  
  res_bootCI90 = BootCIfn(n, b, r, alpha90, bootTheta, dataTheta, thetaTrue)

  res_bootCI95 = BootCIfn(n, b, r, alpha95, bootTheta, dataTheta, thetaTrue)
  
  full = list(thetaHat=thetaHat, thetaBias=thetaBias, 
              thetaVar=thetaVar, MSEthetaHat=MSEthetaHat, 
              res_bootCI90=res_bootCI90, res_bootCI95=res_bootCI95)
  return(full)
}


#> Get Bootstrap CI -------------------------------------------------
BootCIfn = function(n, b, r, alpha, bootTheta, dataTheta, thetaTrue){
  cat("alpha=", alpha, "\n")
  # method 1: percentile interval (intuitive but lack theoretical support)
  upp1 = numeric(r)
  low1 = numeric(r)
  for (j in 1:r){
    sorted = sort(bootTheta[,j])
    low1[j] = sorted[b*(alpha/2)]
    upp1[j] = sorted[b*(1-alpha/2)]
  }
  covrPrCI1 = mean(low1 < thetaTrue & upp1 > thetaTrue)
  avgLenCI1 = mean(upp1 - low1)
  cat("avgLenCI1:", avgLenCI1, ". covrPrCI1:", covrPrCI1, "\n")
  
  # method 2: Larry's Empirical distribution
  low2 = numeric(r)
  upp2 = numeric(r)
  thetaHat = mean(dataTheta)  # change on 2/6/2020
  for(j in 1:r){
    # xbar = mean(bootTheta[,j]);
    IIs=sqrt(n)*(bootTheta[,j]-thetaHat)
    grids=sort(IIs)
    F_hat=rep(0,length(grids))
    for (i in 1:length(grids)){
      F_hat[i]= sum(IIs<=grids[i])/b
    }
    t_u=grids[which.min(F_hat<alpha/2)]
    t_l=grids[which.min(F_hat<1-alpha/2)]
    low2[j]=thetaHat-t_l/sqrt(n)
    upp2[j]=thetaHat-t_u/sqrt(n)
  }
  covrPrCI2 = mean(low2 < thetaTrue & upp2 > thetaTrue) # Approximates to 0.90
  avgLenCI2 = mean(upp2 - low2)
  cat("avgLenCI2:", avgLenCI2, ". coverageCI2:", covrPrCI2, "\n")
  
  # method 3: pivot CI
  upp3 = numeric(r)
  low3 = numeric(r)
  thetaHat = mean(dataTheta)
  for (j in 1:r){
    sorted = sort(bootTheta[,j])
    upp3[j] = 2*thetaHat - sorted[b*(alpha/2)]
    low3[j] = 2*thetaHat - sorted[b*(1-alpha/2)]
  }
  covrPrCI3 = mean(low3<thetaTrue & upp3>thetaTrue)
  avgLenCI3 = mean(upp3 - low3)
  cat("avgLenCI3:", avgLenCI3, ". coverageCI3:", covrPrCI3, "\n")

  # https://ocw.mit.edu/courses/mathematics/18-05-introduction-to-probability-and-statistics-spring-2014/readings/MIT18_05S14_Reading24.pdf
  low4 = numeric(r)
  upp4 = numeric(r)
  for (j in 1:r) {
    delta = bootTheta[,j] - mean(dataTheta)
    sorted_delta = sort(delta)
    low4[j] = mean(dataTheta) - sorted_delta[b*(1-alpha/2)]
    upp4[j] = mean(dataTheta) - sorted_delta[b*(alpha/2)]
  }
  covrPrCI4 = mean(low4 < thetaTrue & upp4 > thetaTrue)
  avgLenCI4 = mean(upp4 - low4)
  cat("avgLenCI4:", avgLenCI4, ". coverageCI4:", covrPrCI4, "\n")
  # alternative (from MIT)
  low5 = numeric(r)
  upp5 = numeric(r)
  for (j in 1:r) {
    delta = bootTheta[,j] - mean(dataTheta)
    low5[j] = mean(dataTheta) - quantile(delta, 1-alpha/2)
    upp5[j] = mean(dataTheta) - quantile(delta, alpha/2)
  }
  covrPrCI5 = mean(low5 < thetaTrue & upp5 > thetaTrue)
  avgLenCI5 = mean(upp5 - low5)
  cat("avgLenCI5:", avgLenCI5, ". coverageCI5:", covrPrCI5, "\n")

  full = list(avgLenCI1=avgLenCI1, covrPrCI1=covrPrCI1, low1=low1, upp1=upp1,
              avgLenCI2=avgLenCI2, covrPrCI2=covrPrCI2, low2=low2, upp2=upp2,
              avgLenCI3=avgLenCI3, covrPrCI3=covrPrCI3, low3=low3, upp3=upp3,
              avgLenCI4=avgLenCI4, covrPrCI4=covrPrCI4, low4=low4, upp4=upp4)
  return(full)
}


