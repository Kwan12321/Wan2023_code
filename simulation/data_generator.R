#===============================================================================

# File name: data_generator.R

# Purpose: Create simulation datasets

# Author: Ke Wan

# R version: R-4.1.2

# Input: n: sample size; p: number of variables; alpha: shape parameter for skewed normal distribution;
#        cor: correlation between continuous variables; pr: probability for generating binary variables;
#        model: model id (First for HTE; Second for main effect); type: dataset type ("rct"/"obv") 

# Output: a list including --> data: simulation dataset; tau: HTE; mu: main effect; 
#                              y: average outcome for each subjects; 
#                              p.score: true propensity score for treatment indicator g

# Require packages: sn, MASS

#===============================================================================

# 1: prognostic effect model (M1 - M3)------------------------------------------

mu <- function(X,model = "model_1"){
  
  if(model == "model_1"){
    
    # mean = 0/sd = 2.6
    mu <- X[,1] - (X[,2]+X[,3]) + 2*X[,4] + X[,5] + abs(X[,6] - X[,7]) - (2*X[,8] + X[,9] + X[,10])
  }
  else if(model == "model_2"){
    
    # mean = 0/sd = 2.6
    mu <- 2*X[,1] - 2*I(X[,2] > 0)*I(X[,3] > 0) - (X[,4] + X[,5]) + abs(X[,6] - X[,7]) + X[,8]*X[,9]*X[,10]
  }
  else if(model == "model_3"){
    
    # mean = 0/sd = 2.6
    mu <- X[,1] + X[,3] - X[,2]*X[,4] - 2.5*sin(X[,5]^3 + X[,6])  + sin(exp(0.5*I(X[,7]>0)*I(X[,8]>0))) + 2*X[,9]*X[,10]
  }
  mu
}

# 2: treatment effect model (T1 - T4)-------------------------------------------

tau <- function(X,model = "model_1"){
  if(model == "model_1"){
    # mean = 3/ sd = 0
    tau <- rep(3,nrow(X)) 
  }
  else if(model == "model_2"){
    # mean = 3/ sd = 3.8
    tau <- 3*X[,1] + 2*X[,2] + X[,3] + X[,4] + X[,5] + X[,6] - X[,7] + X[,8] - X[,9] + X[,10] 
  }
  else if(model == "model_3"){
    # mean = 3/ sd = 3.8
    tau <- 3*I(X[,1]>-1)*I(X[,2] >0) + 3*I(X[,3]< 1)*I(X[,4]>0) + 3*X[,5] + X[,6] + X[,7] + X[,8]*X[,9]*X[,10]
  }
  else if(model == "model_4"){
    # mean = 3/ sd = 3.8
    tau <- 3*X[,1] + 3*X[,2] + sin(3*X[,3]^3 + 3*X[,4] + exp(2*I(X[,5] > 0)*I(X[,6] > 0))) + X[,7]^2 + X[,8] + X[,9]*X[,10]
  }
  tau
}

# 3: data generation function---------------------------------------------------

data_gen <- function(n,p,alpha=0,cor=0,pr=rep(0.5,p/2),model = c("model_1","model_1"),type = "rct"){
  
  # Explanatory variables (continuous)------------------------------------------
  
  if(alpha!=0){
    
    # Skewed normal distribution with scale parameter 1; shape parameter 1
    X <- matrix(sn::rsn(n*p,xi = 0, omega = 1, alpha = alpha),n,p)
  }else{
    
    # Standard normal (cor = 0)/correlation structure (cor = 0.3)
    if(cor == 0){
      sigma_dd <- diag(p)
    }else{
      sigma_dd <- diag(p)
      sigma_dd[upper.tri(sigma_dd)] <- cor
      sigma_dd[lower.tri(sigma_dd)] <- cor
    }		
    mu_dd <- rep(0,p)
    X <- MASS::mvrnorm(n,mu_dd,sigma_dd)  
  }
  
  # Explanatory variables (binary) ---------------------------------------------
  id <- 1
  for (j in 1:p){
    if(j%%2==0){
      X[,j] <- rbinom(n,1,pr[id])
      id <- id + 1
    }	
  }	
  
  # Treatment Effect -----------------------------------------------------------
  
  tau.t <- tau(X,model = model[1])
  print(paste0("tau_",sd(tau.t),"_",mean(tau.t)))
  
  # Prognostic Effect -----------------------------------------------------------
  
  mu.t <- mu(X, model = model[2])
  print(paste0("mu_",sd(mu.t),"_",mean(mu.t)))
  
  # Treatment indicators -------------------------------------------------------
  
  if(type == "rct"){
    p.score <- rep(0.5,n)
    g <- rbinom(n,1,p.score)
  }else{
    p.score <- exp(mu.t - 0.5*tau.t)/(1 + exp(mu.t - 0.5*tau.t))
    g <- rbinom(n,1,p.score)
  }
  
  # Response variables ---------------------------------------------------------
  
  Y.t <- mu.t + (g - 0.5)*tau.t
  Y <- rnorm(n,Y.t,1)
  
  data <- data.frame(y = Y,z=as.numeric(g),X)
  
  return(list(data = data,tau=tau.t,mu = mu.t,y = Y.t,p.score = p.score))
}
