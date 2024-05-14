#===============================================================================

# File name: sim.code.R

# Purpose: Code for each simulation studies

# Author: Ke Wan

# R version: R-4.1.2

# Input: rep: Times to repeat

# Output: a list including --->  sim.res.can: true and estimated HTE for each method;
#                                sim.res: summary of the simulation result, 
#                                        including: RMSE, rbias (Rbias0), variance, spearman's rank correlation (cor.t),
#                                        computation times (time.t);
#                                name: Information about simulation settings;
#                                data: Train data; true_mu_tau.train: The mu, tau & pscore used to create train data;
#                                data.test: Test data; true_mu_tau.test: The mu, tau & pscore used to create test data;
#                                pscore: Estimated propensity score.

# Require packages: bartCause, grf, causalLearning, pre, grpreg.

#===============================================================================

sim <- function(rep){
  
  # Probability for binary variables  
  
  set.seed(rep)
  pr <- runif(p/2,0.25,0.75) 
  
  # Data generation
  
  set.seed(rep + 1)
  dat <- data_gen(n=n,p=p,alpha=alpha,cor=cor,pr=pr,model=model,type = type)
  
  set.seed(rep + 10)
  dat.test <- data_gen(n=n,p=p,alpha=alpha,cor=cor,pr=pr,model=model,type = type)
  
  data <- dat$data
  data.test <- dat.test$data
  
  # BART for treatment effec estimation (Hill et al. 2011) #################################
  
  set.seed(rep + 100)
  start_t.bart <- proc.time()
  bart.fit <- bartCause::bartc(data[,1],data[,2],data[,-c(1:2),drop = FALSE],keepTrees = TRUE,seed = rep + 100)
  end_t.bart <- proc.time()
  pred.bart <- predict(bart.fit,data.test,type = "ite")
  pred.bart <- apply(pred.bart,2,mean)
  rmse.bart <- sqrt(mean((dat.test$tau - pred.bart)^2))
  bias.bart0 <- mean((dat.test$tau - pred.bart)/dat.test$tau)
  bias.bart1 <- median((dat.test$tau - pred.bart)/dat.test$tau)
  var.bart <- var(pred.bart)
  cor.bart <- cor(dat.test$tau,pred.bart, method = "spearman")
  time.bart <- (end_t.bart - start_t.bart)[3]
  
  # ** Calculate the propensity score (only for observational study setting) ###############
  
  if(type == "rct"){
    
    pscore <- rep(0.5,nrow(data))
    
  }else{
    
    pscore <- bart.fit$p.score
    
  }
  
  # Causal forest (Athey et al. 2019) #######################################################
  
  set.seed(rep + 100)
  start_t.cf <- proc.time()
  cf.fit <- grf::causal_forest(X = as.matrix(data[,-c(1:2),drop = FALSE]),Y = data[,1],W = data[,2],tune.parameters = "all")
  end_t.cf <- proc.time()
  pred.cf <- predict(cf.fit,as.matrix(data.test[,-c(1:2),drop = FALSE]))$predictions
  rmse.cf <- sqrt(mean((dat.test$tau - pred.cf)^2))
  bias.cf0 <- mean((dat.test$tau - pred.cf)/dat.test$tau)
  bias.cf1 <- median((dat.test$tau - pred.cf)/dat.test$tau)
  var.cf <- var(pred.cf)
  cor.cf <- cor(dat.test$tau,pred.cf, method = "spearman")
  time.cf <- (end_t.cf - start_t.cf)[3]
  vip.cf <- grf::variable_importance(cf.fit)
  
  # Bagging Causal MARS (Powers et al. 2018) ################################################
  set.seed(rep + 100)
  strata <- causalLearning::stratify(pscore = pscore,tx = data[,2])$stratum
  start_t.bcm <- proc.time()
  bcm.fit <- try(causalLearning::bagged.causalMARS(x = as.matrix(data[,-c(1:2),drop = FALSE]),y = data[,1],tx = data[,2],propensity = TRUE,stratum = strata),silent = TRUE)
  end_t.bcm <- proc.time()
  
  pred.bcm <- try(predict(bcm.fit,newx = as.matrix(data.test[,-c(1:2),drop = FALSE])),silent = TRUE)
  
  if (class(bcm.fit)!="try-error" & class(pred.bcm)!="try-error"){
    rmse.bcm <- sqrt(mean((dat.test$tau - pred.bcm)^2))
    bias.bcm0 <- mean((dat.test$tau - pred.bcm)/dat.test$tau)
    bias.bcm1 <- median((dat.test$tau - pred.bcm)/dat.test$tau)
    var.bcm <- var(pred.bcm)
    cor.bcm <- cor(dat.test$tau,pred.bcm, method = "spearman")
    time.bcm <- (end_t.bcm - start_t.bcm)[3]
  }else{
    pred.bcm <- rep(NA,n)
    rmse.bcm <- NA
    bias.bcm0 <- NA
    bias.bcm1 <- NA
    var.bcm <- NA
    cor.bcm <- NA
    time.bcm <- NA
  }  
  
  # PTO forest (Powers et al. 2018) #########################################################
  set.seed(rep + 100)
  start_t.pto <- proc.time()
  pto.fit <- try(causalLearning::PTOforest(x = as.matrix(data[,-c(1:2),drop = FALSE]),y = data[,1],tx = data[,2],pscore = pscore),silent = TRUE)
  end_t.pto <- proc.time()
  
  if (class(pto.fit)!="try-error"){
    pred.pto <- predict(pto.fit,newx = as.matrix(data.test[,-c(1:2),drop = FALSE]))
    rmse.pto <- sqrt(mean((dat.test$tau - pred.pto)^2))
    bias.pto0 <- mean((dat.test$tau - pred.pto)/dat.test$tau)
    bias.pto1 <- median((dat.test$tau - pred.pto)/dat.test$tau)
    var.pto <- var(pred.pto)
    cor.pto <- cor(dat.test$tau,pred.pto, method = "spearman")
    time.pto <- (end_t.pto - start_t.pto)[3]
  }else{
    pred.pto <- rep(NA,n)
    rmse.pto <- NA
    bias.pto0 <- NA
    bias.pto1 <- NA
    var.pto <- NA
    cor.pto <- NA
    time.pto <- NA
  }
  
  # RuleFit method (Friedman and Popesecu,2010) for transformed outcome #####################
  set.seed(rep + 100)
  trans_y <- data[,2]*(data[,1]/pscore) - (1-data[,2])*(data[,1]/(1-pscore))
  trans_data <- data.frame(y = trans_y,data[,-c(1:2)])
  start_t.rulefit <- proc.time()
  rule.fit <-  pre::pre(y ~ ., data = trans_data)
  end_t.rulefit <- proc.time()
  pred.rulefit <- predict(rule.fit,data.test)
  rmse.rulefit <- sqrt(mean((dat.test$tau - pred.rulefit)^2))
  bias.rulefit0 <- mean((dat.test$tau - pred.rulefit)/dat.test$tau)
  bias.rulefit1 <- median((dat.test$tau - pred.rulefit)/dat.test$tau)
  var.rulefit <- var(pred.rulefit)
  cor.rulefit <- cor(dat.test$tau,pred.rulefit, method = "spearman")
  time.rulefit <- (end_t.rulefit - start_t.rulefit)[3]
  
  #-----------------------------------------------------------------------------
  if(is.na(cor.rulefit)){
    imp.rulefit <- rep(NA,ncol(data)-2)
  }else{
    imp.rulefit <- rep(0,ncol(data)-2)
    names(imp.rulefit) <- colnames(data[,-c(1:2)])
    vip.rulefit <- pre::importance(rule.fit)$varimps
    for(vv in vip.rulefit$varname){
      imp.rulefit[vv] <- vip.rulefit$imp[vip.rulefit$varname == vv]
    }
  }
  vip.rule_fit <- imp.rulefit
  #-----------------------------------------------------------------------------
  
  # Proposed method rule generation #########################################################
  set.seed(rep + 100)
  start_t.prop0 <- proc.time()
  mat_rules <- rules.generate(data = data,meandepth = 3,ntrees = 333,shrinkage = 0.01,pscore = pscore)
  grp.mat <- rules.ensemble(data,rules = mat_rules$rules)
  grp.mat.test <- mat_grp_test(data.test,rules = mat_rules$rules,grp.mat$scale)
  end_t.prop0 <- proc.time()
  time.prop0 <- (end_t.prop0 - start_t.prop0)[3]
  
  # Proposed method (group lasso) with regularization parameter selected by cross-validation #
  set.seed(rep + 100)
  start_t.propcv <- proc.time()
  cv_grpreg <- grpreg::cv.grpreg(X = as.matrix(grp.mat$x_org),y = grp.mat$y,group = grp.mat$index,nfold = 5,returnX = TRUE)
  end_t.propcv <- proc.time()
  t1.grpreg.cv <- predict(cv_grpreg,as.matrix(grp.mat.test$x_trt))
  t0.grpreg.cv <- predict(cv_grpreg,as.matrix(grp.mat.test$x_con))
  tau.grpreg.cv <- t1.grpreg.cv - t0.grpreg.cv
  rmse.cv <- sqrt(mean((dat.test$tau - tau.grpreg.cv)^2))
  bias.cv0 <- mean((dat.test$tau - tau.grpreg.cv)/dat.test$tau)
  bias.cv1 <- median((dat.test$tau - tau.grpreg.cv)/dat.test$tau)
  var.cv <- var(tau.grpreg.cv)
  cor.cv <- cor(dat.test$tau,tau.grpreg.cv, method = "spearman")
  time.propcv <- (end_t.propcv - start_t.propcv)[3] + time.prop0
  
  # variable importance -----------------------------------------------------------------
  fit_cv <- cv_grpreg$fit
  lambda.id <- which(cv_grpreg$lambda == cv_grpreg$lambda.min)
  beta <- fit_cv$beta[,lambda.id]
  vip.cv <- VIP.prop(data,grp.mat,rules = mat_rules$rules,beta)[[2]]
  
  # Proposed method (group lasso) with regularization parameter selected by bic
  set.seed(rep + 100)
  start_t.propbic <- proc.time()
  bic_grpreg <- grpreg::grpreg(X = as.matrix(grp.mat$x_org),y = grp.mat$y,group = grp.mat$index)
  end_t.propbic <- proc.time()
  lambda <- grpreg::select(bic_grpreg,criterion="BIC")$lambda 
  t1.grpreg.bic <- predict(bic_grpreg,as.matrix(grp.mat.test$x_trt),lambda = lambda)
  t0.grpreg.bic <- predict(bic_grpreg,as.matrix(grp.mat.test$x_con),lambda = lambda)
  tau.grpreg.bic <- t1.grpreg.bic - t0.grpreg.bic
  rmse.bic <- sqrt(mean((dat.test$tau - tau.grpreg.bic)^2))
  bias.bic0 <- mean((dat.test$tau - tau.grpreg.bic)/dat.test$tau)
  bias.bic1 <- median((dat.test$tau - tau.grpreg.bic)/dat.test$tau)
  var.bic <- var(tau.grpreg.bic)
  cor.bic <- cor(dat.test$tau,tau.grpreg.bic, method = "spearman")
  time.propbic <- (end_t.propbic - start_t.propbic)[3] + time.prop0
  
  # variable importance ---------------------------------------------------------------------------------
  fit_bic <- bic_grpreg
  lambda.id <- which(bic_grpreg$lambda == lambda)
  beta <- fit_bic$beta[,lambda.id]
  vip.bic <- VIP.prop(data,grp.mat,rules = mat_rules$rules,beta)[[2]]
  
  # ** Calculate adaptive weights (Huang et al. 2009) ##################
  
  start_t.adw <- proc.time()
  mod <- cv_grpreg$fit
  lambda <- cv_grpreg$lambda.min
  beta0 <- mod$beta[,which(mod$lambda == lambda)][-1]
  
  beta0 <- beta0*mod$XG$scale
  beta0[which(is.nan(beta0))] <- 0
  
  grp_weights <- numeric(length(unique(grp.mat$index)))
  for (i in unique(grp.mat$index)) {
    grp_weights[i] <-  1 / sqrt(sum(beta0[which(grp.mat$index == i)]^2))
  }
  end_t.adw <- proc.time()
  time.adw <- (end_t.adw - start_t.adw)[3] + time.propcv
  
  ######################################################################
  
  # Proposed method (adaptive group lasso) with regularization parameter selected by cross-validation
  set.seed(rep + 100)
  start_t.propadcv <- proc.time()
  ad_grpreg.cv <- grpreg::cv.grpreg(X = as.matrix(grp.mat$x_org),y = grp.mat$y,group = grp.mat$index,group.multiplier = grp_weights,nfold = 5,lambda.min = 0.0001,nlambda = 200)
  end_t.propadcv <- proc.time() 
  t1.adgrpreg.cv <- predict(ad_grpreg.cv,as.matrix(grp.mat.test$x_trt))
  t0.adgrpreg.cv <- predict(ad_grpreg.cv,as.matrix(grp.mat.test$x_con))
  tau.adgrpreg.cv <- t1.adgrpreg.cv - t0.adgrpreg.cv
  rmse.ad.cv <- sqrt(mean((dat.test$tau - tau.adgrpreg.cv)^2))
  bias.ad.cv0 <- mean((dat.test$tau - tau.adgrpreg.cv)/dat.test$tau)
  bias.ad.cv1 <- median((dat.test$tau - tau.adgrpreg.cv)/dat.test$tau)
  var.ad.cv <- var(tau.adgrpreg.cv)
  cor.ad.cv <- cor(dat.test$tau,tau.adgrpreg.cv, method = "spearman")
  time.propadcv <- (end_t.propadcv - start_t.propadcv)[3] + time.adw
  
  # variable importance ----------------------------------------------------------------
  fit_adcv <- ad_grpreg.cv$fit
  lambda.id <- which(ad_grpreg.cv$lambda == ad_grpreg.cv$lambda.min)
  beta <- fit_adcv$beta[,lambda.id]
  vip.adcv <- VIP.prop(data,grp.mat,rules = mat_rules$rules,beta)[[2]]
  
  # Proposed method (adaptive group lasso) with regularization parameter selected by bic
  set.seed(rep + 100)
  start_t.propadbic <- proc.time()
  ad_grpreg.bic <- grpreg::grpreg(X = as.matrix(grp.mat$x_org),y = grp.mat$y,group = grp.mat$index,group.multiplier = grp_weights,lambda.min = 0.0001,nlambda = 200)
  end_t.propadbic <- proc.time()
  lambda <- grpreg::select(ad_grpreg.bic,criterion="BIC")$lambda  
  t1.adgrpreg.bic <- predict(ad_grpreg.bic,as.matrix(grp.mat.test$x_trt),lambda = lambda)
  t0.adgrpreg.bic <- predict(ad_grpreg.bic,as.matrix(grp.mat.test$x_con),lambda = lambda)
  tau.adgrpreg.bic <- t1.adgrpreg.bic - t0.adgrpreg.bic
  rmse.ad.bic <- sqrt(mean((dat.test$tau - tau.adgrpreg.bic)^2))
  bias.ad.bic0 <- mean((dat.test$tau - tau.adgrpreg.bic)/dat.test$tau)
  bias.ad.bic1 <- median((dat.test$tau - tau.adgrpreg.bic)/dat.test$tau)
  var.ad.bic <- var(tau.adgrpreg.bic)
  cor.ad.bic <- cor(dat.test$tau,tau.adgrpreg.bic, method = "spearman")
  time.propadbic <- (end_t.propadbic - start_t.propadbic)[3] + time.adw
  
  # variable importance -----------------------------------------------------------
  fit_adbic <- ad_grpreg.bic
  lambda.id <- which(ad_grpreg.bic$lambda == lambda)
  beta <- fit_adbic$beta[,lambda.id]
  vip.adbic <- VIP.prop(data,grp.mat,rules = mat_rules$rules,beta)[[2]]
  
  # save the true treatment effect and each predict treatment effect. 
  
  sim.res.can <- data.frame(true = dat.test$tau,bart = pred.bart,bcm = pred.bcm,pto = pred.pto,cf = pred.cf,rulefit = pred.rulefit,prop.cv = tau.grpreg.cv,prop.bic = tau.grpreg.bic,prop.ad.cv = tau.adgrpreg.cv, prop.ad.bic = tau.adgrpreg.bic)
  
  # simulation result------------------------------------------------------------
  
  # 1 summary of simulation
  
  RMSE <- c(rmse.bart,rmse.cf,rmse.bcm,rmse.pto,rmse.rulefit,rmse.cv,rmse.bic,rmse.ad.cv,rmse.ad.bic)
  Rbias0 <- c(bias.bart0,bias.cf0,bias.bcm0,bias.pto0,bias.rulefit0,bias.cv0,bias.bic0,bias.ad.cv0,bias.ad.bic0)
  Rbias1 <- c(bias.bart1,bias.cf1,bias.bcm1,bias.pto1,bias.rulefit1,bias.cv1,bias.bic1,bias.ad.cv1,bias.ad.bic1)
  Variance <- c(var.bart,var.cf,var.bcm,var.pto,var.rulefit,var.cv,var.bic,var.ad.cv,var.ad.bic)
  cor.t <- c(cor.bart,cor.cf,cor.bcm,cor.pto,cor.rulefit,cor.cv,cor.bic,cor.ad.cv,cor.ad.bic)
  time.t <- c(time.bart,time.cf,time.bcm,time.pto,time.rulefit,time.propcv,time.propbic,time.propadcv,time.propadbic)
  sim.res <- rbind(RMSE,Rbias0,Rbias1,Variance,cor.t,time.t)
  colnames(sim.res) <- c("bart","cf","bcm","pto","rulefit","prop.gl.cv","prop.gl.bic","prop.adgl.cv","prop.adgl.bic")
  
  # 2 variable importance of cf, rule_fit and proposed method 
  
  vip <- cbind(vip.cf,vip.rule_fit,vip.cv,vip.bic,vip.adcv,vip.adbic)
  colnames(vip) <- c("cf","rulefit","prop.gl.cv","prop.gl.bic","prop.adgl.cv","prop.adgl.bic")
  
  # 3 TRUE tau, mu, treatment indicator(z)
  
  true_mu_tau.train <- data.frame(y = dat$y, z = data$z, pscore = dat$p.score, tau = dat$tau, mu = dat$mu)
  true_mu_tau.test <- data.frame(y = dat.test$y, z = data.test$z, pscore = dat.test$p.score, tau = dat.test$tau, mu = dat.test$mu)
  
  ###########################################################################################################################   
  name <- paste0(type,"_",model[2],"_",model[1],"_n_",n,"_p_",p,"_cor_",cor,"_alpha_",alpha,"_rep_",rep,".csv")
  
  res <- list(sim.res.can,sim.res,name,vip,data,true_mu_tau.train,data.test,true_mu_tau.test,pscore)
  
  return(res)
}



