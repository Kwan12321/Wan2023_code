#===============================================================================

# File name: variable_importance.R

# Purpose: Calculate Rule importance & Variable importance

# Author: Ke Wan

# R version: R-4.1.2

# Input: rep: data: Data frame; grp.mat: Object from function "rules.ensemble" (rule_ensemble.R);
#        rules: Created rules; beta: Coefficient of final model. 

# Output: a list including ---> rip: Rule importance; vip: Variables importance

# Require packages: stringr

#===============================================================================

VIP.prop <- function(data,grp.mat,rules,beta){
  
  # Remove intercept
  beta.t <- beta[-1]
  
  # number of coefficients
  coef.num <- length(beta.t)
  con.id <- (1:coef.num)[(1:coef.num)%%2 == 1]
  trt.id <- (1:coef.num)[(1:coef.num)%%2 == 0]
  
  # Rule importance ============================================================
  
  # rule coefficient values
  con.coef <- beta.t[con.id]
  trt.coef <- beta.t[trt.id]
  coef.diff <- trt.coef - con.coef
  names(coef.diff) <- c(colnames(data)[-c(1:2)],rules)
  rule_vars <- rule.predict(data,rules)
  coefr.sd <- apply(rule_vars,2,sd)
  coefv.sd <- grp.mat$scale*0.4
  rip <- abs(coef.diff*c(coefv.sd,coefr.sd))
  rip.res <- data.frame(rip = rip,beta = coef.diff)
  
  # Variable importance=========================================================
  
  #store the variables importance and corresponding variables names  
  vip <- c()
  name.vip <- c()
  
  #calculate the variable importance for each variables
  for(var.id in 1:length(colnames(data)[-c(1:2)])){
    
    var <- colnames(data[,-c(1:2)])[var.id]
    
    #importance of linear terms
    imp.linears <- rip[names(rip)==var]
    
    #Catch the rules that each variable constructs.
    rules.t <- paste0("(",var,"<","|",var,">=",")")
    var_rules <- grep(rules.t,names(rip))
    imp.rules <- rip[var_rules]
    
    #the weight for the importance of each rules 
    rules.weight <- stringr::str_count(names(imp.rules),"&") + 1
    
    #the weighted rule importance
    imp.rules <- sum(imp.rules/rules.weight)
    
    #variable importance
    vip <- c(vip,(imp.linears + imp.rules))
    name.vip <- c(name.vip,var)
  }
  res <- list(rip = rip.res,vip=vip)
  return(res)
}