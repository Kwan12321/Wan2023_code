#===============================================================================

# File name: proposed_method_prediction.R

# Purpose: The prediction function for proposed method 

# Author: Ke Wan

# R version: R-4.1.2

# Input: object: The object from function "prop"; 
#        newdata: Input newdata. 
#        treatment: estimated outcome for treatment group,
#        control: estimated outcome for control group,
#        tau: estimated HTEs

# Require packages: none.

#===============================================================================

prop_predict <- function(object, # The object from function "prop"
                         newdata # Input newdata
                ){
                
                # transform the data for proposed method
                grp.mat.test <- mat_grp_test(data = newdata, rules = object$rules, scale = object$scale_x) # subfunction 1
                treatment.est <- as.matrix(grp.mat.test$x_trt)%*%object$beta[-1] + object$beta[1]
                control.est <- as.matrix(grp.mat.test$x_con)%*%object$beta[-1] + object$beta[1]
                tau <- treatment.est - control.est
                
                return(list(
                  treatment.est = treatment.est,
                  control.est = control.est,
                  tau = tau
                  
                ))
}

##################################################################################################
#Sub_function 1: Transform the test data for proposed method
##################################################################################################
mat_grp_test <- function(data,rules,scale){
  
  # Extract the variable names---------------------------------------------------
  
  y_names <- colnames(data)[1]			# Outcome name
  g_names <- colnames(data)[2]			# Treatment indicator name
  x_names <- colnames(data)[-c(1:2)]		# Covariates names
  
  # Calculate the transformed outcome and create the transformed outcome data----
  
  y <- data[,y_names]						# Outcome
  g <- data[,g_names]						# Treatment indicator
  x <- data[,x_names]						# Covariates
  
  # Scale the covariates
  
  x <- scale(x,center = FALSE, scale = scale)
  data.t <- data.frame(y=data[,y_names],g=data[,g_names],x)
  colnames(data.t) <- c(y_names,g_names,x_names)
  
  # Predict rules
  
  rulevars <- rule.predict(data,rules)
  colnames(rulevars) <- rules 
  
  # Rebuild the data frame
  
  comb.mat <- cbind(data.t,rulevars)
  
  # Transform the data for group lasso
  
  grp.mat <- mat_grp(comb.mat)
  
  grp.mat$scale <- scale
  
  return(grp.mat)
  
}
