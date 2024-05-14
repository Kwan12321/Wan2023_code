#===============================================================================

# File name: Data_transform.R

# Purpose: Reconstruct the data frame for group lasso/adaptive group lasso

# Author: Ke Wan

# R version: R-4.1.2

# Input: data: Data frame; rules: generated rules; 
#        rulevars: rules matrix (transform rules into binary values(0,1))
#        wins: Winsorized size (0.025 is recommeded in Freidman and Popescu, 2008)

# Output: a list including ---> y: outcomes; 
#                               g: treatment indicator;
#                               x_org: Covariates of the reconstruct data frame;
#                               x_con: Covariates of the reconstruct data frame for control group (0 ifelse);
#                               x_trt: Covariates of the reconstruct data frame for treatment group (0 ifelse);
#                               index: group id
#

# Require packages: none

#===============================================================================

rules.ensemble <-  function(data,						# Input data frame
                            rules,						# The rules set
                            rulevars = NULL,			# The rules matrix
                            wins = 0.025				# The winsorized size
){
  
  # Summary of rules and data 
  
  n <- nrow(data)				# The sample size for input data
  p <- ncol(data)				# The number of variables for input data
  num_rules <- length(rules)	# The number of input rules 
  
  # Extract the variable names---------------------------------------------------
  
  y_names <- colnames(data)[1]			# Outcome name
  g_names <- colnames(data)[2]			# Treatment indicator name
  x_names <- colnames(data)[-c(1:2)]		# Covariates names
  
  # Calculate the transformed outcome and create the transformed outcome data----
  
  y <- data[,y_names]						# Outcome
  g <- data[,g_names]						# Treatment indicator
  x <- data[,x_names]						# Covariates
  
  # Get the rules matrix
  
  if (is.null(rulevars)){
    rulevars <- rule.predict(data,rules)
    colnames(rulevars) <- rules 
  }else{
    rulevars <- rulevars
  }
  
  # Get the winsorized covariates(Freideman and Popescu, 2008)------------------
  
  data.win <- winsorize(data,wins) # The winsorized function (Subfunction 1)
  data.win <- data.win[[1]]        # The data with winsorized covariate
  
  # Normalized the covariates---------------------------------------------------
  
  data.norm <- normalize(data.win)
  scale_x <- data.norm$scale		 # Standard deviation for each variables 
  data.norm <- data.norm$data		 # The data with winsorized and normalized covariate
  
  # Combine the rules and linear terms------------------------------------------
  
  comb.mat <- cbind(data.norm,rulevars)
  
  # Transform the data for group lasso
  
  grp.mat <- mat_grp(comb.mat)
  
  # Summary the data
  
  grp.mat$scale <- scale_x
  
  return(grp.mat)
  
}							
###############################################################################


#-----------------------------------------------------------------
#Sub_function 1: Winsorize function
#-----------------------------------------------------------------	
winsorize <- function(data,wins){
  
  # Summary of the input data
  
  y_names <- colnames(data)[1]			# Outcome name
  g_names <- colnames(data)[2]			# Treatment indicator name
  x_names <- colnames(data)[-c(1:2)]		# Covariates names
  data.all <- data
  data <- data[,-c(1:2),drop = FALSE]										
  
  # Winsorize each variables
  data <- sapply(1:ncol(data),function(j,data,wins){
    if(length(unique(data[,j]))> 1){		
      lb <- quantile(data[,j],prob = wins)
      ub <- quantile(data[,j],prob = 1 - wins)
      data[which(data[,j] < lb),j] <- lb  # lower bound 
      data[which(data[,j] > ub),j] <- ub  # upper bound
      data.w <- data[,j]
    }else{									# The variable with binary values do not perform the wisorize 
      data.w <- data[,j]
    }		
    return(data.w)
  },data,wins)
  
  # Get the bound of each variable
  
  lb.point <- apply(data,2,min)
  ub.point <- apply(data,2,max)
  win_cut <- rbind(lb = lb.point,ub = ub.point) 
  
  # Rebuild the data frame
  
  data <- data.frame(y=data.all[,y_names],g=data.all[,g_names],data)
  colnames(data) <- c(y_names,g_names,x_names)
  
  return(list(data = data,win_cut = win_cut))
}						   

#-----------------------------------------------------------------
#Sub_function 2: Normalize function
#-----------------------------------------------------------------	
normalize <- function(data){
  
  # Summary of the input data
  
  y_names <- colnames(data)[1]			# Outcome name
  g_names <- colnames(data)[2]			# Treatment indicator name
  x_names <- colnames(data)[-c(1:2)]		# Covariates names
  data.all <- data
  data <- data[,-c(1:2),drop = FALSE]										
  
  #normalize the linear term
  
  scale_x <- c()
  for(j in 1:ncol(data)){
    scale_x <- c(scale_x,sd(data[,j]) / 0.4) 	
  }			
  
  data <- scale(data,center = FALSE, scale = scale_x)
  
  # Rebuild the data frame
  
  data <- data.frame(y=data.all[,y_names],g=data.all[,g_names],data)
  colnames(data) <- c(y_names,g_names,x_names)
  
  return(list(data = data,scale = scale_x))
}			

#------------------------------------------------------------------
#Sub_function 3: Transform the matrix for group lasso 
#------------------------------------------------------------------
mat_grp <- function(data){
  
  # Summary of the input data
  
  data.all <- data
  names <- colnames(data.all)
  data <- data[,-c(1:2),drop = FALSE]
  
  # Create the grouped variables
  
  index <- rep(1:ncol(data),each = 2)
  data.grp <- data[,index] 
  
  # Rename the grouped variables
  
  names_var <- names[-c(1:2)]
  colnames(data.grp)[1:ncol(data.grp)%%2 == 1] <- paste(names_var,"_c")  # Set odd numbered variables as control group variables ("XXX_c")
  colnames(data.grp)[1:ncol(data.grp)%%2 == 0] <- paste(names_var,"_t")  # Set even numbered variables as treatment group variables ("XXX_t")
  
  # Reconstruct the data frame for applying group lasso
  
  t <- cbind(1 - data.all[,2],data.all[,2])
  data.org <- data.grp*t
  
  # Data with all the treatment indicator equal to 0 
  
  t0 <- cbind(rep(1,nrow(data)),rep(0,nrow(data)))
  data.con <- data.grp*t0
  
  # Data with all the treatment indicator equal to 1
  
  t1 <- cbind(rep(0,nrow(data)),rep(1,nrow(data)))
  data.trt <- data.grp*t1
  
  return(list(y = data.all[,1],g = data.all[,2],x_org = data.org,x_con = data.con,x_trt = data.trt,index = index))
}

#-----------------------------------------------------------------------
#Sub_function 4: Binary matrix for test data
#------------------------------------------------------------------------
mat_grp_test <- function(data,rules,scale){
  
  # Summary of rules and data 
  
  n <- nrow(data)				# The sample size for input data
  p <- ncol(data)				# The number of variables for input data
  num_rules <- length(rules)	# The number of input rules 
  
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





