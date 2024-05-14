#===============================================================================

# File name: proposed_method.R

# Purpose: Code for proposed method

# Author: Ke Wan

# R version: R-4.1.2

# Input: data: Input Data frame; 
#        meandepth: Average depth of each tree-based learner;
#        shrinkage: Shrinkage rate for each boosting steps; ntrees: Maximum number of trees to generate;
#        sampfrac: Size of training sample fraction for each boosting step;
#        pscore: Estimated propensity score for estimating the transformed outcome.
#        wins: Winsorized size (0.025 is recommended in Freidman and Popescu, 2008)
#        nfold: The number of cross-validation folds in adaptive weight calculation

# Output: a list including --->  beta: estimated intercept and coefficient, 
#                                lambda.best: selected tuning parameter lambda, 
#                                rules: generated rules,
#                                scale_x: standard deviation for each variables in train data 

# Require packages: rpart, grpreg.

#===============================================================================

prop  <-      function(data,				      # Input data frame
                       meandepth = 3, 		# The mean depth of each base learner(decision tree)
                       shrinkage = 0.01,  # shrinkage rate for each boosting steps    
                       ntrees = 333, 		  # Maximum number of trees to generate                          
                       sampfrac = NULL,   # Size of training sample in each boosting steps
                       pscore = 0.5,		  # Estimated Propensity score
                       wins = 0.025,      # Winsorized size
                       nfold = 5          # The number of cross-validation folds in adaptive weight calculation
){
  
  # Summary of data-------------------------------------------------------------
  
  n <- nrow(data)							# Sample size
  p <- ncol(data[,-c(1:2)])				# Number of variables
  
  # Extract the variable names--------------------------------------------------
  
  y_names <- colnames(data)[1]			# Outcome name
  g_names <- colnames(data)[2]			# Treatment indicator name
  x_names <- colnames(data)[-c(1:2)]		# Covariates names
  
  # Calculate the transformed outcome and create the transformed outcome data---
  
  y <- data[,y_names]						# Outcome
  g <- data[,g_names]						# Treatment indicator
  x <- data[,x_names]						# Covariates
  
  trans_y <- g*(y/pscore) - (1-g)*(y/(1-pscore)) # Calculate transformed outcome
  
  data.learn <- data.frame(y = trans_y,x)        # Create the transformed outcome data
  colnames(data.learn) <- c(y_names,x_names)     # Rename the variables of transformed outcome data
  
  
  ###########################################################################################
  ## Step 1: Rule generation (refer the code in "Pre" packages)
  ###########################################################################################
  
  # Stochastic Gradient Boosting Trees with random depth base learners (Freidman and Popescu,2008)------
  
  # Determine the training sample size for each boosting steps
  
  if(is.null(sampfrac)){
    size <- min(n/2,100 + 6*sqrt(n))	# If the sampfrac is not specified,it will be simply determined as that in Freidman and Popescu (2008)
  }else{
    size <- sampfrac*n
  }			   					
  
  # Create the sample id for training samples for each boosting steps
  
  subsample <- list()
  subsample <- mapply(function(i)
    subsample[[i]] <- sample(1:n, size = round(size), replace = FALSE)
    ,i=1:ntrees,SIMPLIFY=FALSE)
  
  # Initialize
  
  rules <- c() 								   		   # Initialize the rule
  eta0 <- mean(data.learn[,y_names])					   # Initialize the pseudo residuals
  eta <- rep(eta0, length(data[,y_names]))               # Initialize the memory function eta		
  
  # Determing the depth of each base learners (Freidman and Popescu, 2008)
  
  maxdepth <- ceiling(log(2 + floor(rexp(ntrees, rate = 1/(2^meandepth - 2))), base = 2))
  
  # Calculate the first pseudo residuals
  
  data.learn[,y_names] <- trans_y - eta
  
  # Repeating the boosting procedure "ntrees" times
  
  for(i in 1:ntrees){
    
    # Create the base learners
    
    tree <- rpart::rpart(y~., control = rpart::rpart.control(maxdepth = maxdepth[i]),data = data.learn[subsample[[i]], ])
    
    # Extract rules from the base learner
    paths <- rpart::path.rpart(tree, nodes = rownames(tree$frame), print.it = FALSE, pretty = 0)  # Extract rules from tree
    paths <- unname(sapply(sapply(paths, `[`, index = -1), paste, collapse = " & ")[-1])   # Rename each rule
    
    # Remove the complement rule (root)
    paths <- paths[-1]
    
    # Update the rules set
    rules <- c(rules,paths)
    
    # Update the memory function eta
    
    eta <- eta + shrinkage * predict(tree, newdata = data.learn)
    
    #update the pseudo residual 
    
    data.learn[,y_names] <- trans_y - eta
  }
  
  # Remove the dumplicate and complement rules (Sub_function 1)
  
  fit.rule <- remove.duplicate.complement(data,rules)
  
  # Get the rules set
  
  rules <- fit.rule$rules				 # The set of generated rules 
  rulevars <- rule.predict(data,rules) # Transform rules into binary values to get the rules matrix (Sub_function 2)
  
  colnames(rulevars) <- rules			 # Rename the rules matrix
  
  rules.set <- list(rules = rules,rulevars = rulevars)
              
  ###################################################################################################################            
  ## Reconstruct the data frame for rule ensemble    
  ###################################################################################################################                     
  
  # Get the winsorized covariates(Freideman and Popescu, 2008)------------------
  
  data.win <- winsorize(data,wins) # The winsorized function (Subfunction 3)
  data.win <- data.win[[1]]        # The data with winsorized covariates
  
  # Normalized the covariates---------------------------------------------------
  
  data.norm <- normalize(data.win) # The normalized function (Subfunction 4)
  scale_x <- data.norm$scale		 # Standard deviation for each variables 
  data.norm <- data.norm$data		 # The data with winsorized and normalized covariate
  
  # Combine the rules and linear terms------------------------------------------
  
  comb.mat <- cbind(data.norm,rulevars)
  
  # Transform the data for group lasso
  
  grp.mat <- mat_grp(comb.mat)
  
  # Summary the data
  
  grp.mat$scale <- scale_x   
  
  ####################################################################################################################
  ## Step 2: Rule ensemble
  ####################################################################################################################
  
  cv_grpreg <- grpreg::cv.grpreg(X = as.matrix(grp.mat$x_org),y = grp.mat$y,group = grp.mat$index,nfold = nfold,returnX = TRUE)
  
  # ** Calculate adaptive weights (Huang et al. 2009) ##################
  mod <- cv_grpreg$fit
  lambda <- cv_grpreg$lambda.min
  beta0 <- mod$beta[,which(mod$lambda == lambda)][-1]
  
  beta0 <- beta0*mod$XG$scale
  beta0[which(is.nan(beta0))] <- 0
  
  grp_weights <- numeric(length(unique(grp.mat$index)))
  for (i in unique(grp.mat$index)) {
    grp_weights[i] <-  1 / sqrt(sum(beta0[which(grp.mat$index == i)]^2))
  }
  ######################################################################
  
  ad_grpreg.bic <- grpreg::grpreg(X = as.matrix(grp.mat$x_org),y = grp.mat$y,group = grp.mat$index,group.multiplier = grp_weights)
  lambda <- grpreg::select(ad_grpreg.bic,criterion="BIC")$lambda
  beta <- ad_grpreg.bic$beta[,which(ad_grpreg.bic$lambda == lambda)]
  treatment.est <- as.matrix(grp.mat$x_trt)%*%beta[-1] + beta[1]
  control.est <- as.matrix(grp.mat$x_con)%*%beta[-1] + beta[1]
  tau <- treatment.est - control.est
  
  ### Summary of the estimated parameters ###
  
  return(list(beta = beta,                      # Estimated intercept and coefficient 
              lambda.best = lambda,             # Seleceted tuning parameter lambda 
              rules = rules,                    # Generated rules
              scale_x = scale_x                 # Standard deviation for each variables 
              )
        )
}

##########################################################################################
#Sub_function 1: Remove the dumplicate and complement rules (Refer to the R package "pre")
##########################################################################################				

remove.duplicate.complement <- function(data,rules){
  
  # Remove the duplicate rules--------------------------------------------------
  
  # Creat the rules matrix (F2)
  
  rulevars <- rule.predict(data,rules)
  
  # Remove the rules with same values in rules matrix
  
  duplicates <- duplicated(rulevars,MARGIN = 2)   # Find out the duplicate columns
  duplicates.removed <- rules[duplicates]
  rulevars <- rulevars[,!duplicates,drop = FALSE]
  rules <- rules[!duplicates]
  
  #remove the complement rules--------------------------------------------------
  
  # Find columns with equal variance to reduce the number of comparisons
  
  vars <- apply(rulevars, 2, var_bin) # Calculate the variance for each rules
  
  # Compare the variance equality between each rules
  
  vars_distinct <- lapply(
    unique(vars), function(x) { 
      idx <- which(is_almost_eq(x, vars)) 		# idx is the id of rule that have variance x 
      list(var = x, n = length(idx), idx = idx)	# Summary of the varaince equality (var: variance; n : number of equal variance; idx: Corresponding id) 
    }) 
  
  #store the complements
  complements <- logical(ncol(rulevars))
  
  # Identify the complements
  
  for (va in vars_distinct) {        # Loop for all rules
    
    if(va$n < 2L) next             # If there is no same exist then next
    
    idx <- va$idx
    idx <- setdiff(idx, which(complements))  # Remove the complement and update their id
    
    if(length(idx) < 2)next        # If there is no same left then next
    
    n_idx <- length(idx)		   # Number of the same
    
    # Remove the complement
    
    for(j in 1:(n_idx - 1)){       
      if (complements[idx[j]])next
      this_val <- rulevars[, idx[j]]
      is_compl <- which(apply(rulevars[, idx[(j + 1):n_idx], drop = FALSE], 2, function(x) all(x != this_val))) + j
      if (length(is_compl) > 0)
        complements[idx[is_compl]] <- TRUE
    }
  }
  
  complements <- which(complements)
  if(length(complements)!= 0){
    complements.removed <- rules[complements]
    rules <- rules[-complements]
    rulevars <- rulevars[,-complements,drop = FALSE]
  }
  
  return(list(rules = rules, rulevars = rulevars))
}

# Standard error of rule terms

var_bin <- function(x) {
  p <- mean(x)
  p*(1L-p)
}

# Check near equality

is_almost_eq <- function(x, y, tolerance = sqrt(.Machine$double.eps)) {
  stopifnot(is.numeric(x), length(x) == 1L)
  x_abs <- abs(x)
  xy <- if (x_abs > tolerance) {abs(x - y) / x_abs} else {abs(x - y)}
  xy <= tolerance
}


##################################################################################################
#Sub_function 2: Creat the rules matrix (0: against the rules; 1: obey the rules)
##################################################################################################	
rule.predict <- function(data,rules){
  
  expr <- parse(text = paste0("cbind(", paste0(rules, collapse = ", "), ")"))
  x <- eval(expr,data)
  colnames(x) <- names(rules)
  x <- apply(x,2,as.numeric)
  return(x)
}

##################################################################################################
#Sub_function 3: Winsorize function
##################################################################################################	
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

##################################################################################################
#Sub_function 4: Normalize function
##################################################################################################
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

##################################################################################################
#Sub_function 5: Transform the matrix for adaptive group lasso 
##################################################################################################
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

  
  
  
  
  
  
  
  
  
