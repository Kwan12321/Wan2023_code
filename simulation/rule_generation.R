#===============================================================================

# File name: rules_generation.R

# Purpose: Rules generation algorithm for proposed method

# Author: Ke Wan

# R version: R-4.1.2

# Input: data: Data frame; 
#        meandepth: Average depth of each tree-based learner;
#        shrinkage: Shrinkage rate for each boosting steps; ntrees: Maximum number of trees to generate;
#        sampfrac: Size of training sample fraction for each boosting step;
#        pscore: Estimated propensity score for estimating the transformed outcome.

# Output: a list including ---> rules: generated rules; rulevars: rules matrix (transform rules into binary values(0,1)) 

# Require packages: rpart

#===============================================================================

rules.generate <- function(data,				# Input data frame
                           meandepth = 3, 		# The mean depth of each base learner(decision tree)
                           shrinkage = 0.01,   # shrinkage rate for each boosting steps    
                           ntrees = 500, 		# Maximum number of trees to generate                          
                           sampfrac = NULL,    # Size of training sample in each boosting steps
                           pscore = 0.5		# Propensity score
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
  
  print(3)
  colnames(rulevars) <- rules			 # Rename the rules matrix
  
  rules.set <- list(rules = rules,rulevars = rulevars)
  
  return(rules.set)
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








