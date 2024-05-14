#===============================================================================

# File name: sim.rep.R

# Purpose: Code to perform simulation studies

# Author: Ke Wan

# R version: R-4.1.2

# Require packages: bartCause, grf, causalLearning, pre, grpreg.

#===============================================================================
library(doParallel)

# Save address for 

# Simulation data
source("D:/simulation/data_generator.R")
source("D:/simulation/rule_generation.R")
source("D:/simulation/rule_ensemble.R")
source("D:/simulation/sim.code.R")
source("D:/simulation/variable_importance.R")

# save file
save.result <- "D:/simulation/"

# Parameters for creating the simulation dataset

n.can <- 1000                 # sample size 
p.can <- c(50, 100, 200, 400)    # number of variables
type.can <- c("rct","obv")    # type of dataset RCT or observational study

# parameters for continuous variables

para <- matrix(0,3,2)
para[1,] <- c(0,0)    # independent & standard normal
para[2,] <- c(0,1)    # independent & skewed normal(shape parameter 1)
para[3,] <- c(0.3,0)  # correlated 0.3 & multinormal

# outcome models (different combination of mu & tau)

model.can <- list()
mod.id <- 1
for(mu.id in 1:3){
  for(tau.id in 1:4){
    model.can[[mod.id]] <- c(paste0("model_",tau.id),paste0("model_",mu.id))
    mod.id <- mod.id + 1
  }
}

model.can <- model.can[1:6]
# Start simulation

for(type in type.can){
for(n in n.can){
  for(p in p.can){
    for(i in 1:nrow(para)){
      
      cor <- para[i,1]
      alpha <- para[i,2]
      
      #pr <- runif(p/2,0.25,0.75)
      
      for(model in model.can){
        
        print(paste0(type,"_",model[2],"_",model[1],"_n_",n,"_p_",p,"_cor_",cor,"_alpha_",alpha,"_start"))
        print("start_para")
        cl <- makeCluster(10)
        registerDoParallel(cl)
        res <- foreach(t = 1:100,.export=ls(envir = parent.frame())) %dopar% {sim(t)}
        stopCluster(cl)
        print("end_para")
        
        for(ii in 1:length(res)){
          
          # Output the simulation results 
          sim.res <- res[[ii]][[1]]
          summary.res <-  res[[ii]][[2]]
          vip.res <- res[[ii]][[4]]
          save.name.sim <- paste0(save.result,"res/sim_",res[[ii]][[3]])
          save.name.summary <- paste0(save.result,"res/sum_",res[[ii]][[3]])
          save.name.vip <- paste0(save.result,"res/vip_",res[[ii]][[3]])
          write.csv(sim.res,save.name.sim,row.names = FALSE)
          write.csv(summary.res,save.name.summary,row.names = FALSE)
          write.csv(vip.res,save.name.vip,row.names = FALSE)
          
          # Output the simulation dataset and corresponding parameters
          dat_train <- res[[ii]][[5]]
          dat_para_train <- res[[ii]][[6]]
          dat_test <- res[[ii]][[7]]
          dat_para_test <- res[[ii]][[8]]
          save.name.train <- paste0(save.result,"dat.train/",res[[ii]][[3]])
          save.name.test <- paste0(save.result,"dat.test/",res[[ii]][[3]])
          save.name.para_train <- paste0(save.result,"dat.train/para_",res[[ii]][[3]])
          save.name.para_test <- paste0(save.result,"dat.test/para_",res[[ii]][[3]])
          write.csv(dat_train,save.name.train,row.names = FALSE)
          write.csv(dat_para_train,save.name.para_train,row.names = FALSE)
          write.csv(dat_test,save.name.test,row.names = FALSE)
          write.csv(dat_para_test,save.name.para_test,row.names = FALSE)
          
          # Output the estimated propensity score
          p_score <- res[[ii]][[9]]
          save.name.pscore <- paste0(save.result,"p.score/",res[[ii]][[3]])
          write.csv(p_score,save.name.pscore,row.names = FALSE)
          
        }
        print(paste0(type,"_",model[2],"_",model[1],"_n_",n,"_p_",p,"_cor_",cor,"_alpha_",alpha,"_end"))
      }
    }  
  }   
}
}

