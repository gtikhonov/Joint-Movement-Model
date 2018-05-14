# setwd("D:/HY-data/OVASKAIN/all stuff/manuscripts/Submitted/Danielle Ramos/Rcode")
#.libPaths(c("D:\\HY-data\\OVASKAIN\\all stuff\\software\\R packages", .libPaths()))
setwd("C:/Google Drive/MBG/Hierarchical movement model/Rcode")

library(rgdal)
library(MCMCpack) # for function riwish
library(mvtnorm) # for function dmvnorm and rmvnorm
library(stats)
library(grDevices)
library(graphics)
library(rgl)
library(MASS)
library(fdaPDE)
library(Matrix)
library(doParallel)

source("R_JSMM_functions/trunca.r")
source("R_JSMM_functions/mcmc.r")

case <- 3 #1: simulated; 2: birds; 3: moths

#SIMULATED CASE STUDY
if (case==1){
   load(file = "data\\simulated_data.Rdata")
   source("R_simulated_functions/loglik.r")
}

#BIRD CASE STUDY
if(case==2){
   load(file = "data\\bird_data.Rdata")
   source("R_bird_functions/loglik.r")
}

#MOTH CASE STUDY
if(case==3){
   load(file = "data\\moth_data.Rdata")
   source("R_moth_functions/construct_int.R")
   source("R_moth_functions/construct_mass.R")
   source("R_moth_functions/construct_stiff.R")
   source("R_moth_functions/make_image.R")
   source("R_moth_functions/loglik_single.R")
   parallel <- TRUE
   if(parallel){
      cl <- makeCluster(4)
      registerDoParallel(cl)
      source("R_moth_functions/parallel_loglik.R")
   } else {
      source("R_moth_functions/loglik.R")
   }
}

if(parallel){
   out <- jsmm.mcmc(likelihood = parallel_loglik, data = data, n.iter = 500, n.adapt.iter = 100, n.thin = 10, rotate = TRUE)
} else {
   out <- jsmm.mcmc(likelihood = loglik, data = data, n.iter = 100, n.adapt.iter = 100, n.thin = 1, rotate = TRUE)
}

if(case==1){
   outfile <- "posteriors\\simulated_posterior.Rdata"
}
if(case==2){
   outfile <- "posteriors\\bird_posterior.Rdata"
}
if(case==3){
   outfile <- "posteriors\\moth_posterior.Rdata"
}

save(out, file = outfile)


if(parallel)
   stopCluster(cl)
