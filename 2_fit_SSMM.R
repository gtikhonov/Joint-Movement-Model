# setwd("D:/HY-data/OVASKAIN/all stuff/manuscripts/Submitted/Danielle Ramos/Rcode")
#.libPaths(c("D:\\HY-data\\OVASKAIN\\all stuff\\software\\R packages", .libPaths()))
setwd("E:/Gleb/Joint-movement-model")

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
library(abind)

source("R_SSMM_functions/trunca.r")
source("R_SSMM_functions/mcmc.r")

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
      cl <- makeCluster(16)
      registerDoParallel(cl)
   }
}

out <- ssmm.mcmc(likelihood=loglik.single, data=data, n.iter=50, n.adapt.iter=10, n.thin=1, rotate=TRUE, parallel=parallel)


if(case==1){
   outfile <- "posteriors ssmm\\simulated_posterior.Rdata"
}
if(case==2){
   outfile <- "posteriors ssmm\\bird_posterior.Rdata"
}
if(case==3){
   outfile <- "posteriors ssmm\\moth_posterior.Rdata"
}

save(out, file = outfile)


if(parallel)
   stopCluster(cl)
