set.seed(1)
#.libPaths(c("D:/HY-data/OVASKAIN/all stuff/software/R packages", .libPaths()))
# setwd("D:/HY-data/OVASKAIN/all stuff/manuscripts/Submitted/DCanielle Ramos/Rcode")
setwd("C:/Google Drive/MBG/Hierarchical movement model/Rcode")


library(stats)
library(grDevices)
library(graphics)
library(rgl)
library(MASS)
library(fdaPDE)
library(Matrix)
library(doParallel)

source("R_moth_functions/make_FEM_basis.R")
source("R_moth_functions/construct_int.R")
source("R_moth_functions/construct_mass.R")
source("R_moth_functions/construct_stiff.R")
source("R_moth_functions/make_image.R")
source("R_moth_functions/loglik_single.R")
source("R_moth_functions/loglik.R")
source("R_moth_functions/parallel_loglik.R")

load(file = "data/moth_data.Rdata")

start.time <- Sys.time()
like <- loglik(data$THETAINIT,data)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
print(sum(like))
#Time difference of 2.079317 mins
#-6033.461

cl <- makeCluster(4)
registerDoParallel(cl)

start.time <- Sys.time()
like <- parallel_loglik(data$THETAINIT,data)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
print(sum(like))

stopCluster(cl)
