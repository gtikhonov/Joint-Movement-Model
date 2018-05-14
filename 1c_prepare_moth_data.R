set.seed(1)
# .libPaths(c("D:\\HY-data\\OVASKAIN\\all stuff\\software\\R packages", .libPaths()))
# setwd("D:/HY-data/OVASKAIN/all stuff/manuscripts/Submitted/Danielle Ramos/Rcode")
setwd("C:/Google Drive/MBG/Hierarchical movement model/Rcode")


library(stats)
library(grDevices)
library(graphics)
library(rgl)
library(MASS)
library(fdaPDE)
library(Matrix)

source("R_moth_functions/make_FEM_basis.R")

alldata<-read.table("moth_data/captures.txt",header=FALSE)
effort <- as.matrix(read.table("moth_data/effort.txt",header=FALSE))
FEMbasis<-make_FEM_basis("moth_data/my.tp")

effortTringles = (FEMbasis$MRRM %*% effort)>0
dataList = split.data.frame(alldata, alldata$V1)
dataList = lapply(dataList, function(a){first=apply(a[,-c(1:2)]>0,1,Position,f=function(b) b); a[order(first),,drop=FALSE]})
alldata = do.call("rbind", dataList)

# these lines should be uncommneted once all data is available
# TT<-read.table("moth_data/traits.txt",header=TRUE)
# CC<-read.table("moth_data/C.txt",header=FALSE)

#these lines should be removed once all data is available
load(file = "data\\moth_data.Rdata")
TT = data$TT
CC = data$CC
rm(list="data")

ns<-max(alldata[,1])
res<-matrix(0,ns,3)

SP <- 1:ns
parameternames <- c("k_habitat2","log(mortality)","log(D1)","log(D2)","logit(captureP)")
np <- length(parameternames)

# INITIAL CONDITION THAT RETURNS A NUMERICAL LIKELIHOOD FOR ALL SPECIES
THETAINIT <- matrix(0,ns,5)
for (sp in 1:ns){
  THETAINIT[sp,]<-c(log(0.5),log(0.00001),log(100000),log(100000),log(0.8/(1-0.8)))
}

#PREPARE PARALLEL DATA STARTS
ninds<-rep(0,ns)
for (sp in 1:ns){
  ninds[sp]<-sum(alldata[,1]==sp)
}
spl<-40 #MAX NUMBER OF INDIVIDUALS TO INCLUDE IN EACH COMPUTATIONAL UNIT
nspl<-ceiling(ninds/spl)
ncu<-sum(nspl) #number of computational units for parallel computation
cusp<-matrix(0,ncu,2) #computational units belonging to each species
cu<-0
for (sp in 1:ns){
  for (lcu in 1:nspl[sp]){
    cu<-cu+1
    cusp[cu,]<-c(cu,sp)
  }
}
alldatasp<-alldata[,1]
alldatacu<-alldatasp
for (sp in 1:ns){
  cus<-which(cusp[,2]==sp)
  rows<-which(alldata[,1]==sp)
  blocksize<-round(length(rows)/nspl[sp])
  cus2<-sort(rep(cus,blocksize))
  cus2<-cus2[1:min(length(cus2),ninds[sp])]
  cus2<-c(cus2,rep(cus[nspl[sp]],ninds[sp]-length(cus2)))
  alldatacu[rows]<-cus2
}
alldata2<-alldata
alldata2[,1]<-alldatacu

data <- list(SP = SP, ns = ns,  np = np, alldata = alldata, alldata2 = alldata2, cusp = cusp, effort = effort, effortTringles=effortTringles,
   FEMbasis = FEMbasis, CC = CC, TT = TT, includeTraits = TRUE, includePhylogeny = TRUE,
   parameternames = parameternames, THETAINIT = THETAINIT)
save(data, file = "data\\moth_data.Rdata")
