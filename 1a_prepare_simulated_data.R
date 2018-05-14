set.seed(1)
setwd("D:/HY-data/OVASKAIN/all stuff/manuscripts/Submitted/Danielle Ramos/Rcode")
.libPaths(c("D:\\HY-data\\OVASKAIN\\all stuff\\software\\R packages", .libPaths()))
library(mvtnorm) # for function dmvnorm and rmvnorm

ns <- 50
np <- 3
nt <- 2
ny <- 20
TT <- matrix(rnorm(ns*nt),ns,nt)
TT[,1] <- 1
colnames(TT) <- c("intercept","trait1")
parameternames <- c("p1","p2","p3","p4")
SP <- paste("species_",toString(1),sep="")
for (sp in 2:ns){
  SP <- c(SP,paste("species_",toString(sp),sep=""))
}

CX <- matrix(rnorm(ns*ns),ns,ns)
CC <- t(CX)%*%CX
#TRUECC <- diag(rep(1,ns))
TRUEGAMMA <- matrix(rnorm(nt*np),nt,np)
TRUEMU<-TT%*%TRUEGAMMA
VX <- matrix(rnorm(np*np),np,np)
TRUESIGMA <- t(VX)%*%VX

me<-as.vector(t(TRUEMU)) #FIRST ALL NP PARAMETERS FOR SPECIES 1, THEN FOR SPECIES 2, ...
va<-kronecker(CC,TRUESIGMA) #FIRST ALL NP PARAMETERS FOR SPECIES 1, THEN FOR SPECIES 2, ...
the<-rmvnorm(1,me,va)
TRUETHETA<-t(matrix(the,np,ns))
TRUERHO<-1

X <- matrix(rnorm(ny*np),ny,np)
X[,1] <- 1
Y <- X%*%t(TRUETHETA)+matrix(rnorm(ny*ns),ny,ns)

SP <- 1:ns

data <- list(ns = ns, np = np, X = X, Y = Y, CC = CC, TT = TT, includeTraits = TRUE, includePhylogeny = TRUE,
             parameternames = parameternames, TRUEGAMMA = TRUEGAMMA, TRUESIGMA = TRUESIGMA, TRUETHETA = TRUETHETA, TRUERHO = TRUERHO, SP = SP)
save(data, file = "data\\simulated_data.Rdata")
