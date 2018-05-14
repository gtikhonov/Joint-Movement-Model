loglik <- function(THETA,data){
  ns <- data$ns
  like<-rep(0,ns)
  ndt<-4 #Number of time steps for solving the PDE forward for one day
  for (sp in 1:ns){
    like[sp] <-loglik.single(THETA[sp,],data$alldata[data$alldata[,1]==sp,],data$effort,data$effortTringles,data$FEMbasis,ndt)
    print(like[sp])
  }
  return(like)
}
