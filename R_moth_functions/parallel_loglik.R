parallel_loglik <- function(THETA,data){
  ns <- data$ns
  ndt<-4
  ncu<-max(data$cusp[,1])
  r <- foreach(cu=icount(ncu), .combine=cbind, .export=c('loglik.single','construct_mass','construct_stiff','construct_int'), .packages='Matrix') %dopar% {
    loglik.single(THETA[data$cusp[cu,2],],data$alldata2[data$alldata2[,1]==cu,],data$effort,data$effortTringles,data$FEMbasis,ndt)
  }
  
  like<-rep(0,ns)
  for (sp in 1:ns){
    like[sp]<-sum(r[which(data$cusp[,2]==sp)])
  }
  return(like)
}
