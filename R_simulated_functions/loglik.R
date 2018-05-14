loglik <- function(THETA,data){
  
  X <- data$X
  Y <- data$Y
  ns <- dim(Y)[2]
 
  meY <- X%*%t(THETA)

  like<-rep(0,ns)
  for (k in 1:ns) {
    res <- meY[,k]-Y[,k]
    like[k]<-sum(dnorm(res,mean=0,sd=1,log=TRUE))
  }
  
#  res <- data$TRUETHETA - THETA
#  for (k in 1:ns) {
#    sisi <- diag(np)
#    sisi[1,2] <- 0.99
#    sisi[2,1] <- 0.99
#    like[k] <- dmvnorm(res[k,],mean=c(0,0),sigma=sisi,log=TRUE)
#  }
  return(like)
}
