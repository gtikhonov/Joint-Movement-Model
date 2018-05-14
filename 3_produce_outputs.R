setwd("D:/HY-data/OVASKAIN/all stuff/manuscripts/Submitted/Danielle Ramos/Rcode")
#.libPaths(c("D:\\HY-data\\OVASKAIN\\all stuff\\software\\R packages", .libPaths()))
#library(rgdal)

case <- 1 #1: simulated; 2: birds; 3: moths
if (case==1){
  datafile <-  "data\\simulated_data.Rdata"
  outfile <- "posteriors\\simulated_posterior.Rdata"
  truevalues <- TRUE
}
if (case==2){
  datafile <-  "data\\bird_data.Rdata"
  outfile <- "posteriors\\bird_posterior.Rdata"
  truevalues <- FALSE
}

load(file = datafile)
load(file = outfile)
post<-out$posterior

nt <- dim(post$GAMMA[1,,])[1]
np <- data$np
ns <- data$ns

sp = 1:min(3,ns)
par(mfrow=c(length(sp),np))
for (k in sp){
  for (l in 1:np){
    mi=min(post$THETA[,k,l])
    ma=max(post$THETA[,k,l])
    if (truevalues){
      mi=min(mi,data$TRUETHETA[k,l])
      ma=max(ma,data$TRUETHETA[k,l])
    }
    plot(post$THETA[,k,l],type='l',ylim=c(mi,ma))
    if (truevalues) {abline(a=data$TRUETHETA[k,l],b=0,col='red')}
  }
}

tr=1:min(3,nt)
par(mfrow=c(length(tr),np))
for (k in tr){
  for (l in 1:np){
    mi=min(post$GAMMA[,k,l])
    ma=max(post$GAMMA[,k,l])
    if (truevalues){
      mi=min(mi,data$TRUEGAMMA[k,l])
      ma=max(ma,data$TRUEGAMMA[k,l])
    }
    plot(post$GAMMA[,k,l],type='l',ylim=c(mi,ma))
    if (truevalues) {abline(a=data$TRUEGAMMA[k,l],b=0,col='red')}
  }
}

par(mfrow=c(1,1))
mi=0
ma=1
plot(post$RHO,type='l',ylim=c(mi,ma))
if (truevalues) {abline(a=data$TRUERHO,b=0,col='red')}

par(mfrow=c(np,np))
for (k in 1:np){
  for (l in 1:np){
    mi=min(post$SIGMA[,k,l])
    ma=max(post$SIGMA[,k,l])
    if (truevalues) {
      mi=min(mi,data$TRUESIGMA[k,l])
      ma=max(ma,data$TRUESIGMA[k,l])
    }
    plot(post$SIGMA[,k,l],type='l',ylim=c(mi,ma))
    if (truevalues) {abline(a=data$TRUESIGMA[k,l],b=0,col='red')}
  }
}

#par(mfrow=c(1,1))
#plot(post$THETA[,1,],type='l')
print(out$results_summary)
