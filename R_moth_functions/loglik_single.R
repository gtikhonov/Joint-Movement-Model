loglik.single = function(theta,fdata,effort,effortTringles,FEMbasis,ndt){
   sparse = TRUE
   kval = c(1,exp(theta[1]))
   mval = c(1,1)*exp(theta[2])
   aval = c(exp(theta[3]),exp(theta[4]))
   captureP = exp(theta[5])/(1+exp(theta[5]))
   
   Agrid = FEMbasis$Agrid
   habitat = FEMbasis$habitat
   MRRM = FEMbasis$MRRM
   nDAY = ncol(effort)
   
   ntriangles = nrow(FEMbasis$mesh$triangles)
   nnodes = nrow(FEMbasis$mesh$nodes)
   
   mult = fdata[,2]
   ni = nrow(fdata)
   fdata = as.matrix(fdata[,3:(2+nDAY)])
   StartingDays = which(fdata>0, arr.ind=TRUE)
   StartingDay = min(StartingDays[,2])
   
   habitat1 = habitat==1
   k = kval[1]*(habitat1) + kval[2]*(!habitat1)
   mortality = mval[1]*(habitat1) + mval[2]*(!habitat1)
   diffusion = aval[1]*(habitat1) + aval[2]*(!habitat1)
   
   K = construct_mass(FEMbasis,Agrid,sparse,k,rep(1,ntriangles))
   M = construct_mass(FEMbasis,Agrid,sparse,k,-mortality) + construct_stiff(FEMbasis,Agrid,sparse,k,-diffusion)
   
   int0 = construct_int(FEMbasis,k)
   MRRMint0 = MRRM*matrix(int0,nrow=nrow(MRRM),ncol=ncol(MRRM))
   MRRMint0ColSums = colSums(MRRMint0)
   MRRMcorArea = MRRM / matrix(MRRMint0ColSums, nrow=nrow(MRRM), ncol=ncol(MRRM), byrow=TRUE)
   
   theta = 2/3
   dt = 1/ndt
   
   D1 = t(K) - theta*dt*t(M)
   D2 = t(K) + (1-theta)*dt*t(M)
   CD1 = Cholesky(D1)
   
   like = rep(0,ni)
   v = matrix(NA,nrow=nnodes,ncol=0)
   for (day in  StartingDay:(nDAY-1)){
      ind = which(fdata[,day] > 0)
      if(length(ind)>0){
         StartingTrap = fdata[ind,day]
         vNew = as.matrix(MRRMcorArea[,StartingTrap])
         ind1 = ind<=ncol(v)
         v[,ind[ind1]] = vNew[,ind1]
         v = cbind(v,vNew[,!ind1])
      }
      for (i in 1:ndt){
         u = D2 %*% v
         v = solve(CD1, u) 
      }
      ceffort = effort[,day+1]
      searched = effortTringles[,day+1]
      nsearched = sum(searched)

      if (nsearched>0){
         trapP = as.matrix(crossprod(MRRMint0, v))
         trapP = pmax(trapP, 0)
         CtrapP = as.vector(ceffort %*% trapP)
         for (ind in 1:ni){
            da = fdata[ind,]
            if (da[day]>-1){
               cda = da[day+1]
               if (cda==0){
                  like[ind] = like[ind]+log(1-captureP*CtrapP[ind])
               }
               if (cda>0){
                  like[ind] = like[ind]+log(captureP*trapP[cda,ind])
               }
            }
         }
         smult = ((1-captureP)/(1-captureP*CtrapP))
         nsmult = (1/(1-captureP*CtrapP))
         mul = rbind(nsmult, smult)
         v = mul[searched+1,]*v
      }
   }
   li = sum(mult*like)
   return(li)
}
