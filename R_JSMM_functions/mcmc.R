jsmm.mcmc <- function(likelihood = likelihood, data = data, n.iter = n.iter, n.adapt.iter = iter, n.thin = 1, rotate = TRUE){
  #optional input arguments: muZT, VZT, PSI, nu, rhos, priorRho
  
  tpm0 <- proc.time() #initial time
  ns <- data$ns
  np <- data$np
  includePhylogeny <- data$includePhylogeny
  includeTraits <- data$includeTraits
  CC <- data$CC
  SP <- data$SP
  
  if (!includeTraits){
    TT<-matrix(1,data$ns,1)
    colnames(TT)<- "intercept"
    nt <- 1
  } else {
    TT <- data$TT
    nt <- dim(TT)[2]
  }
  
  # Prior for ZT: ZT ~ N(muZT, VZT)
  muZT <- matrix(0, nrow = np * nt)
  VZT  <- diag(np * nt)
  
  # Prior for SI: SI~IW(PSI,nu)
  PSI <- diag(np)
  nu  <- np
  
  # Pre-compute inverse and determinant of the D-matrix
  
  if(includePhylogeny){
    rhos = seq(0, 1, 0.01)
    nr = length(rhos)
    priorRho = rep(0.5/(nr - 1),nr)
    priorRho[1] = 0.5
    lpriorRho = log(priorRho)
    iDDs = array(dim = c(ns,ns,nr))
    ldetDDs = numeric(nr)
    
    for(i in 1:nr){
      DD = rhos[i]*CC + (1 - rhos[i])*diag(ns)
      iDDs[,,i] = solve(DD)
      ldetDDs[i] = log(det(as.matrix(DD)))
    }
  }
  
  XT = kronecker(TT,diag(np))
  
  #Initial values for THETA parameters to be estimated : THETA, SI, ZT, RH
  if (is.null(data$THETAINIT)){
    THETA = matrix(0,nrow = ns,ncol=np)
  } else {
    THETA = data$THETAINIT
  }
  SI = diag(np) #identity matrix
  iSI = solve(SI) #inverse of matrix SI
  ZT = matrix(0, nrow=nt, ncol=np)
   
  if(includePhylogeny){
    RI = round(nr/2)
    RH = rhos[RI]
    iDD = iDDs[,,RI]
    ldetDD = ldetDDs[RI]
  } else {
    iDD <- diag(ns)
    RH <- 1
  }

  ac = array(0, c(ns, np, 2))  # for acceptance rates
  kk = array(1, c(ns, np))       # sd for THETA proposal distributions
  acr = array(NA, c(ns, np))
  ns1s2 = array(0,c(ns))
  s1 = array(0,c(ns,np))
  s2 = array(0,c(ns,np,np))
  la = array(1, c(ns,np))
  vect = array(0,c(ns,np,np))
  for (k in 1:ns){
    vect[k,,]<-diag(np)
  }
  
  POST_THETA = array(NA,c(n.iter,ns,np))
  POST_ZT = array(NA,c(n.iter,nt,np))
  POST_SI = array(NA,c(n.iter,np,np))
  POST_RHO = array(NA,c(n.iter))
  POST_LIKE = array(NA,c(n.iter,ns))
  
  ##################
  #   MCMC chains  #
  ##################
  
  li1 = likelihood(THETA, data)
  M = TT%*%ZT
  RES <- as.numeric(t(M)) - as.numeric(t(THETA))
  li2 <- -(1/2) * RES %*% kronecker(iDD, iSI) %*% RES
  
  print("inital likelihood: ")
  print(c(li1,li2))
  print("sampling starts")
  
  for (i in 1:(n.iter+n.adapt.iter)){
    print(i)
    for (ii in 1:n.thin){
      #UPDATE THETA
      for (l in 1:np){
        NTHETA = THETA
        for(k in 1:ns){
          nTHETA = THETA[k, ]
          mult <- rnorm(1, mean=0, sd = (kk[k,l]*sqrt(la[k,l])))
          for (l2 in 1:np){
            nTHETA[l2] = nTHETA[l2] + mult*vect[k,l,l2]
          }
          NTHETA[k, ] = nTHETA
        }
        
        nli1 = likelihood(NTHETA, data)
          
        for (k in 1:ns){
          N2THETA <- THETA
          N2THETA[k,] <- NTHETA[k,]
          RES <- as.numeric(t(M)) - as.numeric(t(N2THETA))
          nli2 <- -(1/2) * RES %*% kronecker(iDD, iSI) %*% RES
          ac[k,l,1] = ac[k,l,1] + 1

          if(is.finite(nli1[k]) & is.finite(nli2)){
            if(runif(1) < exp((nli1[k] - li1[k] + nli2 - li2))){
              THETA[k,] = NTHETA[k,]
              li1[k] = nli1[k]
              li2 = nli2
              ac[k,l,2] = ac[k,l,2] + 1
            }
          }
        }
      }    
      
      #UPDATE ZETA
      
      XTHETA = as.vector(t(THETA)) #transform to vector by line
      iXSI = kronecker(iDD,iSI) #iDD - Rho
      Vs = solve(solve(VZT) + t(XT)%*%iXSI%*%XT)
      Vs = (Vs + t(Vs))/2
      mus = Vs%*%(solve(VZT)%*%muZT + t(XT)%*%iXSI%*%XTHETA)
      ##set.seed(42)
      ZT = matrix(rmvnorm(1,mean=mus, sigma = Vs), ncol = np, byrow=TRUE)
      M = TT%*%ZT
      
      #UPDATE SI
      
      RES = THETA - M
      A = t(RES)%*%iDD%*%RES
      PSIA = PSI + A
      PSIA = (PSIA + t(PSIA))/2
      ##set.seed(42)
      SI = riwish((nu+ns), PSIA) #InverseWishartMatrixDistribution;
      SI = (SI + t(SI))/2
      iSI = solve(SI)
      
      #UPDATE RHO
      
      if(includePhylogeny){
        RES = as.numeric(t(M)) - as.numeric(t(THETA))
        likeRho = numeric(nr)
        for(ii in 1:nr){
          likeRho[ii] = (-1/2)*(np*ldetDDs[ii] + RES%*%kronecker(iDDs[,,ii], iSI)%*%RES)
        }
        postRho = lpriorRho + likeRho
        pr = exp(postRho)/sum(exp(postRho))
        RI = sample(seq(1:nr),size=1,prob=pr)
        iDD = iDDs[,,RI]
        ldetDD = ldetDDs[RI]
        RH = rhos[RI]
        RES <- as.numeric(t(M)) - as.numeric(t(THETA))
        li2 <- -(1/2) * RES %*% kronecker(iDD, iSI) %*% RES
      }
      
      #Adaptation
      if (i <= n.adapt.iter){
        for(k in 1:ns){
          q <- 1 + exp(-(i*n.thin)/500)
          w <- 1 - 0.1*exp(-(i*n.thin)/500)
          acr[k,] <- ac[k,,2]/ac[k,,1]
          kk[k,] <- sapply(kk[k,] * q^(acr[k,] - 0.44), trunca)
          s1[k,] <- s1[k,]+w*THETA[k,]
          s2[k,,] <- s2[k,,]+w*(THETA[k,]%*%t(THETA[k,]))
          ns1s2[k] <- ns1s2[k]+w
          if (rotate && ((i*n.thin)>50)){
            cov <- (s2[k,,]-(s1[k,]%*%t(s1[k,])/ns1s2[k]))/(ns1s2[k]-1)
            met <- cov + 10^(-5)*diag(np)
            lavect <- eigen(met)
            la[k,] <- abs(lavect$values)
            vect[k,,] <- lavect$vectors
          }
        }
        ac = ac*w
      }
      
    } #END OF THINNING LOOP
    
    if(i > n.adapt.iter) {
      POST_THETA[i-n.adapt.iter,,] = THETA
      POST_ZT[i-n.adapt.iter,,] = ZT
      POST_SI[i-n.adapt.iter,,] = SI
      POST_RHO[i-n.adapt.iter] = RH
      POST_LIKE[i-n.adapt.iter,] = li1
    }
    
  }
  
  tpm1 <- proc.time() #final time
  
  # Organizing output
  
  dtpm <- (tpm1[3] - tpm0[3])/60 #time elapsed
  
  rnames_summary <- c("Time elapsed (minutes):", "Total number of iterations:",
                      "Number of adaptation iterations:","Length of thinned posterior",
                      "Number of species", "Number of parameters", "Number of traits")
  results_summary <- matrix(c(as.numeric(dtpm), (n.iter+n.adapt.iter)*n.thin, n.adapt.iter*n.thin, n.iter, ns, np, nt),
                            dimnames = list(rnames_summary,"Model summary"))
  
  posterior <- list(THETA = POST_THETA, GAMMA = POST_ZT, SIGMA = POST_SI, RHO = POST_RHO, LIKE = POST_LIKE)
  
  return(list(results_summary = results_summary, posterior = posterior, la = la, vect = vect, ac = ac, kk = kk))
}
