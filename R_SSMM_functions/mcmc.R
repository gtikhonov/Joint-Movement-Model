ssmm.mcmc <- function(likelihood = likelihood, data = data, n.iter = n.iter, n.adapt.iter = iter, n.thin = 1, rotate = TRUE, parallel=TRUE){
   #optional input arguments: muZT, VZT, PSI, nu, rhos, priorRho
   ndt = 4
   
   tpm0 <- proc.time() #initial time
   ns <- data$ns
   np <- data$np
   SP <- data$SP
   
   # Prior for M: M ~ N(muM, VM)
   muM <- rep(0,np)
   VM  <- diag(np)
   iVM = chol2inv(chol(VM))
   
   # Prior for SI: SI~IW(PSI,nu)
   PSI <- diag(np)
   nu  <- np
   
   #Initial values for THETA parameters to be estimated : THETA, SI, ZT, RH
   if (is.null(data$THETAINIT)){
      THETA = matrix(0,ns,np)
   } else {
      THETA = data$THETAINIT
   }
   SI = diag(np) #identity matrix
   iSI = solve(SI) #inverse of matrix SI
   M = rep(0,np)

   ##################
   #   MCMC chains  #
   ##################
   
   print("sampling starts")
   r <- foreach(k=icount(ns), .multicombine=TRUE,  .export=c('construct_mass','construct_stiff','construct_int',"trunca"),
      .packages=c('Matrix',"mvtnorm","MCMCpack")) %dopar% {
   # for(k in 1:ns){
      spData = data$alldata2[data$alldata2[,1]==k,]
      theta = THETA[k, ]
      li1 = likelihood(theta,spData,data$effort,data$effortTringles,data$FEMbasis,ndt)
      print("inital likelihood: ")
      print(li1)
      res <- M - theta
      li2 <- -0.5 * res %*% iSI %*% res
      
      ac = matrix(0,np,2)  # for acceptance rates
      kk = rep(1, np)       # sd for THETA proposal distributions
      acr = rep(NA, np)
      ns1s2 = 0
      s1 = rep(0,np)
      s2 = matrix(0,np,np)
      la = rep(1, np)
      vect = diag(np)
      
      POST_THETA = matrix(NA,n.iter,np)
      POST_LIKE = rep(NA,n.iter)

      for (i in 1:(n.iter+n.adapt.iter)){
         print(i)
         for (ii in 1:n.thin){
            #UPDATE THETA
            for (l in 1:np){
               mult <- rnorm(1, mean=0, sd = (kk[l]*sqrt(la[l])))
               ntheta = theta + mult*vect[l,]
               
               nli1 = likelihood(ntheta,spData,data$effort,data$effortTringles,data$FEMbasis,ndt)
               ac[l,1] = ac[l,1] + 1
               
               RES <- M - ntheta
               nli2 <- -(1/2) * RES %*% iSI %*% RES
               if(is.finite(nli1)){
                  if(runif(1) < exp((nli1 - li1 + nli2 - li2))){
                     theta = ntheta
                     li1 = nli1
                     li2 = nli2
                     ac[l,2] = ac[l,2] + 1
                  }
               }
            }    
            
            #UPDATE M
            Vs = chol2inv(chol(iVM + iSI))
            mus = Vs%*%(iVM%*%muM + iSI%*%theta)
            ##set.seed(42)
            M = as.vector(rmvnorm(1,mean=mus, sigma = Vs))
            
            #UPDATE SI
            
            RES = theta - M
            A = tcrossprod(RES)
            PSIA = PSI + A
            ##set.seed(42)
            SI = riwish((nu+1), PSIA) #InverseWishartMatrixDistribution;
            iSI = chol2inv(chol(SI))

            #Adaptation
            if (i <= n.adapt.iter){
               q <- 1 + exp(-(i*n.thin)/500)
               w <- 1 - 0.1*exp(-(i*n.thin)/500)
               acr <- ac[,2]/ac[,1]
               kk <- sapply(kk * q^(acr - 0.44), trunca)
               s1 <- s1+w*theta
               s2 <- s2+w*tcrossprod(theta)
               ns1s2 <- ns1s2+w
               if (rotate && ((i*n.thin)>50)){
                  cov <- (s2-(s1%*%t(s1)/ns1s2))/(ns1s2-1)
                  met <- cov + 10^(-5)*diag(np)
                  lavect <- eigen(met)
                  la <- abs(lavect$values)
                  vect <- lavect$vectors
               }
            }
            ac = ac*w
            
         } #END OF THINNING LOOP
         
         if(i > n.adapt.iter) {
            POST_THETA[i-n.adapt.iter,] = theta
            POST_LIKE[i-n.adapt.iter] = li1
         }
      }
      list(like=POST_LIKE, theta=POST_THETA)
   }
   
   tpm1 <- proc.time() #final time
   
   # Organizing output
   
   dtpm <- (tpm1[3] - tpm0[3])/60 #time elapsed
   
   rnames_summary <- c("Time elapsed (minutes):", "Total number of iterations:",
      "Number of adaptation iterations:","Length of thinned posterior",
      "Number of species", "Number of parameters")
   results_summary <- matrix(c(as.numeric(dtpm), (n.iter+n.adapt.iter)*n.thin, n.adapt.iter*n.thin, n.iter, ns, np),
      dimnames = list(rnames_summary,"Model summary"))
   
   bind3 = function(...){abind(..., along=3)}
   POST_THETA = do.call(bind3, lapply(r, function(c) c$theta))
   POST_LIKE = do.call(cbind, lapply(r, function(c) c$like))
   
   posterior <- list(THETA = POST_THETA, LIKE = POST_LIKE)
   
   return(list(results_summary=results_summary, posterior=posterior))
}
