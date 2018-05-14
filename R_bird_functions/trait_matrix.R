### Build traits matrix with columns of 0 and 1 for categorical traits and numerical column for continuous variables ###

#input:
  #ns - number of species
  #nt - number of traits
  #transformation for continuous variable
  #vectors of traits (only for the studied species)
  
#add column names for T?

Tmatrix <- function(..., apply.log = TRUE){

  x <- list(...)
  n <- length(x) 
  ns <- length(x[[1]])
  
  T <- matrix(0,nrow=ns)

  for(i in 1:n){
    if(is.character(x[[i]]) | is.factor(x[[i]])){
      Ttmp <-  matrix(0,nrow=ns, ncol=length(unique(x[[i]])))
      colpos = x[[i]]
      for(j in 1:ns){
        Ttmp[j,colpos[j]] = 1 # T matrix does not include the intercept column (tk1=1)
      } 
      T <- cbind(T,Ttmp)
    } else {next}
    
    for(i in 1:n){
      if(is.numeric(x[[i]])){
        if(apply.log == TRUE) Ttmp = log(x[[i]])
        else Ttmp = x[[i]]
        
        T <- cbind(T,Ttmp)
        
      } else {next}
    }
  }
  
  colnames(T) <- NULL
  T <- T[,-1]
  return(T)
  
}
