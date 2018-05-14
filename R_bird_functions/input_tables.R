### Delete species from the tracks table not present in the atrib table or with no traits values 
  # or not in the corr matrix ### Delete species from the atributes and phylogenetic tables not 
  # present in the tracks table

# input:
  # matrix of tracks; columns: specie/track ID/position (pixel number); the position should be 
    # in the correct order of movement and the order of the lines represent the order of movement
  # phylogenetic correlation table
  # traits table (first column are species names; other columns are numerical or categorical 
    # traits)
  # does not require one continuous variable (mass)
  
set.input <- function (tracks = tracks, traits = traits, phy.cor = cor){
  
  # check if tracks table has 3 columns
  if(!ncol(tracks) == 3) {
    stop("Error: the table with tracks should have 3 columns")
  } else {
    # check if 1 column is cathegorical (species names), and if other columns are numerical 
      # (track ID/position (pixel number))
    if(!is.factor(tracks[,1])|!is.numeric(tracks[,2])|!is.numeric(tracks[,3])){
      stop("Error: the first column of the table with tracks should contain the species 
identification, the second the tracks identification, and the third the position")
    } else {
      # check if first line of phy.cor is categorical (species names) ########### to be done
             
        TA1 = NA #tracks with more than 1 point
        SPS1 = NA #species identification for tracks
        #SP1 = NA
        TID1 = NA #track ID
        traits1 = traits[-c(1:nrow(traits)),]
        #ma1 = NA
        #di1 = NA
        cindex = NA
        SP = NA
        prob1 = prob2 = prob3 = NA #problems found
        
        sps <- as.vector(unique(tracks[,1]))
        ns <- length(sps)
        csp <- as.vector(unlist(phy.cor[1,]))
        cor <- phy.cor[-1,]
        mass <- traits[,2]
        
        for(i in 1:ns){
          po = traits[,1][which(traits[,1] %in% sps[i])]
          poc = which(csp==sps[i])
          if(length(po)==0) prob1 = c(prob1,sps[i])
          if(length(poc)==0) prob3 = c(prob3,sps[i])
          if(length(po)!=0 & length(poc)!=0) {
            lma = mass[which(traits[,1] %in% sps[i])]
            if(!is.na(lma)) {
              TA1 = c(TA1,tracks[,3][tracks[,1]==sps[i]])
              cindex = c(cindex, poc[1])
              SPS1 = c(SPS1,as.vector(tracks[,1][tracks[,1]==sps[i]]))
              TID1 = c(TID1, tracks[,2][tracks[,1]==sps[i]])
              #SP1 = c(SP1,sps[i])
              traits1 = rbind.data.frame(traits1,traits[which(traits[,1] %in% sps[i]),])
              SP = c(SP,sps[i])
            }
            
            else prob2 = c(prob2,sps[i])
          }
        }
                
        prob1 = matrix(prob1, dimnames = list(NULL,"species with no traits:"))
        prob2 = matrix(prob2, dimnames = list(NULL,"species with no trait values:"))
        prob3 = matrix(prob3, dimnames = list(NULL,"species with no phylogeny:"))
        results_probs <- list(prob1, prob2, prob3)
            
        TA = TA1[!is.na(TA1)] #só tracks com pelo menos 2 pontos
        #SP = SP1[!is.na(SP1)] #espécies observadas
        TID = TID1[!is.na(TID1)]
        SPS = SPS1[!is.na(SPS1)]
        SP = SP[!is.na(SP)]
        cindex = cindex[!is.na(cindex)]
        CC= cor[cindex,cindex] #same order for species in the correlation table
        
        return(list(results_probs,SP, tracks = data.frame(SPS,TID,TA), traits = traits1, CC = CC))
      }
    }
  }