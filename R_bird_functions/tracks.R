### Build tracks table only with tracks of at least 2 positions, inside the landscape ###

#input: 
  # matrix of tracks; columns: specie/track ID/position (pixel number)
  # minimum length of the track (default:2 positions)

strack <- function (tracks = tracks, tlength = 2, lsize = 3600){
  
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
      # check if tlength is integer
      if(!is.wholenumber(tlength)){
        stop("Error: tlength should be wholenumber")
      } else{
        
        TA = NA #tracks with more than 1 point
        SPS = NA #species identification for tracks
        #SP = NA
        TID = NA #track ID
        
        sp = tracks[,1]
        sps = sp #list of observed species
        tID = tracks[,2]
        tIDs = unique(tID)
        z = tracks[,3]
        
        for (i in 1:length(tIDs)){
          take = tID[tID==tIDs[i]]
          
          if(length(take)>(tlength-1)){
            tmp = z[which(tID==tIDs[i])] #column z
            tr = numeric()
            for(j in 1:length(tmp)){
              if(tmp[j]<3601) tr = c(tr,tmp[j])
            }
            if(length(tr)>1){
              TA = c(TA,tr)
              #SP[length(SP)+1] = as.vector(unique(sp[tracks$ID==inds[i]]))
              SPS= c(SPS,as.vector(sp[tID==tIDs[i] & z<(lsize+1)]))
              TID = c(TID,tID[tID==tIDs[i] & z<(lsize+1)])
              if(length(unique(sp[tID==tIDs[i]])) > 1) print("problem with species names") 
            }
          }
        }
        
        TA = TA[!is.na(TA)] #só tracks com pelo menos 2 pontos
        #SP = SP[!is.na(SP)] #espécies por track
        TID = TID[!is.na(TID)] #identificação do track
        SPS = SPS[!is.na(SPS)] #vetor de espécies por ponto
        
        return(data.frame(SPS,TID,TA))
      }
    }
  }
}
