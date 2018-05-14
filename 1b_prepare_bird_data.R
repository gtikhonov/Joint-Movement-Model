setwd("D:/HY-data/OVASKAIN/all stuff/manuscripts/Submitted/Danielle Ramos/Rcode")

source("R_bird_functions/iswholenumber.r")
source("R_bird_functions/tracks.r")
source("R_bird_functions/input_tables.r")
source("R_bird_functions/trait_matrix.r")

xp <- rep(c(1:60), 60) # x coordinates
yp <- rep(c(1:60), each=60) # y coordinates
xv = cbind(xp,yp)
di <- as.matrix(dist(xv,upper=T,diag=T))

landscape = "pa"

habitat = read.table(paste("bird_data/type_", landscape, ".csv", sep=""), sep = ",") #landscape data
h = habitat
h[h==3|h==4] = 2 #corridors,isolated trees and groups of trees are corridors (or stepping-stones)
h[h==5|h==6] = 3 #forest patches are forest
hh=h

corridor = h #matrix; each pixel has a value =1 for corridor pixels or =0 for patch or pasture pixels
corridor[corridor==1|corridor==3]=0
corridor[corridor==2]=1
corridor <- unlist(corridor)

patch = h #matrix; each pixel has a value =1 for patch pixels or =0 for corridor or pasture pixels
patch[patch==1|patch==2]=0
patch[patch==3]=1
patch <- unlist(patch)

tracks = read.table(paste("bird_data/tracks_", landscape, ".csv", sep=""), sep = ",", header = TRUE) #movement data
tracks = tracks[,c(2,3,5)]
tracks = tracks[!tracks$sp=="Guira guira",]
tracks <- strack(tracks = tracks, tlength = 2, lsize = 3600)

cor = read.table("bird_data/correlation_mean.csv", sep = ",", header=FALSE)

# species traits
traits = read.table("bird_data/traits.csv", sep = ",", header=TRUE) #species traits data
traits <- traits[,c(5,10,11)]

input <- set.input(tracks = tracks, traits = traits, phy.cor = cor)
problems <- input[[1]]
tracks <- input$tracks
traits <- input$traits
CC <- apply(input$CC,FUN = as.numeric, MARGIN = 1) ################## correct this in the code set.input
TT <- Tmatrix(diet = traits[,2], mass = traits[,3], apply.log = TRUE) #### changed T to TT

SPS <- tracks[ , 1]
SP <- unique(SPS)
ns <- length(SP)
parameternames <- c("log distance","corridor_aff", "forest_aff")
np <- length(parameternames)

data <- list(SP = SP, ns = ns,  np = np, corridor = corridor, patch = patch, di = di,  CC = CC, TT = TT, tracks = tracks, includeTraits = TRUE, includePhylogeny = TRUE,
             parameternames = parameternames)
save(data, file = "data\\bird_data.Rdata")
