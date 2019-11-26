# Get MCMC data:
begyear <- 1990
endyear <- 2011
predyear <- 2012

setwd(paste("/Users/blakemcshane/Documents/Academic Research/Baseball/trajectory/SplineModels_",
	begyear,endyear,"_",predyear,sep=""))
load("mcmcdat.Rdata")
mydat <- BaseballData
myyear <- max(mydat[,"Year"])
sel <- mydat[,"Year"]==myyear
myid <- as.character(sort(unique(mydat[sel,"Name"])))

# Create predicted data:
team <- rep(NA, length(myid))
pos <- rep(NA, length(myid))
age <- rep(NA, length(myid))
ab <- rep(NA, length(myid))
hand <- rep(NA, length(myid))
park <- rep(NA, length(myid))
for(i in 1:length(myid)){

  sel <- mydat[,"Name"]==myid[i] & mydat[,"Year"]==myyear

  if(sum(sel)<1){print("ERROR!!!!!")}

  if(sum(sel)==1){
    team[i] <- as.character(mydat[sel,"Team"])
    pos[i] <- as.character(mydat[sel,"Position"])
    age[i] <- mydat[sel,"Age"]+1
    ab[i] <- mydat[sel,"AB"]
    hand[i] <- as.character(mydat[sel,"Hand"])
    park[i] <- as.character(mydat[sel,"Park"])
  }

  if(sum(sel)>1){
    team[i] <- as.character(mydat[sel,"Team"][sum(sel)])
    pos[i] <- as.character(mydat[sel,"Position"][sum(sel)])
    age[i] <- mydat[sel,"Age"][sum(sel)]+1
    ab[i] <- sum(mydat[sel,"AB"])
    hand[i] <- as.character(mydat[sel,"Hand"][sum(sel)])
    park[i] <- as.character(mydat[sel,"Park"][sum(sel)])
  }
}

# Create pred data:
preddat1 <- data.frame(Name=myid, Team=team, Position=pos, Year=myyear+1, Age=age, AB=ab, Hand=hand, Park=park, Homeruns=-1)
preddat1 <- rbind(BaseballData, preddat1)
SortKey <- paste( preddat1[, "Name"], preddat1[, "Year"], sep="_" )
SortOrder <- order(SortKey)
SortKey <- SortKey[SortOrder]
preddat1 <- preddat1[SortOrder, ]
preddat <- preddat1
save(preddat, file="preddat_lastab.Rdata")


preddat2 <- data.frame(Name=myid, Team=team, Position=pos, Year=myyear+1, Age=age, AB=500, Hand=hand, Park=park, Homeruns=-1)
preddat2 <- rbind(BaseballData, preddat2)
SortKey <- paste( preddat2[, "Name"], preddat2[, "Year"], sep="_" )
SortOrder <- order(SortKey)
SortKey <- SortKey[SortOrder]
preddat2 <- preddat2[SortOrder, ]
preddat <- preddat2
save(preddat, file="preddat_500ab.Rdata")


