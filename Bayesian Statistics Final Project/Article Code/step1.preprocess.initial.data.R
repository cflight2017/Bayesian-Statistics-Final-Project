# Get player IDs and restrict years:
begyear <- 1990
endyear <- 2011
predyear <- 2012
lahman <- "lahman591-csv"

setwd(paste("/Users/blakemcshane/Documents/Academic Research/Baseball/trajectory/SplineModels_",
	begyear,endyear,"_",predyear,sep=""))
data <- read.csv(paste(lahman,"/Batting.csv",sep=""))
sel <- data[,"yearID"]>=begyear & data[,"yearID"]<=endyear
data <- data[sel,]
tmp <- unique( data[, c("playerID","yearID")] )
pid <- as.character(tmp[,"playerID"])
yid <- tmp[,"yearID"]

# Get position data and eliminate P and OF:
data <- read.csv(paste(lahman,"/Fielding.csv",sep=""))
#sel <- data[,"POS"]=="OF" | data[,"POS"]=="P"
sel <- data[,"POS"]=="P"
data <- data[!sel,]
sel <- pid %in% as.character(data[,"playerID"])
pid <- pid[sel]
yid <- yid[sel]

# Get position:
pos <- rep(NA, length(pid))
for(i in 1:length(pid)){
  sel <- data[,"playerID"]==pid[i] & data[,"yearID"]==yid[i]
  if(sum(sel)==1){idx <- 1}else{idx <- which(data[sel,"G"]==max(data[sel,"G"]))}
  if(sum(sel)>0){pos[i] <- as.character(data[sel,"POS"][idx])}
  if(i%%1000==0){print(paste("Complete: ", i, "/", length(pid), sep=""))}
}

# Get birthday/age:
data <- read.csv(paste(lahman,"/Master.csv",sep=""))
sel <- data[,"playerID"]==""
data <- data[!sel,]
rownames(data) <- data[,"playerID"]
tmp <- paste( data[pid,"birthYear"], "-", data[pid,"birthMonth"], "-", data[pid,"birthDay"], sep="" )
bday <- as.Date(tmp)
tmp <- as.Date(paste(yid, "-04-01", sep=""))
age <- floor(as.numeric(tmp-bday)/365)

# Get hand:
hand <- as.character(data[pid,"bats"])

# Output data:
procdat <- data.frame(playerID=pid, yearID=yid, Age=age, Position=pos, Hand=hand)
save(procdat, file="procdat.Rdata")
