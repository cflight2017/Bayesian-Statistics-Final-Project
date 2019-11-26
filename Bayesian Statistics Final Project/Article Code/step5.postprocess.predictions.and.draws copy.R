# Get posterior draws and predictions:
filename = "R.hr.glm.mh.hmm"
begyear <- 1990
endyear <- 2011
predyear <- 2012
lahman <- "lahman591-csv"

setwd(paste("/Users/blakemcshane/Documents/Academic Research/Baseball/trajectory/SplineModels_",
	begyear,endyear,"_",predyear,sep=""))
load("mcmcdat.rdata")
PositionList <- as.character(sort(unique(BaseballData[,"Position"])))
burn <- 1000
thin <- 10
v00list <- list(NULL)
v11list <- list(NULL)
Betaslist <- list(NULL)
Pred500 <- list(NULL)
PredRate500 <- list(NULL)
Predlast <- list(NULL)
PredRatelast <- list(NULL)

# Load MCMC Draws
NUMDRAW <- 5
for(i in 1:NUMDRAW){
	if(i==1){
		load(paste("output/",filename,"_output.Rdata",sep=""))
		nn <- dim(v00)[1]
		idx <- seq(burn+1, nn, by=thin)
	}else{ load(paste("output/",filename,"_output_",i,".Rdata",sep="")) }
	v00list[[i]] <- v00[idx,]
	v11list[[i]] <- v11[idx,]
	Betaslist[[i]] <- Betas[idx,]
}

# Load predictions:
for(i in 1:NUMDRAW){
	if(i==1){ load(paste("output/",filename,"__predictions_lastab.Rdata",sep=""))
	}else{ load(paste("output/",filename,"__predictions_lastab_",i,".Rdata",sep="")) }
	Predlast[[i]] <- Predictions
	PredRatelast[[i]] <- PredictionsRate
}

# Plot v00
NumDraws <- dim(v00list[[1]])[1]
name <- paste("output/fig_final_v00.pdf", sep="_")
pdf(name, width=8, height=8)
par(mfrow=c(2, 2))
for(j in 1:9){
  plot( 1:NumDraws, v00list[[1]][,j], type="l", xlab="Iterations", ylim=c(0,1), col=2,
        ylab=paste("v00 - ", PositionList[j], sep="")) 
  lines( 1:NumDraws, v00list[[2]][,j], col=3 )
  lines( 1:NumDraws, v00list[[3]][,j], col=4 )
}
dev.off()

# Plot v11
name <- paste("output/fig_final_v11.pdf", sep="_")
pdf(name, width=8, height=8)
par(mfrow=c(2, 2))
for(j in 1:9){
  plot( 1:NumDraws, v11list[[1]][,j], type="l", xlab="Iterations", ylim=c(0,1), col=2,
        ylab=paste("v11 - ", PositionList[j], sep="")) 
  lines( 1:NumDraws, v11list[[2]][,j], col=3 )
  lines( 1:NumDraws, v11list[[3]][,j], col=4 )
}
dev.off()

# Plot Betas
name <- paste("output/fig_final_betas.pdf", sep="_")
pdf(name, width=8, height=8)
par(mfrow=c(3,3))
for(i in 1:dim(Betaslist[[1]])[2]){
  m <- min( Betaslist[[1]][,i], Betaslist[[2]][,i], Betaslist[[3]][,i] )
  M <- max( Betaslist[[1]][,i], Betaslist[[2]][,i], Betaslist[[3]][,i] )
  plot( 1:NumDraws, Betaslist[[1]][,i], type="l", xlab="Iterations", ylim=c(m,M), col=2, ylab=colnames(Betaslist[[1]])[i]) 
  lines( 1:NumDraws, Betaslist[[2]][,i], col=3 )
  lines( 1:NumDraws, Betaslist[[3]][,i], col=4 )
}
dev.off()

# Plot Predictions
name <- paste("output/fig_final_pred_lastab.pdf", sep="_")
pdf(name, width=8, height=8)
par(mfrow=c(3,3))
for(i in 1:dim(Predlast[[1]])[1]){
  m <- min( Predlast[[1]][i,], Predlast[[2]][i,], Predlast[[3]][i,] )
  M <- max( Predlast[[1]][i,], Predlast[[2]][i,], Predlast[[3]][i,] )
  plot( 1:NumDraws, Predlast[[1]][i,], type="l", xlab="Iterations", ylim=c(m,M), col=2, ylab=rownames(Predlast[[1]])[i]) 
  lines( 1:NumDraws, Predlast[[2]][i,], col=3 )
  lines( 1:NumDraws, Predlast[[3]][i,], col=4 )
}
dev.off()

# Plot Prediction Rate
name <- paste("output/fig_final_predrate_lastab.pdf", sep="_")
pdf(name, width=8, height=8)
par(mfrow=c(3,3))
for(i in 1:dim(PredRatelast[[1]])[1]){
  m <- min( PredRatelast[[1]][i,], PredRatelast[[2]][i,], PredRatelast[[3]][i,] )
  M <- max( PredRatelast[[1]][i,], PredRatelast[[2]][i,], PredRatelast[[3]][i,] )
  plot( 1:NumDraws, PredRatelast[[1]][i,], type="l", xlab="Iterations", ylim=c(m,M), col=2, ylab=rownames(PredRatelast[[1]])[i]) 
  lines( 1:NumDraws, PredRatelast[[2]][i,], col=3 )
  lines( 1:NumDraws, PredRatelast[[3]][i,], col=4 )
}
dev.off()

outPredlast <- Predlast[[1]]
outPredRatelast <- PredRatelast[[1]]
if(NUMDRAW > 1){
	for(i in 2:NUMDRAW){
		outPredlast <- cbind(outPredlast, Predlast[[i]])
		outPredRatelast <- cbind(outPredRatelast, PredRatelast[[i]])
	}
}
save(BaseballData_ForecastAssumptions, outPredlast, outPredRatelast, file="final_pred.Rdata")




data <- read.csv(paste(lahman,"/Master.csv",sep=""))
sel <- data[,2]==""
data <- data[!sel,]
rownames(data) <- as.character(data[,"playerID"])
myplayers <- rownames(outPredlast)
csvout <- data.frame(ID=myplayers, FirstName=data[myplayers,"nameFirst"], LastName=data[myplayers,"nameLast"],
					 BaseballData_ForecastAssumptions[,c(-1,-9)],
	 				 HR=rowMeans(outPredlast), HRRate=rowMeans(outPredRatelast))
write.csv(csvout, "final_pred.csv", row.names=FALSE, quote=FALSE)