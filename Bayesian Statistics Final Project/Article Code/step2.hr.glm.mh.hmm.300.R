##########################################################
# Baseball Hitting Model
# MCMC Implementation
#
# Author: Blake McShane
# Date: August 2006
# (Revisions through November 2006)
#
# R Libraries Required:
# splines
#
#
# Purpose: This code is designed to draw points from the
# posterior distribution of the model specified in
# "Bayesball".  For more information, please consult that
# paper or the author.
#
##########################################################

##########################################################
# Set Filename, Working Directory, and Formula
#filename = "R.hr.glm.mh.hmm.nohand"
#setwd("/Users/blakemcshane/Documents/Academic Research/Baseball/trajectory/SplineModels_19902010_2011")
#model_formula <- paste( "YMatrix ~ BaseballData[, 'Park'] + BaseballData[, 'Position']",
#                        " + BaseballData[, 'Position'] * Elite ",
#                        " + BaseballData[, 'Position'] * bs(BaseballData[, 'Age'], df = 4)", sep="")

filename = "R.hr.glm.mh.hmm"
lahman <- "lahman591-csv"
begyear <- 1990
endyear <- 2011
predyear <- 2012

setwd(paste("/Users/blakemcshane/Documents/Academic Research/Baseball/trajectory/SplineModels_",
	begyear,endyear,"_",predyear,sep=""))
model_formula <- paste( "YMatrix ~ BaseballData[, 'Park']*BaseballData[,'Hand'] + BaseballData[, 'Position']",
                        " + BaseballData[, 'Position'] * Elite ",
                        " + BaseballData[, 'Position'] * bs(BaseballData[, 'Age'], df = 4)", sep="")

# End: Filename, Working Directory, and Formula
##########################################################

##########################################################
# Call External Libraries
library(splines)
# End: Call External Libraries
##########################################################

##########################################################
# Functions written for this program:

cnt_fcn <- function(n, e){

  # This function takes as inputs a list of names and a 0/1
  # "elite" vector corresponding to each name.  It then computes
  # the number of transitions from state to state for each name
  # and sums across names.

  n00 <- 0
  n01 <- 0
  n10 <- 0
  n11 <- 0
  un <- sort(unique(n))

  for(i in 1:length(un)){
    tmpe <- e[n==un[i]]

    if( length(tmpe)>=2 ){
      for(j in 2:length(tmpe)) {
        if( tmpe[j-1]==0 & tmpe[j]==0 ) { n00 <- n00+1 }
        if( tmpe[j-1]==0 & tmpe[j]==1 ) { n01 <- n01+1 }
        if( tmpe[j-1]==1 & tmpe[j]==0 ) { n10 <- n10+1 }
        if( tmpe[j-1]==1 & tmpe[j]==1 ) { n11 <- n11+1 }
      }
    }

  }

  n <- c(n00, n01, n10, n11)
  names(n) <- c("n00", "n01", "n10", "n11")
  return(n)
}

MLR.pred <- function(X, beta, cat){

  # This function takes as inputs a multinomial logistic model
  # matrix, a vector of coefficients, and a number of categories
  # and returns the predicted probabilities for each category.
  # It was written to be quicker than the canned routine because
  # MCMC requires these predictions often.
  
  k <- cat-1
  n <- dim(X)[1] / k
  p <- dim(X)[2] / k

  tmp <- matrix( X %*% beta, n, k, byrow=TRUE )
  phat <- exp( tmp ) / ( 1 + rowSums(exp(tmp)) )
  phat <- cbind( phat, 1 - rowSums(phat) )
  phat
}

MLR.LL <- function(y, X, beta, cat){

  # This function takes as inputs a multinomial logistic model
  # matrix, a vector of coefficients, and a number of categories
  # and returns the log likelihood for the regression.
  # It was written to be quicker than the canned routine because
  # MCMC requires this calculation often.

  phat <- MLR.pred(X, beta, cat)
  LL <- sum(y * log(phat))
  LL
}

MLR.modelmatrix <- function(M, cat){

  # This function takes as inputs a *binomial* logistic model
  # matrix a number of categories and returns a *multinomial*
  # logisitic model matrix with the appropriate number of categories
  # It was written to be quicker than the canned routine because
  # the canned routine performs maximum likelihood estimation before
  # returning the model matrix.  This is often unnecessary for 
  # our MCMC procedure.

  kk <- cat - 1
  n <- dim(M)[1]
  p <- dim(M)[2]

  #Build Basic Matrix
  X <- matrix(NA, kk * dim(M)[1], kk * dim(M)[2])
  for(i in 1:dim(M)[2]){
    for(j in 1:kk){
      tmp <- NULL
      for(k in 1:kk){

        if (j==k){                    # This fills the appropriate values into
          tmp <- c(tmp, M[, i])       # the MLR model matrix.  The appropriate
        }else{                        # values are either the values from the
          tmp <- c(tmp, rep(0, n))    # BLR model matrix or 0s.
        }

      }
      X[, (i-1) * kk + j] <- tmp
    }
  }

  #Organize Columns
  for(i in 1:dim(X)[2]){
    tmp <- matrix( X[, i], n, kk )
    X[, i] <- as.vector( t(tmp) )
  }
  X
}

BLR.modelmatrix <- function(Elite){

  # This function takes as an input an elite status matrix
  # and returns a binomial logistic regression matrix.

  return(model.matrix( as.formula(model_formula) ))
}


NuPosteriorDraw <- function(beta_a, beta_b, beta_c, beta_d,Data,E){

  # This function takes as an input an set of parameters for the prior
  # distribution on the elite status transition probability matrix,
  # a data matrix with positions and names, and an elite path matrix. It
  # returns a draw from the posterior distribution for the transition
  # probability matrix.

  P <- length(levels(Data[,"Position"]))
  nu <- list(NULL)
  nu[[1]] <- rep(NA, P)     #nu00
  nu[[2]] <- rep(NA, P)     #nu11

 for(j in 1:P){
    TempSelection <- Data[, "Position"]==PositionList[j] 
    n <- cnt_fcn( Data[TempSelection, "Name"], E[TempSelection] )
    n00 <- n[1]
    n01 <- n[2]
    n10 <- n[3]
    n11 <- n[4]

    # The posterior distribution is conjugate beta.  More details
    # can be found on page XXX of the paper.
    nu[[1]][j] <- rbeta( 1, beta_a + n00, beta_b + n01 )
    nu[[2]][j] <- rbeta( 1, beta_c + n11, beta_d + n10 )
  }

  return(nu)
}

# End: Functions written for this program
##########################################################


##########################################################
# Load and Organize Data:

# Read in the raw data from BBRF file.
# Sort it and place in a dataframe.
#data <- read.table("trajdata.1990-2006.txt", sep=",", header=T)
#data <- read.csv("lahman58-csv/Batting.csv")
data <- read.csv(paste(lahman,"/Batting.csv",sep=""))

# Augment with processed data:
load("procdat.Rdata")
sel <- rowSums(is.na(procdat)) > 0
procdat <- procdat[!sel,]
sel <- data[,"playerID"] %in% procdat[,"playerID"]
data <- data[sel,]
sel <- data[,"yearID"] >= min(procdat[,"yearID"])  &  data[,"yearID"] <= max(procdat[,"yearID"])
data <- data[sel,]

# Get indices and augment data with processed data:
tmpidx <- rep(NA, dim(data)[1])
for(i in 1:length(tmpidx)){
  sel <- as.character(data[i,"playerID"])==as.character(procdat[,"playerID"]) & data[i,"yearID"]==procdat[,"yearID"]
  if(sum(sel)>0){ tmpidx[i] <- which(sel) }
  if(i%%1000==0){print(paste("Complete: ", i, "/", length(tmpidx), sep=""))}
}
table(is.na(tmpidx))
data <- data[!is.na(tmpidx),]
tmpidx <- tmpidx[!is.na(tmpidx)]
data <- data.frame(data, pos=procdat[tmpidx,"Position"], age=procdat[tmpidx,"Age"], hand=procdat[tmpidx,"Hand"])

# Sort data:
SortKey <- paste( data[, "playerID"], data[, "yearID"], data[,"stint"], sep="_" )
SortOrder <- order(SortKey)
SortKey <- SortKey[SortOrder]
data <- data[SortOrder, ]

BaseballData <- data.frame(   Name = data[, "playerID"], 
                              Team = data[, "teamID"], 
                              Position = data[, "pos"], 
                              Year = data[, "yearID"], 
                              Age = data[, "age"], 
                              AB = data[, "AB"],
                              Hand = data[,"hand"],
                              Park = data[,"teamID"],        # <----- park is temporarily given by team
                              Homeruns = data[, "HR"] )

# Change B Hand to R:
sel <- BaseballData[,"Hand"]=="B"
BaseballData[sel,"Hand"] <- "R"

# End: Load and Organize Data
##########################################################

##########################################################
# Toy Data for Debugging:
#BaseballData <- BaseballData[1:20,]
#BaseballData[,"Name"] <- c( rep("a",5), rep("b",5), rep("c",5), rep("d",5) )
#BaseballData[,"Position"] <- as.factor(c( rep("1B",5), rep("1B",5), rep("2B",5), rep("2B",5) ))
#BaseballData[,"Year"] <- rep( 2000:2004, 4)
#BaseballData[,"Age"] <- rep(1:10,2)
#BaseballData[,"Park"] <- as.factor(c( rep("P1",5), rep("P2",5), rep("P1",5), rep("P2",5) ))
#BaseballData[,"AB"] <- rep(100,20)
#BaseballData[,"Homeruns"] <- c(5,5,5,10,10,7,7,7,7,12, 6,6,6,15,15, 3,3,10,3,3)
# End: Toy Data for Debugging
##########################################################

##########################################################
# Check for Missing Data and Exclude as appropriate.
# Also, select data to be used in analysis.

# Eliminate Missing Data:
for(i in 1:dim(BaseballData)[2]){
  TempSelection <- is.na(BaseballData[, i])
  BaseballData <- BaseballData[!TempSelection, ]
}

# Eliminate Players with 0 At Bats:
TempSelection <- BaseballData[, "AB"]<=0
BaseballData <- BaseballData[!TempSelection, ]

# Eliminate Data before begyear:
TempSelection <- BaseballData[, "Year"]>=begyear
BaseballData <- BaseballData[TempSelection, ]

# Eliminate Data after endyear:
TempSelection <- BaseballData[, "Year"]<=endyear
BaseballData <- BaseballData[TempSelection, ]

# Eliminate Pitchers since we do not care about their hitting:
TempSelection <- BaseballData[, "Position"]=="P"
BaseballData <- BaseballData[!TempSelection, ]

# Reset the levels for the categorical variables:
# We must do this so the MLR will not complain about "empty" X categories
for(i in 1:dim(BaseballData)[2]){
  if( is.factor(BaseballData[, i]) ){
    BaseballData[, i] <- as.factor(as.character(BaseballData[, i]))
  }
}

# End: Check for Missing Data and Exclude as appropriate.
#      Also, select data to be used in analysis.
##########################################################

##########################################################
# Eliminate players who lack a single 1HR/40AB and 300 AB season.

sel <- BaseballData[,"Homeruns"]/BaseballData[,"AB"] >= 1/40 & BaseballData[,"AB"] >= 300
tmpname <- sort(unique( BaseballData[sel,"Name"] ))
sel <- BaseballData[,"Name"] %in% tmpname
BaseballData <- BaseballData[sel,]

# End: Eliminate players who lack a single 1HR/40AB and 300 AB season.
##########################################################

##########################################################
# Save final dataset.
save(BaseballData, file="mcmcdat.Rdata")
# End: Save final dataset.
##########################################################


##########################################################
# Create Y matrix:

YMatrix <- cbind( BaseballData[, "Homeruns"],
                  BaseballData[,"AB"]-BaseballData[,"Homeruns"] )
colnames(YMatrix) <- c("Homeruns", "Failures")

# End: Create Y matrix:
##########################################################

##########################################################
# Set up variables for MCMC:

# Parameter Values:
NumDraws <- 6000     # The number of MCMC draws from the posterior.
beta_a <- 1         # These are the hyper-parameters for the elite transition
beta_b <- 1         # probability matrix.
beta_c <- 1         # We use a Uniform (ie, Beta(1,1)) prior.
beta_d <- 1
gamma2 <- 10000       # This is the variance for the MLR coefficient prior.
                      # We use a N(0,10000) as our relatively diffuse prior.
P <- length(levels(BaseballData[,"Position"]))

# Create and Initialize the elite matrix:
# We start by setting the top 10% of players from each position as elite
# for each of the 4 hitting events.
Elite <- rep(NA, dim(YMatrix)[1])
PositionList <- levels(BaseballData[, "Position"])
c <- rep(NA, P)
names(c) <-  PositionList
for(j in 1:P){
  tmp <- BaseballData[, "Homeruns"] / BaseballData[, "AB"]
  TempSelection <- BaseballData[, "Position"]==PositionList[j]
  c[j] <- quantile(tmp[TempSelection], .9)
  Elite[TempSelection & tmp>c[j]] <- 1 * 
       (TempSelection & tmp>c[j])[TempSelection & tmp>c[j]]
}
Elite[is.na(Elite)] <- 0

# Create v00, v01, v10, v11:
# These are the sampled elite probability transitions for each position.
v00 <- matrix(NA, NumDraws, P)   # Each "v" is a NumDraws x 9 matrix.
v01 <- matrix(NA, NumDraws, P)   
v10 <- matrix(NA, NumDraws, P)   
v11 <- matrix(NA, NumDraws, P)
colnames(v00) <- PositionList
colnames(v01) <- PositionList
colnames(v10) <- PositionList
colnames(v11) <- PositionList

# Initialize v00, v01, v10, v11:
# In order to do this, we impose a Beta(1,1) prior on the transition probabilities
# and compute the posterior based on the empirical number of "counts".  We draw
# from the posterior using a function above and then put the values into the
# appropriate variable.
tmp <- NuPosteriorDraw(beta_a,beta_b,beta_c,beta_d,BaseballData,Elite)
v00[1, ] <- tmp[[1]]
v01[1, ] <- 1 - tmp[[1]]
v11[1, ] <- tmp[[2]]
v10[1, ] <- 1 - tmp[[2]]

# Create and Initialize betas:
# The starting values for the betas will be the MLE.
lm <- glm( as.formula(model_formula), family="binomial")
Betas <- matrix(NA, NumDraws, length(coef(lm)))
colnames(Betas) <- names(lm$coef)
Betas[1,] <- lm$coef

# Initialize Model Matrix:
XX <- model.matrix(lm)                                                       

# Create Elite/Non-elite Model Matrices:
# For sampling elite paths, we will need to calculate MLR fitted probabilities
# for when the player is elite and when he is not.  Hence, we need two extra
# model vectors, one with all player-seasons set to elite and one with none
# set to elite.
Elite1 <- rep(1,length(Elite))
Elite0 <- rep(0,length(Elite))
XX1 <- BLR.modelmatrix(Elite1)
XX0 <- BLR.modelmatrix(Elite0)

# Finally, we need a list of unique player names:                                                                                                                 
UniqueNames <- sort(unique(BaseballData[, "Name"]))

# End: Set up variables for MCMC:
##########################################################

##########################################################
# MCMC:
date()
for(i in 2:NumDraws){

  # Draw Elite Paths: Setup
  Elite <- rep(NA, dim(YMatrix)[1])            # We must discard the Elite Matrix on
                                               # each draw due to memory constraints.
  tmp1 <- XX1 %*% Betas[i-1,]                  # This gives the MLR phats for the
  tmp0 <- XX0 %*% Betas[i-1,]                  # player-season when elite and non-elite
  p1   <- YMatrix[,1] * tmp1 - rowSums(YMatrix) * log(1+exp(tmp1))
  p0   <- YMatrix[,1] * tmp0 - rowSums(YMatrix) * log(1+exp(tmp0))
  if(is.na(sum(tmp1))){ print(c("tmp1", i)); stop() }
  if(is.na(sum(tmp0))){ print(c("tmp0", i)); stop() }

  # Draw Elite Paths: Player Loop
  # See the paper page XXX for documentation of the Hidden Markov Model.
  # Variable names used below correspond exactly to variables names
  # used in the paper or to those used above.
  # The comments below refer to sections of the paper documentation.
  for(j in 1:length(UniqueNames)) {

    # Set-up: Retrieve pertinent variables for a given player.
    TempSelection <- BaseballData[, "Name"]==UniqueNames[j]
    tmpindex <- (1:dim(BaseballData)[1])[TempSelection]
    tmpT <- length(tmpindex)
    #tmppos <- rep(BaseballData[TempSelection, "Position"][tmpT], tmpT)
    tmppos <- BaseballData[TempSelection, "Position"]
    tmpp0 <- p0[tmpindex]
    tmpp1 <- p1[tmpindex]
    tmpr <- rep(NA, tmpT)

    # Declare forward probability vector for each state and sample elite vector
    a0 <- rep(NA, length(tmpindex))
    a1 <- rep(NA, length(tmpindex))
    e.samp <- rep(NA, length(tmpindex))

    # Forward Probabilites: Initialization
    # This corresponds to "Step 1" of the HMM Forward Probability section.   
    # It uses the obvious result that joint probability of the initial state
    # emitting the first observation is simply the product of the initial
    # distribution vector and the emission probabilities
    loga0 <- tmpp0[1] + log(v00[i-1, tmppos[1]])
    loga1 <- tmpp1[1] + log(v01[i-1, tmppos[1]])
    logM <- max(loga0, loga1)
    a0[1] <- exp(loga0 - logM)
    a1[1] <- exp(loga1 - logM)

    # Forward Probabilites: Induction
    # This corresponds to "Step 2" of the HMM Forward Probability section.
    # This step allows us to calculate the joint probability of:
    #  i)  Observations 1,2,...,t
    #  ii) State at time t is State i
    # given the results for t-1.
    if( tmpT>=2 ){
      for(k in 2:tmpT) {
        loga0 <- tmpp0[k] +
                      log( a0[k-1]*v00[i-1, tmppos[k]]
                           + a1[k-1]*v10[i-1, tmppos[k]] )
        loga1 <- tmpp1[k] +
                      log( a0[k-1]*v01[i-1, tmppos[k]] 
                           + a1[k-1]*v11[i-1, tmppos[k]] )
        logM <- max(loga0, loga1)
        a0[k] <- exp(loga0 - logM)
        a1[k] <- exp(loga1 - logM)
      }
    }

    # Sample: Final Elite State
    # This corresponds to "Step 1" of the HMM Elite Sampling section.   
    # Given the results above, we know the joint probability of:
    #  i)  Observations 1,2,...,T
    #  ii) State at time T is State i
    # Hence, we draw a random number to sample the elite state for the
    # final year.
    r <- a1[tmpT] / (a0[tmpT] + a1[tmpT])
    u <- runif(1, 0, 1)
    if( u<=r ) { e.samp[tmpT] <- 1 }
    if( u>r ) { e.samp[tmpT] <- 0 }
    tmpr[tmpT] <- r

    # Sample: Reverse Induction
    # This corresponds to "Step 2" of the HMM Elite Sampling section.   
    # Conditional upon the draw we made above, we can perform a simple reverse
    # induction and then draw a random number to sample the elite state for
    # year t.
    if( tmpT>=2 ){
      for(k in (tmpT-1):1) {

        if( e.samp[k+1] == 0 ) {
          tmp_v0k <- v00[i-1, tmppos[k]]
          tmp_v1k <- v10[i-1, tmppos[k]]
        }
        if( e.samp[k+1] == 1 ) {
          tmp_v0k <- v01[i-1, tmppos[k]]
          tmp_v1k <- v11[i-1, tmppos[k]]
        }

        r <- a1[k] * tmp_v1k / (a0[k] * tmp_v0k + a1[k] * tmp_v1k)
        u <- runif(1, 0, 1)
        if( u<=r ) { e.samp[k] <- 1 }
        if( u>r ) { e.samp[k] <- 0 }
        tmpr[k] <- r
      }
    }
    Elite[tmpindex] <- e.samp
  }

  # Draw Nus:
  # We sample the Nus in the same fashion that we initialized.  Namely,
  # we calculate the sampled number of elite transitions and make a draw
  # from the posterior,
  tmp <- NuPosteriorDraw(beta_a,beta_b,beta_c,beta_d,BaseballData,Elite)
  v00[i, ] <- tmp[[1]]
  v01[i, ] <- 1 - tmp[[1]]
  v11[i, ] <- tmp[[2]]
  v10[i, ] <- 1 - tmp[[2]]
  if(is.na(sum(v00[i,]))){ print(c("v00", i)); stop() }
  if(is.na(sum(v11[i,]))){ print(c("v11", i)); stop() }

  # In order to sample from the posterior distribution for the regression
  # coefficients, we use Metropolis-Hastings.  Again, variable names correspond
  # to those set out in the paper.
  
  # Metropolis-Hastings Step 1: Proposal:
  # Our proposal distribution is a normal distribution centered at the previous
  # value with a standard deviation of one.

  # THESE FIRST FEW LINES WERE ADDED IN 2012 SO THAT ALL ELITE ARE NOT ZERO OR ONE
  # FOR A GIVEN POSITION.  THERE ARE FEW POS=RF IN THE NEWER DATA SO IT IS KLUGE.
  tmptab <- table(BaseballData[,"Position"], Elite)
  tmptab <- tmptab/rowSums(tmptab)
  if( any(tmptab==0) ){
		tmpidx <- which(tmptab==0, arr.ind=TRUE)
		for(alpha in 1:nrow(tmpidx)){
			sel <- BaseballData[,"Position"]==rownames(tmptab)[tmpidx[alpha,1]]
			tmpy <- YMatrix[,1] / (YMatrix[,1]+YMatrix[,2])
			rr <- range(tmpy[sel])
			Elite[sel & tmpy==rr[1]] <- 0
			Elite[sel & tmpy==rr[2]] <- 1
		}
  }

  tmplm <- glm( as.formula(model_formula), family="binomial")  
  beta_prop <- tmplm[[1]]
  beta_sd <- summary(tmplm)$coefficients[,2]
  if(is.na(sum(beta_prop))){ print(c("beta_prop", i)); stop() }

  # Rewrite the model matrix based on the new Elite draw:
  XX <- BLR.modelmatrix(Elite)

  # Metropolis-Hastings Step 2: Accept or Reject proposal:
  for(j in 1:dim(Betas)[2]) {
    tmpbeta <- Betas[i-1, ]
    beta.prop <- rnorm(1, beta_prop[j], beta_sd[j])
    #beta.prop <- rnorm(1, Betas[i-1,j], 1)
    tmpbeta[j] <- beta.prop
    Xold <- XX %*% Betas[i-1,]
    Xnew <- XX %*% tmpbeta

    pold <- (YMatrix[,1] * Xold 
             - rowSums(YMatrix) * log(1 + exp(Xold))
             - (1/2) * log(gamma2)
             - (1/2) * (Betas[i-1, j]^2) / gamma2)
    pnew <- (YMatrix[,1] * Xnew
             - rowSums(YMatrix) * log(1 + exp(Xnew))
             - (1/2) * log(gamma2)
             - (1/2) * (beta.prop^2) / gamma2)
    gold <- log( dnorm(Betas[i-1, j], beta_prop[j], beta_sd[j]) )
    gnew <- log( dnorm(beta.prop,     beta_prop[j], beta_sd[j]) )
    #gold <- log( dnorm(Betas[i-1, j], Betas[i-1, j], 1) )
    #gnew <- log( dnorm(beta.prop,     Betas[i-1, j], 1) )

    logr <- sum(pnew)-sum(pold)+gold-gnew
    logu <- log(runif(1, 0, 1))
    if( logu<=logr ) { Betas[i, j] <- beta.prop }
    if( logr<logu ) { Betas[i, j] <- Betas[i-1, j] }
    if(is.na(Betas[i,j])){ print(c("Betas", i, j)); stop() }
  }

  # Print to the screen every 100 iterations:
  if(i%%100==0){ print(i) }
}
date()
# End: MCMC:
##########################################################

##########################################################
# Write output graphs and tables:

#v00 pdf:
name <- paste(filename, "v00.pdf", sep="_")
pdf(name, width=8, height=8)
par(mfrow=c(2, 2))
for(j in 1:9) {
  plot( 1:NumDraws, v00[, j], type="l", xlab="Iterations", ylim=c(0,1),
        ylab=paste("v00 - ", colnames(YMatrix)[i], " ", PositionList[j], sep="")) 
}
dev.off()

#v11 pdf:
name <- paste(filename, "v11.pdf", sep="_")
pdf(name, width=8, height=8)
par(mfrow=c(2, 2))
for(j in 1:9) {
  plot( 1:NumDraws, v11[, j], type="l", xlab="Iterations", ylim=c(0,1),
        ylab=paste("v00 - ", colnames(YMatrix)[i], " ", PositionList[j], sep="")) 
}
dev.off()

#BETAs pdf:
name <- paste(filename, "beta.pdf", sep="_")
pdf(name, width=8, height=8)
par(mfrow=c(3, 3))
for(i in 1:dim(Betas)[2] ) {
  plot( 1:NumDraws, Betas[, i], type="l",
        xlab="Iterations", ylab=colnames(Betas)[i])
}
dev.off()

# Print percent of unique betas:
ubetas <- rep(NA, dim(Betas)[2])
for( i in 1:dim(Betas)[2] ) {
  ubetas[i] <- length(unique(Betas[, i]))/NumDraws
}
print( ubetas )

# Write tables:
name <- paste(filename, "output.Rdata", sep="_")
save(v00, v11, Betas, file=name)

# End: Write output graphs and tables:
##########################################################
