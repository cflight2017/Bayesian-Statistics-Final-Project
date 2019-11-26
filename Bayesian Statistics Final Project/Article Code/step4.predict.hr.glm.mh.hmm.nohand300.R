##########################################################
# Baseball Hitting Model
# MCMC Prediction
#
# Author: Blake McShane
# Date: September 2007
#
# R Libraries Required:
# splines
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
#SortKey <- paste( data[, "Name"], data[, "Year"], sep="_" )
#SortOrder <- order(SortKey)
#SortKey <- SortKey[SortOrder]
#data <- data[SortOrder, ]

#BaseballData <- data.frame(   Name = data[, "id"], 
#                              Team = data[, "team"], 
#                              Position = data[, "pos"], 
#                              Year = data[, "year"], 
#                              Age = data[, "age"], 
#                              AB = data[, "AB"],
#                              Hand = data[,"hand"],
#                              Park = data[,"park"],    
#                              Homeruns = data[, "HR"] )

load("mcmcdat.Rdata")
mcmcData <- BaseballData
load("preddat_lastab.Rdata")
BaseballData <- preddat

# Change B Hand to R:
sel <- BaseballData[,"Hand"]=="B"
BaseballData[sel,"Hand"] <- "R"

# Load old runs:
load(paste("output/",filename,"_output.Rdata",sep=""))      #*# <--- CHANGE to output1, output2, etc.
nn <- dim(v00)[1]

# Burn old runs:
burn <- 1000
v00 <- v00[(burn+1):nn,]
v11 <- v11[(burn+1):nn,]
Betas <- Betas[(burn+1):nn,]

# Thin old runs:
thin <- 10
idx <- seq(1, dim(v00)[1], thin)
v00 <- v00[idx,]
v11 <- v11[idx,]
Betas <- Betas[idx,]
v01 <- 1-v00
v10 <- 1-v11
kk <- dim(v00)[1]

# Kill Parks That Are New for Predict Year
#Parks <- read.table( "parks.unique.1970-2006.txt", header=T, sep=",")
#TempSelection <- Parks[,"firstyear"]==2006
#BadParks <- Parks[TempSelection,"parkid"]
#for(i in length(BadParks)){
#  TempSelection <- as.character(BaseballData[,"Park"])==as.character(BadParks[i])
#  BaseballData <- BaseballData[!TempSelection,]
#}

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

# Eliminate Data before 1990:
TempSelection <- BaseballData[, "Year"]>=begyear
BaseballData <- BaseballData[TempSelection, ]

### Eliminate Data after 2010:
### Note: 2010 is the last year for which we have data so we will have
### one full out of sample year
###TempSelection <- BaseballData[, "Year"]<=2010
###BaseballData <- BaseballData[TempSelection, ]

# Eliminate Pitchers since we do not care about their hitting:
TempSelection <- BaseballData[, "Position"]=="P"
BaseballData <- BaseballData[!TempSelection, ]

# Kill Players Who Do Not Play in Predict Year
UniqueNames <- sort(unique(BaseballData[,"Name"]))
Keep <- rep(1, dim(BaseballData)[1])
for(i in 1:length(UniqueNames)){
  TempSelection1 <- BaseballData[,"Name"]==UniqueNames[i] & BaseballData[,"Year"]==max(BaseballData[,"Year"])
  TempSelection2 <- BaseballData[,"Name"]==UniqueNames[i]
  if(sum(TempSelection1)==0){
    Keep <- Keep - TempSelection2
  }
}
Keep <- as.logical(Keep)
BaseballData <- BaseballData[Keep,]

# Kill Players Who Only Play in Predict Year
UniqueNames <- sort(unique(BaseballData[,"Name"]))
Keep <- rep(1, dim(BaseballData)[1])
for(i in 1:length(UniqueNames)){
  TempSelection1 <- BaseballData[,"Name"]==UniqueNames[i] & BaseballData[,"Year"]==predyear
  TempSelection2 <- BaseballData[,"Name"]==UniqueNames[i]
  if(all(TempSelection1==TempSelection2)){
    Keep <- Keep - TempSelection1
  }
}
Keep <- as.logical(Keep)
BaseballData <- BaseballData[Keep,]

###### KILL PLAYERS WHO ARE NOT IN INTERSECTION OF MARCEL PECOTA & HIGHHR
###tmp <- read.table("names_new_intersection300.txt", header=T)
###tmp <- tmp[,1]
###sel <- BaseballData[,"Name"] %in% tmp
###BaseballData <- BaseballData[sel,]

# Reset the levels for the categorical variables:
# We must do this so the MLR will not complain about "empty" X categories
for(i in 1:dim(BaseballData)[2]){
  if( is.factor(BaseballData[, i]) ){
    BaseballData[, i] <- factor(as.character(BaseballData[, i]), levels=as.character(sort(unique(mcmcData[,i]))))
  }
}

#### Fix Park levels so the columns of Betas matches the columns of the XX model matrix
###p <- read.table("parksXX.txt", header=T)[,1]
###BaseballData[,"Park"] <- factor( as.character(BaseballData[,"Park"]), levels=as.character(p) )

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
# Create Y matrix:

YMatrix <- cbind( BaseballData[, "Homeruns"],
                  BaseballData[,"AB"]-BaseballData[,"Homeruns"] )
colnames(YMatrix) <- c("Homeruns", "Failures")

# End: Create Y matrix:
##########################################################

##########################################################
# Set up variables for MCMC:

# Parameter Values:
NumDraws <- kk     # The number of MCMC draws from the posterior.
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

# Create and Initialize linear model:
# The starting values for the betas will be the MLE.
#lm <- glm( as.formula(model_formula), family="binomial")

# Initialize Model Matrix:
XX <- model.matrix(as.formula(model_formula)) #model.matrix(lm)                                                       

# Create Elite/Non-elite Model Matrices:
# For sampling elite paths, we will need to calculate MLR fitted probabilities
# for when the player is elite and when he is not.  Hence, we need two extra
# model vectors, one with all player-seasons set to elite and one with none
# set to elite.
Elite1 <- rep(1,length(Elite))
Elite0 <- rep(0,length(Elite))
XX1 <- BLR.modelmatrix(Elite1)
XX0 <- BLR.modelmatrix(Elite0)

# We need a list of unique player names:                                                                                                                 
UniqueNames <- sort(unique(BaseballData[, "Name"]))

# Finally, we need to store prediction variables:
TempSelection <- BaseballData[,"Year"]==max(BaseballData[,"Year"])
PredictionNames <- BaseballData[TempSelection,"Name"]
Predictions <- matrix(NA, length(PredictionNames), kk)
PredictionsRate <- matrix(NA, length(PredictionNames), kk)

# End: Set up variables for MCMC:
##########################################################

##########################################################
# MCMC:
sum(Betas)
date()
for(i in 1:NumDraws){

  # Draw Elite Paths: Setup
  Elite <- rep(NA, dim(YMatrix)[1])          # We must discard the Elite Matrix on
                                             # each draw due to memory constraints.
  tmp1 <- XX1 %*% Betas[i,]                  # This gives the MLR phats for the
  tmp0 <- XX0 %*% Betas[i,]                  # player-season when elite and non-elite
  p1   <- YMatrix[,1] * tmp1 - rowSums(YMatrix) * log(1+exp(tmp1))
  p0   <- YMatrix[,1] * tmp0 - rowSums(YMatrix) * log(1+exp(tmp0))

  # Draw Elite Paths: Player Loop
  # See the paper page XXX for documentation of the Hidden Markov Model.
  # Variable names used below correspond exactly to variables names
  # used in the paper or to those used above.
  # The comments below refer to sections of the paper documentation.
  for(j in 1:length(UniqueNames)) {

    # Set-up: Retrieve pertinent variables for a given player.
    TempSelection <- BaseballData[, "Name"]==UniqueNames[j] & BaseballData[,"Year"]<=(max(BaseballData[,"Year"])-1)
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
    loga0 <- tmpp0[1] + log(v00[i, tmppos[1]])
    loga1 <- tmpp1[1] + log(v01[i, tmppos[1]])
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
                      log( a0[k-1]*v00[i, tmppos[k]]
                           + a1[k-1]*v10[i, tmppos[k]] )
        loga1 <- tmpp1[k] +
                      log( a0[k-1]*v01[i, tmppos[k]] 
                           + a1[k-1]*v11[i, tmppos[k]] )
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
          tmp_v0k <- v00[i, tmppos[k]]
          tmp_v1k <- v10[i, tmppos[k]]
        }
        if( e.samp[k+1] == 1 ) {
          tmp_v0k <- v01[i, tmppos[k]]
          tmp_v1k <- v11[i, tmppos[k]]
        }

        r <- a1[k] * tmp_v1k / (a0[k] * tmp_v0k + a1[k] * tmp_v1k)
        u <- runif(1, 0, 1)
        if( u<=r ) { e.samp[k] <- 1 }
        if( u>r ) { e.samp[k] <- 0 }
        tmpr[k] <- r
      }
    }
    Elite[tmpindex] <- e.samp
    
    # Now predict Elite forward to Predict Year:
    tmpindex <- BaseballData[, "Name"]==UniqueNames[j] & BaseballData[,"Year"]==max(BaseballData[,"Year"])
    tmppos <- tmppos[tmpT]
    tmpelite <- e.samp[tmpT]
    if( tmpelite == 0){
      tmp_vk0 <- v00[i, tmppos]
    }else{
      tmp_vk0 <- v10[i, tmppos]
    }
    u <- runif(1)
    if(u<=tmp_vk0){
      Elite[tmpindex] <- 0
    }else{
      Elite[tmpindex] <- 1
    }

  }

  # Rewrite the model matrix based on the new Elite draw and Predict:
  XX <- BLR.modelmatrix(Elite)
  LogOdds <- XX %*% Betas[i,]
  p <- exp(LogOdds) / (1+exp(LogOdds))
  Sample <- rbinom( dim(BaseballData)[1], BaseballData[,"AB"], p)
  TempSelection <- BaseballData[,"Year"]==max(BaseballData[,"Year"])
  Predictions[,i] <- Sample[TempSelection]
  PredictionsRate[,i] <- p[TempSelection]

  if(i%%25==0){ print(i) }
}
date()
# End: MCMC:
##########################################################

TempSelection <- BaseballData[,"Year"]==max(BaseballData[,"Year"])
BaseballData_ForecastAssumptions <- BaseballData[TempSelection,]
rownames(Predictions) <- as.character(PredictionNames)
rownames(PredictionsRate) <- as.character(PredictionNames)
save(BaseballData_ForecastAssumptions, Predictions, PredictionsRate, file= paste("output/",filename,"__predictions_lastab.Rdata",sep=""))        #*# <--- CHANGE to lastab1, lastab2, etc.