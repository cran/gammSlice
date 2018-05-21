########### R function: gSlc ##########

# Fit a generalised additive mixed model using
# slice sampling. 

# Last changed: 19 MAY 2018

# Set flag for type of likelihood

gSlc <- function(formula,data = NULL,random = NULL,family,control = gSlc.control())
{
   # Extract control variables:

   nBurn <- control$nBurn
   nKept <- control$nKept
   nThin <- control$nThin
   fixedEffPriorVar <- control$fixedEffPriorVar  
   sdPriorScale <- control$sdPriorScale
   numBasis <- control$numBasis
   preTransfData <- control$preTransfData
   msgCode <- control$msgCode

   # Extract relevant information for the formula:

   formularInformation <- readFormula(formula)

   # Construct matrices and dimension variables.

   responseName <- formularInformation$response   
   linPredNames <- formularInformation$lvar
   splPredNames <- formularInformation$svar

   basisSizes <- formularInformation$nbas
   numLinCompons <- length(linPredNames)
   numSplCompons <- length(splPredNames) 

   if (is.null(data)) y <- get(responseName,pos=1)
   if (!is.null(data)) y <- data[,responseName]
  
   if (length(linPredNames)>0) 
   {
      XlinPreds <- NULL  ;  namesXlinPreds <- NULL
      for (jLin in 1:numLinCompons)
      {
         if (is.null(data))  currPred <- get(linPredNames[jLin],pos=1)
         if (!is.null(data)) currPred <- data[,linPredNames[jLin]]

         if (!is.character(currPred)) 
         {
            XlinPreds <- cbind(XlinPreds,currPred)  
            namesXlinPreds <- c(namesXlinPreds,linPredNames[jLin])
         }
         if (is.character(currPred))
         {
            suqcurrPred <- sort(unique(currPred)) ; baseChar <- suqcurrPred[1] ; 
            numLevs <- length(suqcurrPred)
            if (length(suqcurrPred)==1) stop("Factor has all levels equal.")
            XcurrPred <- NULL ; namesXcurrPred <- NULL
            for (jLev in 2:numLevs)
            {
               XcurrPred <- cbind(XcurrPred,as.numeric(currPred==suqcurrPred[jLev]))
               namesXcurrPred <- c(namesXcurrPred,paste(suqcurrPred[jLev],"VS",suqcurrPred[1],sep=""))
            }
            XcurrPred <- as.data.frame(XcurrPred)
            XlinPreds <- cbind(XlinPreds,XcurrPred)
            namesXlinPreds <- c(namesXlinPreds,namesXcurrPred)
         }
      }
      names(XlinPreds) <- namesXlinPreds
      XlinPreds <- as.matrix(XlinPreds)
      numLinCompons <- ncol(XlinPreds)
      linPredNames <- namesXlinPreds
   }

   if (length(linPredNames)==0) XlinPreds <- NULL

   if (length(splPredNames)>0) 
   {
      XsplPreds <- NULL
      for (jSpl in 1:length(splPredNames))
      {
         if (is.null(data))  XsplPreds <- cbind(XsplPreds,get(splPredNames[jSpl],pos=1))
         if (!is.null(data)) XsplPreds <- cbind(XsplPreds,data[,splPredNames[jSpl]])
      }
      XsplPreds <- as.matrix(XsplPreds)
   }
   if (length(splPredNames)==0) XsplPreds <- NULL
 
   if (preTransfData)
   {

      # Replace 'XlinPreds' and 'XsplPreds' by unit interval transformed
      # values for fitting:

      XlinPredsOrig <- XlinPreds   
      XsplPredsOrig <- XsplPreds

      if (numLinCompons>0)
         for (jLin in 1:numLinCompons)
         {  
            minCurr <- min(XlinPreds[,jLin]) ; maxCurr <- max(XlinPreds[,jLin])
            XlinPreds[,jLin] <- (XlinPreds[,jLin] - minCurr)/(maxCurr - minCurr)
         }

      if (numSplCompons>0)
         for (jSpl in 1:numSplCompons)
         {  
            minCurr <- min(XsplPreds[,jSpl]) ; maxCurr <- max(XsplPreds[,jSpl])
            XsplPreds[,jSpl] <- (XsplPreds[,jSpl] - minCurr)/(maxCurr - minCurr)
         }
   }

   # Construct numBasis vector to be used:

   if (!is.null(numBasis))
   {
      if (length(numBasis)==1) numBasis <- rep(numBasis,numSplCompons)
      if (length(numBasis)!=numSplCompons)
         stop("The length of numSplCompons must match the number of smooth functions.\n")
   }
   if (is.null(numBasis)) 
   {
      numBasis <- rep(NA,numSplCompons)
      for (jSpl in 1:numSplCompons) 
         numBasis[jSpl] <- min(25,round(length(unique(XsplPreds[,jSpl]))/4))
   }
 
   # Find the random variable name:

   if (!is.null(random)) 
   {
      idnumPresent <- TRUE
      if (is.list(random)) 
      {
         r.names <- names(random)
         idnumName <- c(unlist(lapply(random,function(x) all.vars(formula(x)))),r.names)
         if (is.null(data)) idnumOrig <- get(idnumName,pos=1)
         if (!is.null(data)) idnumOrig <- data[,idnumName]
      } 
      else  
         idnumOrig <- NULL

      # Form the natural identification number vector 

      idnum <- match(idnumOrig,unique(idnumOrig))
      numGps <- length(unique(idnum))
   } 
   else 
   {   
      idnum <- NULL
      idnumPresent <- FALSE
   }

   # Build X matrix:

   X <- cbind(1,XlinPreds,XsplPreds)

   # Obtain first component of the `blockInds' list:

   blockInds <- list(1:ncol(X))
   blkIndStt <- ncol(X) + 1

   # Build Zspl matrix:
   
   if (numSplCompons == 0) 
      Zspl <- NULL

   if (numSplCompons == 0) 
   {
      range.x.list <- NULL
      intKnots.list <- NULL
      ncZspl <- NULL
   }
   if (numSplCompons > 0) 
   {
      range.x.list <- vector("list",numSplCompons)
      intKnots.list <- vector("list",numSplCompons)
      ncZspl <- NULL   ;   Zspl <- NULL
      for (jSpl in 1:numSplCompons)
      {
         xCurr <- XsplPreds[,jSpl]
         range.xVal <- c(min(xCurr),max(xCurr))
         numIntKnotsVal <- numBasis[jSpl] - 2
         intKnots <-  as.numeric(quantile(unique(xCurr),seq(0,1,length=numIntKnotsVal+2)
                      [-c(1,numIntKnotsVal+2)]))
         Zcurr <- ZOSull(xCurr,intKnots=intKnots,range.x=range.xVal)
         range.x.list[[jSpl]] <- attr(Zcurr,"range.x")
         intKnots.list[[jSpl]] <- attr(Zcurr,"intKnots")

         ncZspl <- c(ncZspl,ncol(Zcurr)) 
         Zspl <- cbind(Zspl,Zcurr)    

         # Update the blockInds list:

         blkIndEnd <- blkIndStt + ncol(Zcurr) - 1
         blockInds <- c(blockInds,list(blkIndStt:blkIndEnd))
         blkIndStt <- blkIndEnd + 1
      }   
   }
 
   # Append final component to `blockInds' if applicable:

   if (idnumPresent)
      blockInds <-  c(blockInds,list(blkIndStt:(blkIndStt+numGps-1)))

   # Form Cspl matrix:

   Cspl <- cbind(X,Zspl)

   # Obtain Markov chain Monte Carlo samples:

   CsplTy <- as.numeric(crossprod(Cspl,y))
   lenbki <- length(blockInds)
   lenssq <- lenbki - 1
   Ahyp <- rep(sdPriorScale,lenssq)
   ssqbeta <- fixedEffPriorVar
   sttEntry <- function(x) return(x[1])
   endEntry <- function(x) return(x[length(x)])
   ibkstt <- unlist(lapply(blockInds,sttEntry))
   ibkend <- unlist(lapply(blockInds,endEntry))
   nMCMC <- nBurn + nKept
   numObs <- length(y)
   ncCspl <- ncol(Cspl)

   idPres <- idnumPresent
      
   if (idPres==0)
   {
      idnum <- 1   ;   numGps <- 1
      nVec <- 1    ;   lennu <- ncCspl
      idStt <- 1   ;   idEnd <- 1
   } 
   if (idPres==1)
   {
      numGps <- length(unique(idnum))
      nVec <- as.numeric(table(idnum))
      lennu <- ncCspl + numGps
      idEnd <- c((1:numObs)[diff(idnum)==1],numObs)
      idStt <- idEnd - nVec + 1
   } 

   xnu <- matrix(0,lennu,nMCMC)
   sigsq <- matrix(1,lenssq,nMCMC)
   xnucur <- rep(0,lennu)
   lnumn1 <- lennu - 1
   if (family=="binomial") ifam <- 1
   if (family=="poisson")  ifam <- 2
   gparm <- rep(0,2)
   xncnoj <- rep(0,lnumn1)
   buSpl  <- rep(0,ncCspl)
   idxSpl <- rep(0,ncCspl)
   betauO <- rep(0,lennu)
   uSbj  <- rep(0,numGps)
   uSbjO <- rep(0,numGps)
   idxSbj <- rep(0,numGps)
   aux <- matrix(1,lenssq,nMCMC)
   ssqnu <- c(rep(ssqbeta,ncol(X)),rep(1,lennu-ncol(X)))

   cat("\n")
   cat("=====================\n")    
   cat("Welcome to gammSlice.\n")
   cat("=====================\n\n")    
   cat("Rudimentary progress reports available from\n")
   cat("inside the Fortran 77 'engine room' according\n")
   cat("to the control parameter 'msgCode'.\n\n")
   if (msgCode==0)
      cat("Currently msgCode = 0 leading to no reporting.\n")
   if (msgCode==1)
   {
      cat("Currently msgCode = 1 (the default) leading to\n")
      cat("reporting for percentages 1,2,...,10\n")
      cat("and then for percentages 20,30,...,100.\n")
      Sys.sleep(3)
   }
   if (msgCode==2)
   {
      cat("Currently msgCode = 2 leading to reporting\n")
      cat("for percentages 1,2,...,100.\n")
      Sys.sleep(3)
   }
   if (msgCode==3)
   {
      cat("Currently msgCode = 3 leading to reporting\n")
      cat("for percentages 10,20,...,100.\n")
      Sys.sleep(3)
   }

   cat("\n")
   F77obj <- .Fortran("gslcmc",as.double(y),as.double(Cspl),as.double(CsplTy),
                       as.double(Ahyp),as.integer(ibkstt),as.integer(ibkend),
                       as.integer(nMCMC),as.integer(numObs),as.integer(ncCspl),
                       as.integer(idPres),as.integer(idnum),as.integer(numGps),
                       as.integer(idStt),as.integer(idEnd),as.integer(nVec),
                       xnu=as.double(xnu),sigsq=as.double(sigsq),as.double(xnucur),
                       as.integer(lennu),as.integer(lenbki),as.integer(lenssq),
                       as.integer(lnumn1),as.integer(ifam),as.double(gparm),
                       as.double(xncnoj),as.double(buSpl),as.double(betauO),
                       as.integer(idxSpl),as.double(uSbj),as.double(uSbjO),
                       as.integer(idxSbj),as.double(aux),as.double(ssqnu),
                       as.integer(msgCode))
   cat("\n")

   nuMCMC  <- matrix(F77obj$xnu,lennu,nMCMC)
   sigsqMCMC <- matrix(F77obj$sigsq,lenssq,nMCMC)

   # Remove burnin samples:

   nuMCMC <- nuMCMC[,-(1:nBurn)]
   sigsqMCMC <- sigsqMCMC[,-(1:nBurn)]
   if (is.vector(sigsqMCMC)) sigsqMCMC <- t(as.matrix(sigsqMCMC))

   # Perform thinning:

   nKept <- ncol(nuMCMC) ; nKeptPostThin <- floor(nKept/nThin)
   indsFromThinning <- 1 + nThin*(0:nKeptPostThin)
   if (max(indsFromThinning)>nKept)
      indsFromThinning <- indsFromThinning[-length(indsFromThinning)]

   nuMCMC <- nuMCMC[,indsFromThinning]
   sigsqMCMC <- sigsqMCMC[,indsFromThinning]
   if (is.vector(sigsqMCMC)) sigsqMCMC <- t(as.matrix(sigsqMCMC))

   # Prepare nu and sigmasquared MCMC output matrices:

   nuOut <- t(nuMCMC)
   sigmaSquaredOut <- t(sigsqMCMC) 

   # Create betaOut if relevant:

   betaOut <- NULL
   if (numLinCompons>0)
   {
      betaOut <- nuOut[,2:(1+numLinCompons)]
      if (is.vector(betaOut)) betaOut <- t(as.matrix(betaOut))
      if (nrow(betaOut)==1) betaOut <- t(betaOut)
   }

   # Determine whether the model is one of the six following types:
   # "linOnly", "linAdd", "pureAdd", "linOnlyMix", "linAddMix", "pureAddMix".

   if (!idnumPresent)  # Determine if "linOnly", "linAdd" or "pureAdd":
   {
      if (numSplCompons==0) modelType <- "linOnly"
      if ((numLinCompons==0)&(numSplCompons>0)) modelType <- "pureAdd"
      if ((numLinCompons>0)&(numSplCompons>0)) modelType <- "linAdd"
   }
   if (idnumPresent) # Determine if "linOnlyMix", "linAddMix" or "pureAddMix":
   {
      if (numSplCompons==0) modelType <- "linOnlyMix"
      if ((numLinCompons==0)&(numSplCompons>0)) modelType <- "pureAddMix"
      if ((numLinCompons>0)&(numSplCompons>0)) modelType <- "linAddMix"
   }

   if (preTransfData)
   {
      # Reset predictor data to original values:

      XlinPreds <- XlinPredsOrig
      XsplPreds <- XsplPredsOrig
   }

   # Return output list:
  
   outObj <- list(nu=nuOut,beta=betaOut,sigmaSquared=sigmaSquaredOut,
                  y=y,XlinPreds=XlinPreds,linPredNames=linPredNames,
                  XsplPreds=XsplPreds,splPredNames=splPredNames,Zspl=Zspl,ncZspl=ncZspl,
                  range.x.list=range.x.list,intKnots.list=intKnots.list,family=family,
                  modelType=modelType,preTransfData=preTransfData)

   class(outObj) <- "gSlc"

   return(outObj)
}

########## End of gSlc ##########
