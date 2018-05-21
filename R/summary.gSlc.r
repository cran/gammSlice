########## R function: summary.gSlc ##########

# Produce a summary from a gSlc() fit object.

# Last changed: 12 MAY 2018

summary.gSlc <-function(object,colour=TRUE,paletteNumber=1,...) 
{
   fitObject <- object

   # Check legality of seconday inputs:

   if (!any(colour==c(TRUE,FALSE))) stop("colour must be TRUE or FALSE.\n")
   if (!any(paletteNumber==c(1,2))) stop("paletteNumber must be 1 or 2.\n")

   # Extract model type:

   modelType <- fitObject$modelType

   # Extract dimension type variables:

   idnumPresent <- FALSE
   if (any(modelType==c("linOnlyMix","linAddMix","pureAddMix")))
      idnumPresent <- TRUE

   if (is.null(fitObject$XlinPreds)) numLinCompons <- 0
   if (!is.null(fitObject$XlinPreds)) numLinCompons <- ncol(fitObject$XlinPreds)
   if (is.null(fitObject$XsplPreds)) numSplCompons <- 0
   if (!is.null(fitObject$XsplPreds)) numSplCompons <- ncol(fitObject$XsplPreds)

   ncZspl <- fitObject$ncZspl

   # Extract variable names:

   linPredNames <- fitObject$linPredNames 
   splPredNames <- fitObject$splPredNames

   # Extract `nu' and `sigsq' chains:
 
   nuMCMC <- t(fitObject$nu)
   sigsqMCMC <- t(fitObject$sigmaSquared)

   # Extract `beta' chains if relevant:

   if (any(modelType==c("linOnly","linAdd","linOnlyMix","linAddMix")))
      betaMCMC <- t(nuMCMC[1:(numLinCompons+1),])

   # Extract random effect variance chain if relevant:    

   if (idnumPresent)
      sigsqRanEffMCMC <- sigsqMCMC[nrow(sigsqMCMC),]

   # Extract pre-transformation flag:
 
   preTransfData <- fitObject$preTransfData

   # Perform back-transformation of betaMCMC and sigsqRanEffMCMC if relevant:

   if (preTransfData&any(modelType==c("linOnly","linAdd","linOnlyMix","linAddMix")))
   {
      XlinPreds <- fitObject$XlinPreds

      minVec <- apply(XlinPreds,2,min)   ;   maxVec <- apply(XlinPreds,2,max)
      denVec <- maxVec - minVec

      betaStarMCMC <- betaMCMC

      # First do the adjustment for the non-intercept fixed effect 
      # regression coefficients:

      for (jLin in 2:(1+numLinCompons))
         betaMCMC[,jLin] <- betaMCMC[,jLin]/denVec[jLin-1]

      # Next, do the adjustment for the fixed effect intercept:

      for (jLin in 2:(1+numLinCompons))
         betaMCMC[,1] <- betaMCMC[,1] - betaStarMCMC[,jLin]*minVec[jLin-1]/denVec[jLin-1]
   }

   # Obtain `edf' chains if relevant:

   if (numSplCompons>0)
   {
      # Extract required design matrices:

      XsplPreds <- fitObject$XsplPreds
      Zspl <- fitObject$Zspl
      Cspl <- cbind(XsplPreds,Zspl)

      # Obtain the `wMCMC' matrix:

      etaMCMC <- crossprod(t(Cspl),nuMCMC[1:ncol(Cspl),])

      if (fitObject$family=="binomial") wMCMC <- 0.5/(1 + cosh(etaMCMC))

      if (fitObject$family=="poisson") wMCMC <- exp(etaMCMC)

      # Obtain the `diagonalEDFMAT' matrix:

      nMCMCfinal <- ncol(nuMCMC)
      edfMCMC <- rep(NA,nMCMCfinal)
      Dmat <- diag(c(rep(0,ncol(XsplPreds)),rep(1,ncol(Zspl))))
      diagonalEDFMAT <- matrix(NA,ncol(Cspl),nMCMCfinal)
      for (iMCMC in 1:nMCMCfinal)
      {
         CTWCcurr <- crossprod(Cspl*wMCMC[,iMCMC],Cspl)

         diagonalLambda <- rep(0,numSplCompons)
         for (jSpl in 1:numSplCompons) 
            diagonalLambda <- c(diagonalLambda,(1/sigsqMCMC[jSpl,iMCMC])*rep(1,ncZspl[jSpl]))
      
         diagonalEDFMAT[,iMCMC] <- diag(solve(CTWCcurr + diag(diagonalLambda),CTWCcurr))
      }

      # Extract the edfMCMC matrices:
     
      edfMCMC <- vector("list",numSplCompons)
 
      indsZstt <- numSplCompons + 1 
      for (jSpl in 1:numSplCompons)
      {
         indsZend <- indsZstt + ncZspl[jSpl] - 1
         indsZcurr <- indsZstt:indsZend
         indsCcurr <- c(jSpl,indsZcurr)
         edfMCMC[[jSpl]] <- colSums(diagonalEDFMAT[indsCcurr,])
         indsZstt <- indsZend + 1
      }

      # Make intecept adjustment for the pure nonparametric regression
      # or nonparametric regression with random intecept special case:

      if (any(modelType==c("pureAdd","pureAddMix"))&(numSplCompons==1))
         edfMCMC[[1]] <- edfMCMC[[1]] + 1
   }

   MCMClistLength <- numLinCompons + numSplCompons + as.numeric(idnumPresent)

   beta0included <- any(modelType==c("linOnly","linOnlyMix"))

   if (beta0included)   # beta0 included:
      MCMClistLength <- MCMClistLength + 1

   parNamesVal <- vector("list",MCMClistLength)

   MCMCmat <- NULL
   if (any(modelType==c("linOnly","linAdd","linOnlyMix","linAddMix")))
   {
      if (any(modelType==c("linOnly","linOnlyMix")))    # beta0 included:
      {   
         MCMCmat <- cbind(MCMCmat,betaMCMC)
         parNamesVal[[1]] <- expression(beta[0])
         indStt <- 1
      }

      if (!any(modelType==c("linOnly","linOnlyMix")))   # beta0 excluded:
      {   
         MCMCmat <- cbind(MCMCmat,betaMCMC[,-1])
         indStt <- 0
      }

      if (numLinCompons>0)
         for (jLin in 1:numLinCompons)
            parNamesVal[[indStt+jLin]] <- eval(bquote(expression(beta[.(linPredNames[jLin])])))
   }

   if (numSplCompons>0)
   {
      for (jSpl in 1:numSplCompons)
      {
         MCMCmat <- cbind(MCMCmat,edfMCMC[[jSpl]])
         parNamesVal[[numLinCompons+jSpl]] <- eval(bquote(expression(edf[.(splPredNames[jSpl])])))
      }
   }

   if (idnumPresent)
   {
      MCMCmat <- cbind(MCMCmat,sigsqRanEffMCMC)
      parNamesVal[[as.numeric(beta0included)+numLinCompons+numSplCompons+1]] <- expression(sigma^2)
   }

   # Do MCMC summary plot:

   summMCMC(list(MCMCmat),colourVersion=colour,paletteNum=paletteNumber,parNames=parNamesVal)
}

############ End of summary.gSlc ############
