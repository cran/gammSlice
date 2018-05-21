########## R function: plot.gSlc ##########

# Produce plots from a gSlc() fit object.

# Last changed: 06 MAY 2018

plot.gSlc <- function(x,gridSize=401,colour = TRUE,responseScale = FALSE,
                      rug = TRUE,rugColour="dodgerblue",curveColour = "darkgreen",
                      varBandPolygon = TRUE,varBandColour = "palegreen",
                      xlab = NULL,ylab = NULL,bty = "l",cex.axis = 1,
                      cex.lab = 1,...)
{
   fitObject <- x

   # Check legality of secondary inputs:

   if (!any(varBandPolygon == c(TRUE,FALSE)))
   {
      warning("varBandPolygon must be either TRUE or FALSE.\n
               The default value of TRUE was used instead.")
      varBandPolygon <- TRUE
   }

   # Set or overwrite colours for black and white version:

   if (!colour)
   {
      curveColour <- "black"
      varBandColour <- "grey70"
      rugColour <- "black"
   } 
 
   # Extract model type:

   modelType <- fitObject$modelType

   # Create flag for axis label extraction:

   extractXlabs <- is.null(xlab)
   extractYlabs <- is.null(ylab)

   # Extract pre-transformation flag:
 
   preTransfData <- fitObject$preTransfData

   # Obtain `XlinPreds' and `XsplPreds' since they help
   # determine plotting grids:

   XlinPreds <- fitObject$XlinPreds   
   XsplPreds <- fitObject$XsplPreds

   if (is.null(XlinPreds)) numLinCompons <- 0 
   if (!is.null(XlinPreds)) numLinCompons <- ncol(XlinPreds)

   if (is.null(XsplPreds)) numSplCompons <- 0 
   if (!is.null(XsplPreds)) numSplCompons <- ncol(XsplPreds)

   # Save original versions for back-transformation purposes:   
 
   XlinPredsOrig <- XlinPreds
   XsplPredsOrig <- XsplPreds

   if (preTransfData)
   {

      # Replace 'XlinPreds' and 'XsplPreds' by unit interval transformed
      # values to match the situation in gSlc():
 
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

   # Print non-plot message for models without smooth function components:

   if (any(modelType==c("linOnly","linOnlyMix")))
      cat("There are no smooth function components to plot.\n")

   # Now deal with the case of smooth function components being present:

   if (!any(modelType==c("linOnly","linOnlyMix")))
   {
      # Extract dimension type variables:

      idnumPresent <- FALSE
      if (any(modelType==c("linAddMix","pureAddMix")))
         idnumPresent <- TRUE

      ncZspl <- fitObject$ncZspl
      
      # Extract `nu' chain:
 
      nuMCMC <- t(fitObject$nu)
      if (!idnumPresent) nuMCMCnonRE <- nuMCMC
      if (idnumPresent) 
      {
         nonREinds <- 1:(1 + numLinCompons + numSplCompons + sum(ncZspl))
         nuMCMCnonRE <- nuMCMC[nonREinds,] 
      }

      # Extract linear and spline predictors matrices and spline basis 
      # function  information:

      range.x.list <- fitObject$range.x.list
      intKnots.list <- fitObject$intKnots.list

      # Construct full design matrix corresponding to average predictor values:

      CsplMeang <- rep(1,gridSize)    # For beta0.
       
      if (numLinCompons>0)
         for (jLin in 1:numLinCompons)
            CsplMeang <- cbind(CsplMeang,rep(mean(XlinPreds[,jLin]),gridSize))

      for (jSpl in 1:numSplCompons)
      {      
         xgCurr <- rep(mean(XsplPreds[,jSpl]),gridSize)
         CsplMeang <- cbind(CsplMeang,xgCurr)
      }
      for (jSpl in 1:numSplCompons)
      {      
         xgCurr <- rep(mean(XsplPreds[,jSpl]),gridSize)
         ZgCurr <- ZOSull(xgCurr,intKnots=intKnots.list[[jSpl]],
                          range.x=range.x.list[[jSpl]])
         CsplMeang <- cbind(CsplMeang,ZgCurr)
      }
      
      # Obtain plot for each smooth function with all other predictors
      # set to their average values:

      indsZstt <- 1 + numLinCompons + numSplCompons + 1 

      Csplg <- CsplMeang

      for (jSpl in 1:numSplCompons)
      {
         # Set up the current gridwise design matrices:

         xCurr <- XsplPreds[,jSpl]
         xLow <- min(xCurr) 
         xUpp <- max(xCurr) 
         
         xgCurr <- seq(xLow,xUpp,length=gridSize)
         
         ZgCurr <- ZOSull(xgCurr,intKnots=intKnots.list[[jSpl]],
                          range.x=range.x.list[[jSpl]])

         # Insert current gridwise design matrix into the average
         # value design matrix:

         indXspl <- 1 + numLinCompons + jSpl 
         Csplg[,indXspl] <- xgCurr

         indsZend <- indsZstt + ncZspl[jSpl] - 1

         Csplg[,indsZstt:indsZend] <- ZgCurr

         # Obtain current curve MCMC matrix:

         curveMCMCg <- crossprod(t(Csplg),nuMCMCnonRE)

         indsZstt <- indsZend + 1

         if (responseScale)
         {
            if (fitObject$family=="binomial") curveMCMCg <- plogis(curveMCMCg)
            if (fitObject$family=="poisson")  curveMCMCg <- exp(curveMCMCg)
         }

         # Do current plot: 

         meang <- apply(curveMCMCg,1,mean)
         lowerg <- apply(curveMCMCg,1,quantile,0.025)
         upperg <- apply(curveMCMCg,1,quantile,0.975)

         if (extractXlabs) xlabVal <- fitObject$splPredNames[jSpl]
         if (!extractXlabs) xlabVal <- xlab[jSpl]

         if (extractYlabs) 
         {
            if (fitObject$family=="binomial")       
            {
               if (responseScale) effectWording <- "probability"
               if (!responseScale) effectWording <- "logit(probability)"
            }
            if (fitObject$family=="poisson")       
            {
               if (responseScale) effectWording <- "mean response"
               if (!responseScale) effectWording <- "log(mean response)"
            }
            ylabVal <- paste("effect of ",xlabVal," on ",effectWording,sep="")
         }
         if (!extractYlabs) ylabVal <- ylab[jSpl]

         # Obtain versions of the horizontal data and plotting grid 
         # corresponding the original predictor units:

         if (!preTransfData)
         {
            xOrigCurr <- xCurr
            xgOrigCurr <- xgCurr
         }
         if (preTransfData)
         {
            minCurr <- min(XsplPredsOrig[,jSpl])  ;    maxCurr <- max(XsplPredsOrig[,jSpl])
            denCurr <- maxCurr - minCurr 
            xOrigCurr <- minCurr + xCurr*denCurr
            xgOrigCurr <- minCurr + xgCurr*denCurr
         }

         plot(0,type="n",xlim=range(xgOrigCurr),ylim=range(c(lowerg,upperg)),
              xlab = xlabVal, ylab = ylabVal, bty = bty,cex.axis = cex.axis,
              cex.lab = cex.lab)

         if (varBandPolygon)
            polygon(c(xgOrigCurr,rev(xgOrigCurr)),c(lowerg,rev(upperg)),border=FALSE,col=varBandColour)

         if (!varBandPolygon)
         { 
            varBandColour <- curveColour                          
            lines(xgOrigCurr,lowerg,col=curveColour,lwd=2,lty=2)
            lines(xgOrigCurr,upperg,col=curveColour,lwd=2,lty=2)
         }

         lines(xgOrigCurr,meang,col=curveColour,lwd=2)

         if (rug) rug(xOrigCurr,col=rugColour)

         if (jSpl<numSplCompons) readline("Hit Enter to continue.\n")

         Csplg <- CsplMeang
      }
   }
   invisible()
}

############ End of plot.gSlc #################
