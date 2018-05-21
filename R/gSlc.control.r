########## R function: gSlc.control ##########

# For control of gSlc().

# Last changed: 03 MAY 2018

gSlc.control <- function(nBurn = 5000,nKept = 5000,nThin = 5,
                         fixedEffPriorVar = 1e10,sdPriorScale = 1e5,
                         numBasis = NULL,preTransfData = TRUE,msgCode = 1)
{  
   # Guard against non-integer and negative inputs of MCMC sample size inputs:

   nBurn <- round(nBurn)
   nKept <- round(nKept)
   nThin <- round(nThin)

   if (nBurn < 0) 
   {
      warning("The inputted number of burnin iterations is negative.\n
               The default value of 5000 was used instead.")
      nBurn <- 5000
   }

   if (nKept < 0) 
   {
      warning("The inputted number of kept iterations is negative.\n
               The default value of 5000 was used instead.")
      nKept <- 5000
   }

   if (nThin < 0) 
   {
      warning("The inputted thinning factor is negative.\n
               The default value of 5 was used instead.")
      nThin <- 5
   }

   if (nThin > nKept) 
   {
      warning("The inputting thinning factor is larger than the\n
               number of kept iterations. The default value of\n
               these inputs were used instead.")
      nThin <- 5
      nKept <- 5000
   }

   # Guard against the negative inputs of the hyperparameters:

   if (fixedEffPriorVar < 0) 
   {
      warning("The inputted fixed effect prior variance is negative.\n
               The default value of 1e10 was used instead.")
      fixedEffPriorVar <- 1e10
   }

   if (sdPriorScale < 0) 
   {
      warning("The inputted standard deviation prior scale is negative.\n
               The default value of 1e5 was used instead.")
      sdPriorScale <- 1e5
   }
 
   # Guard against non-integer and low number of basis functions:
  
   if (!is.null(numBasis))
   {
      numBasis <- round(numBasis)
      if (min(numBasis)<3)
      {
          warning("An inputted number of basis functions was less than 3.
                   The default value was used instead.") 
          numBasis <- NULL
      }
   }

   # Guard against illegal message code value:

   if (!any(msgCode==0:3))
   {
      warning("The inputted message code was neither 0,1,2 nor 3.\n
               The default value of 1 was used instead.") 
      msgCode <- 1 
   }   

   outList <- list(nKept=nKept,nBurn = nBurn, nThin = nThin,fixedEffPriorVar = fixedEffPriorVar, 
                   sdPriorScale = sdPriorScale,numBasis = numBasis,preTransfData = preTransfData,
                   msgCode = msgCode)
   return(outList)
}

############ End of gSlc.control ############

