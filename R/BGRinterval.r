########## R function: BGRinterval ##########

# Brooks-Gelman-Rubin "Rhat_interval" diagnostic.

# Last changed: 22 FEB 2018 by M.P. Wand

BGRinterval <- function(xList,sttPos = 5,alpha = 0.05)
{
   chainLen <- length(xList[[1]])
   intervLen <- function(x,alpha = 0.05)
   {
      ans <- quantile(x,(1-alpha/2))-quantile(x,(alpha/2))
      names(ans) <- NULL
      return(ans)
   }
   firsti <- function(x,i)
   {
     if (i>length(x)) stop("i too large for this array") 
     return(x[1:i])
   }
   RhatIntNumer <- rep(NA,(chainLen-sttPos))
   RhatIntDenom <- rep(NA,(chainLen-sttPos))

   for (i in sttPos:chainLen)
   {
      xListCurr <- lapply(xList,firsti,i)
      RhatIntNumer[i-sttPos+1] <- intervLen(unlist(xListCurr))
      RhatIntDenom[i-sttPos+1] <- mean(unlist(lapply(xListCurr,intervLen)))
   }   

   return(list(x = (sttPos:chainLen),numer = RhatIntNumer,
              denom = RhatIntDenom))
}

########## End of BGRinterval ##########

