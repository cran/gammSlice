\name{gSlc.control}
\alias{gSlc.control}
\title{Controlling generalized additive mixed model fitting via slice sampling}
\description{Function for optional use in calls to \code{gSlc()} to control Markov chain Monte Carlo sample sizes values and other specifications for slice sampling-based fitting of generalized additive mixed models.}
%
\usage{
gSlc.control(nBurn=5000,nKept=5000,nThin=5,fixedEffPriorVar=1e10,
             sdPriorScale=1e5,numBasis=NULL,preTransfData=TRUE,msgCode=1)
}
%
\arguments{
\item{nBurn}{The length of the Markov chain Monte Carlo burnin. The first \code{nBurnin} Markov chain Monte Carlo samples are discarded. The default value of \code{nBurnin} is 5000.}
\item{nKept}{The number of kept Markov chain Monte Carlo samples after the burnin period. The default value of \code{nKept} is 5000.}
\item{nThin}{Thinning factor applied to the retained Markov chain Monte Carlo samples. Setting \code{nThin} to be an intege greater than 1 results in every \code{nThin}th value in the post-burnin samples being retained. The final Markov chain Monte Carlo sample size is an integer close to \code{nIter} divided by \code{nIter}. The default value of \code{nThin} is 5.}
\item{fixedEffPriorVar}{The variance in the independent zero mean Normal priors of the fixed effect parameters after the data of each predictor have been transformed to the interval [0,1]. The default value of \code{fixedEffPriorVar} is 1e10.
}
\item{sdPriorScale}{The scale parameter in the Half Cauchy priors on standard deviation parameters after the data of each predictor have been transformed to the interval [0,1]. The default value of \code{sdPriorScale} is 1e5.
}
\item{numBasis}{Vector of positive integers specifying the number of spline basis functions to be used for each smooth function component.}
\item{preTransfData}{Boolean flag:\cr
                     TRUE = pre-transform each of the predictors to unit interval for Bayesian analysis with the priors specified by \code{fixedEffPriorVar} and \code{sdPriorScale} (the default),\cr
                     FALSE = do not perform any pre-transformation of the predictors. 
}
\item{msgCode}{A code for specification of the nature of messages printed concerning progress of the Markov chain Monte Carlo sampling:\cr
     0 = no messages printed,\cr
     1 = percentages 1,2,...,10 and then 20,30,...,100\cr 
         (the default),\cr
     2 = percentages 1,2,...,100,\cr
     3 = percentages 10,20,...,100.
}
}
\references{
Pham, T. and Wand, M.P. (2018). Generalized additive mixed model analysis via \code{gammSlice}. \emph{Australian and New Zealand Journal of Statistics}, \bold{60}, 279-300.

Zhao, Y., Staudenmayer, J., Coull, B.A. and Wand, M.P. (2003). General design Bayesian generalized linear mixed models. \emph{Statistical Science}, \bold{21}, 35-51.
}
\author{
Tung Pham \email{tungstats@gmail.com} and Matt Wand \email{matt.wand@uts.edu.au}.
}
\seealso{\code{gSlc}
}
\examples{
\dontrun{
library(gammSlice)
set.seed(39402) ; m <- 100 ; n <- 2
beta0True <- 0.5 ; betaxTrue <- 1.7 ; sigsqTrue <- 0.8
idnum <- rep(1:m,each=n) ; x <- runif(m*n)
U <- rep(rnorm(m,0,sqrt(sigsqTrue)),each=n)
mu <- 1/(1+exp(-(beta0True+betaxTrue*x+U)))
y <- rbinom((m*n),1,mu)
fit <- gSlc(y ~ x,random = list(idnum = ~1),family = "binomial")
summary(fit)

# Illustration of user-specified priors:

fitMyPriors <- gSlc(y ~ x,random = list(idnum = ~1), 
                    family = "binomial", 
                    control = gSlc.control(fixedEffPriorVar=1e13, 
                                           sdPriorScale=1e3))
summary(fitMyPriors)

# Illustration of specification of larger Markov chain Monte Carlo samples:

fitBigMCMC <- gSlc(y ~ x,random = list(idnum = ~1),family = "binomial",
                   control = gSlc.control(nBurn=10000,nKept=8000,nThin=10))
summary(fitBigMCMC)
}
}

