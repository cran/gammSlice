\name{summary.gSlc}
%
\alias{summary.gSlc}
%
\title{Summary of the generalized additive mixed model fit
       produced by \code{gSlc}}
%
\description{A graphical table showing,for key model parameters,
  the Markov chain Monte Carlo (MCMC)
  samples, diagnostic plots and meaningful summaries.
}
%
\usage{
\method{summary}{gSlc}(object,pages = 1,...)
}
%
\arguments{
   \item{object}{A \code{gSlc} object.}
   \item{pages}{The number of pages used to display summary. 
        The default value of \code{pages} is 1.}
   \item{\dots}{Other parameters.}
}
%
\details{
  The columns of the graphical table are:\cr
  
   1. parameter name,\cr
   2. trace plot of the MCMC sample,\cr
   3. plot of MCMC sample against its lag 1 sample,\cr
   4. sample autocorrelation function,\cr
   5. kernel density estimate of posterior density function,\cr
   6. numerical summaries of posterior density function.
 }
%
\references{
  Pham, T. and Wand, M.P. (2012).\cr
  Generalized additive mixed model analysis via \code{gammSlice}.\cr
  Submitted.
}
%
\author{
Tung Pham \email{tung.pham@epfl.ch} and Matt Wand \email{matt.wand@uts.edu.au}. 
}
%
\seealso{
\code{\link{gSlc}}, \code{\link{plot.gSlc}}
}
%
\examples{
\dontrun{
dat1 <- gSlcSim(eg=3, num.rep = 10, num.gps = 100, family = "binomial",
                random.factor = FALSE)
fit1 <- gSlc(y~s(x1),family = "binomial", data = dat1)
summary(fit1)

dat2 <- gSlcSim(eg=6,family = "poisson", random.factor = TRUE)
fit2 <- gSlc(y~x1 + x2 + s(x3) + s(x4), random = list(idnum=~1),
             family = "poisson", data = dat2)
summary(fit2,pages=2)
}
}
%
\keyword{models}
\keyword{regression}
