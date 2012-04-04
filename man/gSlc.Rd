\name{gSlc}
%
\alias{gSlc}
%
\title{
Generalized additive mixed model analysis via slice sampling
}
%
\description{
  Use slice sampling-based Markov chain Monte Carlo
  to fit a generalized additive mixed model.
}
%
\usage{
gSlc(formula, data = NULL, random = NULL, family, control = gSlc.control())
}
%
\arguments{
  \item{formula}{
Formula describing the generalized additive mixed model.}
  \item{data}{Data frame containing the input data.}
  \item{random}{List describing random effects structure.
               This argument is optional.}
\item{family}{Distribution family of the response variable.
Options are "binomial" and "poisson".   
}
  \item{control}{Control options specified by \link{gSlc.control}.
}
}
%
\details{A Bayesian generalized additive mixed model is
  fitted to the input data according to specified formula.
  Such models are special cases of the general design
  generalized linear mixed models of Zhao, Staudenmayer,
  Coull and Wand (2003). Markov chain Monte Carlo, with slice
  sampling for the fixed and random effects, is used to obtain samples
  from the posterior distributions of the model parameters.
  Full details of the sampling scheme are in the appendix
  of Pham and Wand (2012). 
}
%
\value{
  \item{nu}{Matrix containing the MCMC samples for the combined fixed
    effects and random effects vectors. Each column of \code{nu} is
    a separate MCMC sample.}
  \item{beta}{Matrix containing the MCMC samples for the fixed effects
    vector.}
  \item{u}{Matrix containing the MCMC samples for the random effects
    vector. If the model contains smooth function components then
    \code{u} includes both random intercept and spline coefficient
    MCMC samples.
   }
  \item{sigmaSquared}{Matrix contain of variances.}
  \item{scaledData}{The scaled data set was used to fit in.}
  \item{formulaInfor}{Information obtained from the formula.}
  \item{timeTaken}{Time in seconds taken by the MCMC sampling.}
  \item{Xmin}{The minimum values of each predictor variable.}
  \item{Xmax}{The maximum values of each predictor variable.}
  \item{Xrange}{The difference between Xmax and Xmin.}
}
%
\references{
  Neal, R.M. (2003).\cr
  Slice sampling (with discussion).\cr
  \emph{The Annals of Statistics}, \bold{31}, 705-767.

  Pham, T. and Wand, M.P. (2012).\cr
  Generalized additive mixed model analysis via \code{gammSlice}.\cr
  Submitted.

  Zhao, Y., Staudenmayer, J., Coull, B.A. and Wand, M.P. (2003).\cr
  General design Bayesian generalized linear mixed models.\cr
  \emph{Statistical Science}, \bold{21}, 35-51.
}
%
\author{
Tung Pham \email{tung.pham@epfl.ch} and Matt Wand \email{matt.wand@uts.edu.au}
}
%
\seealso{
\code{\link{gSlc.control}}, \code{\link{plot.gSlc}}, \code{\link{summary.gSlc}}
}
%
\examples{
\dontrun{
library(mgcv)
dat0 <- gamSim(eg=1, n=500, scale = 0.2, dist = "poisson")
fit0 <- gSlc( y~s(x0) + s(x1) + s(x2) + s(x3), family = "poisson", data = dat0)
plot(fit0,pages = 1)
summary(fit0)

dat1 <- gamSim(eg=6, n = 400,scale = 0.1, dist = "poisson")
fit1 <- gSlc(y ~ s(x0) + s(x1) + s(x2) + s(x3), family = "poisson",
             data = dat1, random = list(fac=~1))
plot(fit1,pages=2)
summary(fit1)

dat2 <- gSlcSim(eg = 2, numGrp = 200, family = "poisson",
                randomFactor = FALSE)
fit2 <- gSlc(y~x1 + x2, family = "poisson", data = dat2)
summary(fit2)

dat3 <- gSlcSim(eg = 3,numGrp = 1000, family = "binomial",
                randomFactor = FALSE)
fit3 <- gSlc(y~s(x1),family = "binomial", data = dat3)
plot(fit3)
summary(fit3)

fit3a <- gSlc(y~s(x1,nBasis=10),family = "binomial",
                 data = dat3)
plot(fit3a)
summary(fit3a)


dat4 <- gSlcSim(eg = 4,  numGrp = 400, family = "poisson",
               randomFactor = FALSE)
fit4 <- gSlc(y~x1 + s(x2), family = "poisson", data = dat4)
plot(fit4)
summary(fit4)

dat5 <- gSlcSim(eg=6,family = "poisson", randomFactor = TRUE)
fit5 <- gSlc(y~x1 + x2 + s(x3) + s(x4), random = list(idnum=~1),
             family = "poisson", data = dat5)
plot(fit5)
summary(fit5)
}
}
\keyword{models}
\keyword{regression}
