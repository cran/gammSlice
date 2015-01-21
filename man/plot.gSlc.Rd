\name{plot.gSlc}
%
\alias{plot.gSlc}
%
\title{Plot smooth function components of \code{gSlc} fits}
%
\description{Smooth function components of generalized
  additive mixed model fits obtained via \code{gSlc}
  are plotted.
}
%
\usage{
\method{plot}{gSlc}(x,pages=0,responseScale = FALSE,xlab = NULL,ylab = NULL,main=NULL,
                    bty = NULL,...)
}
%
\arguments{
  \item{x}{\code{gSlc} object.}
  \item{pages}{The number of pages used to display the smooth
    function fits. The default value of \code{pages} is 0.
    Set to 0 to have the routine leave the graphic settings
    as they are.}
  \item{responseScale}{Boolean flag indicating whether or the smooth
    function fits are plotted on the same scale as the response
    variable. If \code{responseScale} is set to \code{FALSE} then
    the the smooth function fits are plotted on the link scale.
    The default value of \code{responseScale} is \code{FALSE}.}
  \item{xlab}{They are used as the x label for plots if supplied. }
  \item{ylab}{They are used as the y labels for plots if supplied.}
  \item{main}{They are used as the title for plots if supplied.}
  \item{bty}{A character string which determined the type of \link{box} which is drawn about plots. See \link{par} for reference.}
  \item{\dots}{Other graphic parameters.}
}
\details{For each smooth function component of the generalized
  additive mixed model specified in the call to \code{gSlc} the
  pointwise posterior mean is plotted along with a shaded polygon
  corresponding to pointwise 95\% credible sets.
}
\references{
  Pham, T. and Wand, M.P. (2012).\cr
  Generalized additive mixed model analysis via \code{gammSlice}.\cr
  Submitted.
}
\author{
  Tung Pham \email{tung.pham@epfl.ch} and Matt Wand
 \email{matt.wand@uts.edu.au}.
}
%
\seealso{\code{\link{gSlc}}, \code{\link{summary.gSlc}}}
\examples{
\dontrun{
dat1 <- gSlcSim(eg=3, num.rep = 10, num.gps = 100, family = "binomial",
                random.factor = FALSE)
fit1 <- gSlc(y~s(x1),family = "binomial", data = dat1)
plot(fit1)
plot(fit1,responseScale=TRUE)

dat2 <- gSlcSim(eg=6,family = "poisson", random.factor = TRUE)
fit2 <- gSlc(y~x1 + x2 + s(x3) + s(x4), random = list(idnum=~1),
             family = "poisson", data = dat2)
plot(fit2,pages=2)
}
}
% 
\keyword{models}
\keyword{regression}
