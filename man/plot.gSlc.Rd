\name{plot.gSlc}
\alias{plot.gSlc}
\title{Plot smooth function components of \code{gSlc()} fits}
\description{Smooth function components of generalized additive mixed model fits obtained via \code{gSlc} are plotted.
}
%
\usage{
\method{plot}{gSlc}(x,gridSize=401,colour = TRUE,responseScale = FALSE,
                    rug = TRUE,rugColour="dodgerblue",curveColour = "darkgreen",
                    varBandPolygon = TRUE,varBandColour = "palegreen",
                    xlab = NULL,ylab = NULL,bty = "l",cex.axis = 1,
                    cex.lab = 1,...)
}
%
\arguments{
  \item{x}{\code{gSlc()} fit object.}
  \item{gridSize}{Number of grid points used in graphical display of smooth function fits.} 
  \item{colour}{Boolean flag:\cr
                TRUE = produce colour plots\cr
                FALSE = produce black and white plots.}
  \item{responseScale}{Boolean flag:\cr
        TRUE  = the smooth function fits are plotted on the response scale\cr
        FALSE = the smooth function fits are plotted on the link scale (the default).}
  \item{rug}{Boolean flag:\cr
        TRUE  = add rug graphics to the base of each smooth function plot showing the predictor data (the default),\cr
        FALSE = do not add rug graphs.}
  \item{rugColour}{colour of the rug graphics. The default value is "dodgerblue".}
  \item{curveColour}{colour of the curves in the smooth function display. The default value is "darkgreen".}
  \item{varBandPolygon}{Boolean flag:\cr
                      TRUE = display the variability band as a polygon (the default),\cr
                      FALSE = display the variability band using dashed curves.}
  \item{varBandColour}{colour of the variability band polygon in the smooth function display. The default value is "palegreen".}
  \item{xlab}{optional argument: character string vector for horizontal labels for smooth function plots.}
  \item{ylab}{optional argument: character string vector for vertical labels for smooth function plots.}
  \item{bty}{character string which specifies the type of box which is drawn about plots. See \code{help(par)} for details. The default value is "l".}
  \item{cex.axis}{positive number specifying the factor by which numbers along the axes are expanded.}
  \item{cex.lab}{positive number specifying the factor by which characters in the axis labels are expanded.}
  \item{...}{place-holder for other graphic parameters.}
}
\details{For each smooth function component of the generalized additive mixed model specified in the call to \code{gSlc} the pointwise posterior mean is plotted along with a shaded polygon corresponding to pointwise 95\% credible sets.
}
\references{Pham, T. and Wand, M.P. (2018). Generalized additive mixed model analysis via \code{gammSlice}. \emph{Australian and New Zealand Journal of Statistics}, \bold{60}, 279-300.
}
\author{
Tung Pham \email{tungstats@gmail.com} and Matt Wand \email{matt.wand@uts.edu.au}.
}
%
\seealso{\code{\link{gSlc}}, \code{\link{summary.gSlc}}}
\examples{
library(gammSlice)
set.seed(53902)
n <- 400 ; x <- runif(n)
fTrue <- function(x) return(cos(4*pi*x) + 2*x - 1)
mu <- exp(fTrue(x)) ; y <- rpois(n,mu)
fit <- gSlc(y~s(x),family="poisson",control=gSlc.control(nBurn=200,nKept=200,nThin=1,msgCode=0))
plot(fit)
plot(fit,responseScale=TRUE,rug=FALSE)
points(x,y,col="dodgerblue")
}
% 

