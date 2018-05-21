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
\item{formula}{Formula describing the generalized additive mixed model.}
\item{data}{Data frame containing the input data. This argument is optional.}
\item{random}{List describing random effects structure.This argument is optional.}
\item{family}{Distribution family of the response variable. The options are "binomial" and "poisson".}
\item{control}{Control options specified by \link{gSlc.control}.
}
}
%
\details{A Bayesian generalized additive mixed model is fitted to the input data according to specified formula. Such models are special cases of the general design generalized linear mixed models of Zhao, Staudenmayer, Coull and Wand (2003). Markov chain Monte Carlo, with slice sampling for the fixed and random effects, is used to obtain samples from the posterior distributions of the model parameters. Full details of the sampling scheme are in the appendix of Pham and Wand (2018). 
}
%
\value{
An object of \code{class} \code{"gScl"}. The functions \code{summary()} and \code{plot()} are used to obtain a summary and plot of the fits. The object is a list with the following components:

\item{nu}{Matrix containing Markov chain Monte Carlo samples of the entire nu=(beta,u) vector. The rows correspond to Markov chain Monte Carlo replicates and the columns correspond to entries of the nu=(beta,u) vector.}
\item{beta}{Matrix containing Markov chain Monte Carlo samples of the beta vector corresponding to the linear components of the model. The rows correspond to Markov chain Monte Carlo replicates and the columns correspond to entries of the beta vector.}
\item{sigmaSquared}{Matrix containing Markov chain Monte Carlo samples of the entire sigma squared vector. The rows correspond to Markov chain Monte Carlo replicates and the columns correspond to entries of the sigmaSquared vector.}
\item{y}{Response data vector.}
\item{XlinPreds}{Matrix containing predictors that are purely linear components of the model.}
\item{linPredNames}{Names of \code{XlinPreds}.}
\item{XsplPreds}{Matrix containing predictors that are penalised spline components of the model.}
\item{splPredNames}{Names of \code{XsplPreds}.}
\item{Zspl}{Horizontal concatenation of each of the spline basis "Z" matrices used for smooth function components.}
\item{ncZspl}{Vector giving the numbers of columns in the horizontal partition of \code{Zspl} corresponding to each smooth function component.}
\item{range.x.list}{List containing values of the \code{range.x} input to the internal \code{ZOSull()} function.}
\item{intKnots.list}{List containing values of the \code{intKnots} input to the internal \code{ZOSull()} function.}
\item{family}{Character string indicating the family of the fitted model; either "binomial" or "poisson".}
\item{modelType}{Charater string indicating the type of model fitted.}

}
%
\references{
  Neal, R.M. (2003).\cr
  Slice sampling (with discussion).\cr
  \emph{The Annals of Statistics}, \bold{31}, 705-767.

  Pham, T. and Wand, M.P. (2018).\cr
  Generalized additive mixed model analysis via \code{gammSlice}.\cr
  \emph{Australian and New Zealand Journal of Statistics}, in press.

  Zhao, Y., Staudenmayer, J., Coull, B.A. and Wand, M.P. (2003).\cr
  General design Bayesian generalized linear mixed models.\cr
  \emph{Statistical Science}, \bold{21}, 35-51.
}
%
\author{
Tung Pham \email{tungstats@gmail.com} and Matt Wand \email{matt.wand@uts.edu.au}
}
%
\seealso{
\code{\link{gSlc.control}}, \code{\link{plot.gSlc}}, \code{\link{summary.gSlc}}
}
%
\examples{

\dontrun{
# Example 1 of Pham & Wand (2018):

set.seed(39402)
m <- 100 ; n <- 2
beta0True <- 0.5 ; betaxTrue <- 1.7 
sigsqTrue <- 0.8 ; idnum <- rep(1:m,each=n)
x <- runif(m*n)
U <- rep(rnorm(m,0,sqrt(sigsqTrue)),each=n)
mu <- 1/(1+exp(-(beta0True+betaxTrue*x+U)))
y <- rbinom((m*n),1,mu)
fit1 <- gSlc(y ~ x,random = list(idnum = ~1),family = "binomial")
summary(fit1)
}

\dontrun{
# Example 2 of Pham & Wand (2018):

set.seed(53902)
n <- 400 ; x <- runif(n)
fTrue <- function(x) return(cos(4*pi*x) + 2*x - 1)
mu <- exp(fTrue(x)) ; y <- rpois(n,mu)
fit2 <- gSlc(y~s(x),family="poisson")
summary(fit2)
plot(fit2)
}

\dontrun{
# Example 3 of Pham & Wand (2018):

set.seed(981127)
n <- 500 ; betax1True <- 0.5;  x1 <- sample(c(0,1),n,replace=TRUE) 
x2 <- runif(n) ; fTrue <- function(x) return(sin(2*pi*x))
mu <- 1/(1+exp(-(betax1True*x1+fTrue(x2)))) ; y <- rpois(n,mu)
y <- rbinom(n,1,mu)
fit3 <- gSlc(y ~ x1 + s(x2),family="binomial")
summary(fit3)
plot(fit3)
}

\dontrun{
# Example 4 of Pham & Wand (2018):

set.seed(2966703)
m <- 100 ; n <- 10;  x1 <- runif(m*n);   x2 <- runif(m*n)
idnum <- rep(1:m,each=n) ; sigsqTrue <- 1
U <- rep(rnorm(m,0,sqrt(sigsqTrue)),each=n)
mu <- exp(U + cos(4*pi*x1) + 2*x1 + sin(2*pi*x2^2)) ; y <- rpois(m*n,mu)
fit4 <- gSlc(y ~ s(x1) + s(x2),random = list(idnum=~1),family = "poisson")
summary(fit4)
plot(fit4)
}
}

