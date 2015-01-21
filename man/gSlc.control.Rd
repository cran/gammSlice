\name{gSlc.control}
%
\alias{gSlc.control}
%
\title{Modification of the \code{gSlc} defaults}
%
\description{The fitting of a generalized additive mixed model
  in \code{gammSlice} via the core function \code{gSlc} entails
  various fitting options to be set to default values. The
  \code{gScl.control} function supports modification of these
  options.}
%
\usage{
gSlc.control(nBurnin = 5000,nIter = 5000,nThin = 5,fixedEffPriorVar = 1e10,
                                                   sdPriorScale = 1e5)
}
%
\arguments{
\item{nBurnin}{The length of the burnin. The first \code{nBurnin} Markov
chain Monte Carlo samples are discarded. Its default value of
\code{nBurnin} is 5000.
}
\item{nIter}{The number of retained Markov chain Monte Carlo samples
after the burnin period. The default value of \code{nIter} is 5000.
}
\item{nThin}{Thinning factor applied to the retained Markov chain
  Monte Carlo samples. Setting \code{nIter} to be an integer greater
  than 1 results in every \code{nIter}th value in the
  post-burnin samples being retained. The final Markov chain Monte Carlo sample
  size is an integer close to \code{nIter} divided by \code{nIter}. 
  The default value of \code{nThin} is 5.
}
\item{fixedEffPriorVar}{The variance in the independent zero mean Normal
  priors of the fixed effect parameters after the predictor data have
  been transformed to the interval [0,1]. The default value of
  \code{fixedEffPriorVar} is 1e10.
}
\item{sdPriorScale}{The scale parameter in the Half Cauchy priors on
  standard deviation parameters after the predictor data have
  been transformed to the interval [0,1]. The default value of
  \code{sdPriorScale} is 1e5.
}
}
\references{
  Pham, T. and Wand, M.P. (2012).\cr
  Generalized additive mixed model analysis via \code{gammSlice}.\cr
  Submitted.

  Zhao, Y., Staudenmayer, J., Coull, B.A. and Wand, M.P. (2003).\cr
  General design Bayesian generalized linear mixed models.\cr
  \emph{Statistical Science}, \bold{21}, 35-51.
}
\author{
  Tung Pham \email{tung.pham@epfl.ch} and Matt Wand \email{matt.wand@uts.edu.au}.
}

\seealso{\code{\link{gSlc}}
}
\examples{
\dontrun{
dat <- gSlcSim(eg=3, num.rep = 10, num.gps = 100, family = "binomial",
                random.factor = FALSE)

fitDefault <- gSlc(y~s(x1), family = "binomial", data = data)

fitMyPriors <- gSlc(y~s(x1), family = "binomial", data = data,
                   control = gSlc.control(fixedEffPriorVar=1e13,
                   sdPriorScale=1e3))

fitBigMCMC <- gSlc(y~s(x1), family = "binomial", data = data,
                   control = gSlc.control(nBurnin=10000,nIter=8000,
                   nThin=10))
}
}
\keyword{models}
\keyword{regression}
