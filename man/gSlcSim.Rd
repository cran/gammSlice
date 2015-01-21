\name{gSlcSim}
%
\alias{gSlcSim}
%
\title{
Simulate data from a generalized additive mixed model
}
\description{Six types of generalized additive mixed model
  structures can be simulated. The user has the option
  to specify the response distribution and sample sizes.
}
\usage{
gSlcSim(eg, family, numRep = 10, numGrp = 100, randomFactor)
}
%
\arguments{
  \item{eg}{An integer between 1 and 6 that specifies the type
    of example. Details are given below.
}
  \item{family}{Distribution family of the response variable.
Options are "binomial" and "poisson".
}
  \item{numRep}{The number of repeated measures within each group.
}
  \item{numGrp}{The number of groups.
}
\item{randomFactor}{Boolean flag that indicates whether or not
  the the model has a random intercept. If \code{randomFactor}
  is \code{FALSE} then the random intercept is omitted and
  \code{numRep} is automatically set to 1.
}
}
\details{
The example types are (as specified by the value of \code{eg}:\cr
   1. Example contains only one linear function.\cr
   2. Example contains two linear functions.\cr
   3. Example has only one smooth function.\cr
   4. Example contains one linear function and one smooth function.\cr
   5. Example contains only two smooth functions.\cr
   6. Example contains two linear functions and two smooth functions.
}
\value{A data frame containing the response and predictors as
  columns. If \code{randomFactor} is true then a column of
  containing group labels, and named \code{idnum}, is also
  included.
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
\seealso{
  \code{\link{gSlc}}
}
\examples{
dat1 <- gSlcSim(eg = 6, family = "binomial", randomFactor = TRUE)

dat2 <- gSlcSim(eg = 4, numGrp = 500, family = "poisson",
                randomFactor = FALSE)
}
%
\keyword{models}
\keyword{regression}
