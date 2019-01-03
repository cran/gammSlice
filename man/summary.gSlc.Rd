\name{summary.gSlc}
\alias{summary.gSlc}
\usage{
\method{summary}{gSlc}(object,colour=TRUE,paletteNumber=1,...)}
\title{Summary of the generalized additive mixed model fit produced by \code{gSlc}}
\description{A graphical table showing, for key model parameters, the Markov chain Monte Carlo samples, diagnostic plots and numerical summaries.
}
\arguments{
   \item{object}{A \code{gSlc()} fit object.}
   \item{colour}{Boolean flag:\cr
                TRUE = produce a colour graphical table,\cr
                FALSE = produce a black and white graphical table .}
   \item{paletteNumber}{If \code{colour = TRUE} then there there are two possible colour palettes. These are determined by whether \code{paletteNumber} is set to 1 or 2. The default is \code{paletteNumber=1}.} 
   \item{...}{place-holder for additional arguments.}
}
\details{
  The columns of the graphical table are:\cr
  
   1. parameter name,\cr
   2. trace plot of the Markov chain Monte Carlo sample,\cr
   3. plot of Markov chain Monte Carlo sample against its lag 1 sample,\cr
   4. sample autocorrelation function,\cr
   5. kernel density estimate of the posterior density function,\cr
   6. posterior mean and 95\% credible interval.
 }
\references{
Pham, T.H. and Wand, M.P. (2018). Generalized additive mixed model analysis via \code{gammSlice}. \emph{Australian and New Zealand Journal of Statistics}, \bold{60}, 279-300. 
}
\author{
Tung Pham \email{tungstats@gmail.com} and Matt Wand \email{matt.wand@uts.edu.au}. 
}
\seealso{
\code{gSlc}, \code{plot.gSlc}
}
\examples{
library(gammSlice)
set.seed(39402) ; m <- 100 ; n <- 2
beta0True <- 0.5 ; betaxTrue <- 1.7 ; sigsqTrue <- 0.8
idnum <- rep(1:m,each=n) ; x <- runif(m*n)
U <- rep(rnorm(m,0,sqrt(sigsqTrue)),each=n)
mu <- 1/(1+exp(-(beta0True+betaxTrue*x+U)))
y <- rbinom((m*n),1,mu)
fit1 <- gSlc(y ~ x,random = list(idnum = ~1),family = "binomial",
             control = gSlc.control(nBurn=150,nKept=100,nThin=1))
summary(fit1)
summary(fit1,paletteNumber = 2)
summary(fit1,colour = FALSE)

\dontrun{
# Re-fit with higher Markov chain Monte Carlo sample:

fit2 <- gSlc(y ~ x,random = list(idnum = ~1),family = "binomial")
summary(fit2)
summary(fit2,paletteNumber = 2)
summary(fit2,colour = FALSE)
}
}



