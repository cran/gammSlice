\name{indonRespir}
\alias{indonRespir}
\usage{data(indonRespir)}
\title{Eespiratory infection in Indonesian children}
\description{Indonesian Children's Health Study of respiratory infections for a cohort of 275 Indonesian children. The data are longitudinal with each child having between 1 and 6 repeated measurements.
}
\format{
  A data frame with 1200 observations on the following 12 variables:
  \describe{
    \item{\code{idnum}}{child identification number.}
    \item{\code{respirInfec}}{indicator of presence of resipiratory infection.}
    \item{\code{age}}{age of the child in years.}
    \item{\code{vitAdefic}}{indicator of Vitamin A deficiency:\cr
                               1 = the child had Vitamin A deficiency,\cr
                               0 = the child did not have Vitamin A deficiency.}
    \item{\code{female}}{indicator of child being female:\cr
                               1 = the child is female,\cr
                               0 = the child is male.}
    \item{\code{height}}{height of the child in centimeters.}
    \item{\code{stunted}}{indicator of the child being "short for his/her age":\cr
                               1 = the child is "short for his/her age",\cr
                               0 = the child is not "short for his/her age"}
    \item{\code{visit2}}{indicator that the child had exactly 2 clinical visits:\cr
                               1 = the exact number of clinical visits was 2,\cr
                               0 = the exact number of clinical visits was not 2.}
    \item{\code{visit3}}{indicator that the child had exactly 3 clinical visits:\cr
                               1 = the exact number of clinical visits was 3,\cr
                               0 = the exact number of clinical visits was not 3.}
    \item{\code{visit4}}{indicator that the child had exactly 4 clinical visits:\cr
                               1 = the exact number of clinical visits was 4,\cr
                               0 = the exact number of clinical visits was not 4.}
    \item{\code{visit5}}{indicator that the child had exactly 5 clinical visits:\cr
                               1 = the exact number of clinical visits was 5,\cr
                               0 = the exact number of clinical visits was not 5.}
    \item{\code{visit6}}{indicator that the child had exactly 6 clinical visits:\cr
                               1 = the exact number of clinical visits was 6,\cr
                               0 = the exact number of clinical visits was not 6.}
  }
}
\source{Sommer, A. (1982). \emph{Nutritional Blindness}. New York: Oxford University Press.
}
\references{Diggle, P., Heagerty, P., Liang, K.-L. and Zeger, S. (2002). \emph{Analysis of Longitudinal Data (Second Edition).} Oxford: Oxford University Press.
}
\examples{
library(gammSlice) ; data(indonRespir)
plot(indonRespir$age,jitter(indonRespir$respirInfec))
}

