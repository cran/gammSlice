\name{toenail}
\alias{toenail}
\usage{data(toenail)}
\title{Toenail infection clinical trial}
\description{Data from a clinical trial in which two anti-fungal treatments for toenail infection are compared.
}
\format{
  A data frame with 1908 observations on the following 5 variables:
  \describe{
    \item{\code{idnum}}{patient identification number.}
    \item{\code{onycholysis}}{indicator concerning the severity of onycholysis:\cr
                               1 = moderate or severe onycholysis,\cr
                               0 = no or mild onycholysis.}
    \item{\code{terb}}{indicator of whether the treatement was terbinafine:\cr
                               1 = treatment was terbinafine,\cr
                               0 = treatment was itraconazole.}
    \item{\code{months}}{time in months since the the start of the trial when clinical visit took place.}
    \item{\code{visit}}{visit number.}
  }
}
\references{De Backer, M., De Vroey, C., Lesaffre, E., Scheys, I. and De Keyser, P. (1998). Twelve weeks of continuous oral therapy for toenail onychomycosis caused by dermatophytes: a double-blind comparative trial of terbinafine 250 mg/day versus itraconazole 200 mg/day. \emph{Journal of the American Academy of Dermatology}, \bold{38}, S57-S63.
}
\examples{
library(gammSlice) ; data(toenail)
plot(jitter(toenail$terb),jitter(toenail$onycholysis))
}

