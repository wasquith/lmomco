\name{f2fpds}
\alias{f2fpds}
\title{Conversion of Annual Nonexceedance Probability to Partial Duration Nonexceedance Probability}
\description{
This function takes an annual exceedance probability and converts it to a \dQuote{partial-duration series} (a term in Hydrology) nonexceedance probability through a simple assumption that the Poisson distribution is appropriate for arrive modeling. The relation between the cumulative distribution function \eqn{G(x)} for the partial-duration series is related to the cumulative distribution function \eqn{F(x)} of the annual series (data on an annual basis and quite common in Hydrology) by

\deqn{G(x) = [\log(F(x)) + \eta]/\eta\mathrm{.}}

The core assumption is that successive events in the partial-duration series can be considered as \emph{independent}. The \eqn{\eta} term is the arrival rate of the events. For example, suppose that 21 events have occurred in 15 years, then \eqn{\eta = 21/15 = 1.4} events per year.

A comprehensive demonstration is shown in the example for \code{\link{fpds2f}}. That function performs the opposite conversion. Lastly, the cross reference to \code{\link{x2xlo}} is made because the example contained therein provides another demonstration of partial-duration and annual series frequency analysis.
}
\usage{
f2fpds(f, rate=NA)
}
\arguments{
  \item{f}{A vector of annual nonexceedance probabilities.}
  \item{rate}{The number of events per year.}
}
\value{
  A vector of converted nonexceedance probabilities.
}
\references{
Stedinger, J.R., Vogel, R.M., Foufoula-Georgiou, E., 1993, Frequency analysis of extreme events: \emph{in} Handbook of Hydrology, ed. Maidment, D.R., McGraw-Hill, Section 18.6 Partial duration series, mixtures, and censored data, pp. 18.37--18.39.
}
\author{W.H. Asquith}
\seealso{\code{\link{fpds2f}}, \code{\link{x2xlo}}, \code{\link{f2flo}}, \code{\link{flo2f}}  }
\examples{
# See examples for fpds2f().
}
\keyword{utility (probability)}
