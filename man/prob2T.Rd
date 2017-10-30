\name{prob2T}
\alias{prob2T}
\title{Convert a Vector of Annual Nonexceedance Probabilities to T-year Return Periods}
\description{
This function converts a vector of annual nonexceedance probabilities \eqn{F} to \eqn{T}-year return periods
\deqn{T = \frac{1}{1 - F}\mbox{,}}
where \eqn{0 \le F \le 1}.
}
\usage{
prob2T(f)
}
\arguments{
  \item{f}{A vector of annual nonexceedance probabilities.}
}
\value{
   A vector of \eqn{T}-year return periods.
}
\author{ W.H. Asquith}
\seealso{\code{\link{T2prob}}, \code{\link{nonexceeds}}, \code{\link{add.lmomco.axis}}, \code{\link{prob2grv}}, \code{\link{prob2lrv}}
}
\examples{
F <- nonexceeds()
T <- prob2T(F)
}
\keyword{utility (probability)}
