\name{T2prob}
\alias{T2prob}
\title{Convert a Vector of T-year Return Periods to Annual Nonexceedance Probabilities}
\description{
This function converts a vector of \eqn{T}-year return periods to annual nonexceedance probabilities \eqn{F}
\deqn{F = 1 - \frac{1}{T}\mbox{,}}
where \eqn{0 \le F \le 1}.
}
\usage{
T2prob(T)
}
\arguments{
  \item{T}{A vector of \eqn{T}-year return periods.}
}
\value{
   A vector of annual nonexceedance probabilities.
}
\author{ W.H. Asquith}
\seealso{\code{\link{prob2T}}, \code{\link{nonexceeds}}, \code{\link{add.lmomco.axis}}}
\examples{
T <- c(1, 2, 5, 10, 25, 50, 100, 250, 500)
F <- T2prob(T)
}
\keyword{utility (probability)}
