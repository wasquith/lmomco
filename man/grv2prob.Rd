\name{grv2prob}
\alias{grv2prob}
\title{Convert a Vector of Gumbel Reduced Variates to Annual Nonexceedance Probabilities}
\description{
This function converts a vector of Gumbel reduced variates (\eqn{grv}) to annual nonexceedance probabilities \eqn{F}
\deqn{F = \exp(-\exp(-grv))\mbox{,}}
where \eqn{0 \le F \le 1}.
}
\usage{
grv2prob(grv)
}
\arguments{
  \item{grv}{A vector of Gumbel reduced variates.}
}
\value{
   A vector of annual nonexceedance probabilities.
}
\author{W.H. Asquith}
\seealso{\code{\link{prob2grv}}, \code{\link{prob2T}}
}
\examples{
T <- c(1, 2, 5, 10, 25, 50, 100, 250, 500); grv <- prob2grv(T2prob(T))
F <- grv2prob(grv)
}
\keyword{utility (probability)}
\concept{Gumbel reduced variate}
\concept{Gumbel Reduced Variate}
