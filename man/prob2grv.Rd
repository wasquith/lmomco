\name{prob2grv}
\alias{prob2grv}
\title{Convert a Vector of Annual Nonexceedance Probabilities to Gumbel Reduced Variates}
\description{
This function converts a vector of annual nonexceedance probabilities \eqn{F} to Gumbel reduced variates (GRV, \eqn{grv}; Hosking and Wallis [1997, p. 92]) 
\deqn{grv = -\log(-\log(F))\mbox{,}}
where \eqn{0 \le F \le 1}. The Gumbel distribution (\code{\link{quagum}}), which is a special case of the Generalized Extreme Value (\code{\link{quagev}}), will plot as a straightline when the horizontal axis is GRV transformed.
}
\usage{
prob2grv(f)
}
\arguments{
  \item{f}{A vector of annual nonexceedance probabilities.}
}
\value{
   A vector of Gumbel reduced variates.
}
\author{W.H. Asquith}
\references{
Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\seealso{\code{\link{grv2prob}}, \code{\link{prob2T}}
}
\examples{
F <- nonexceeds()
grv <- prob2grv(F)
}
\keyword{utility (probability)}
\concept{Gumbel reduced variate}
\concept{Gumbel Reduced Variate}
