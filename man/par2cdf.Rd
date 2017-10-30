\name{par2cdf}
\alias{par2cdf}
\title{Cumulative Distribution Function of the Distributions}
\description{
This function acts as a front end or dispatcher to the distribution-specific cumulative distribution functions.
}
\usage{
par2cdf(x, para, ...)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
  \item{...}{The additional arguments are passed to the cumulative distribution function such as \code{paracheck=FALSE} for the Generalized Lambda distribution (\code{\link{cdfgld}}).}
}
\value{
  Nonexceedance probability (\eqn{0 \le F \le 1}) for \code{x}.
}
\author{ W.H. Asquith}
\seealso{\code{\link{par2pdf}}, \code{\link{par2qua}} }
\examples{
lmr       <- lmoms(rnorm(20))
para      <- parnor(lmr)
nonexceed <- par2cdf(0,para)
}
\keyword{cumulative distribution function}
