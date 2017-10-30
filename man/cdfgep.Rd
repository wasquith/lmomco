\name{cdfgep}
\alias{cdfgep}
\title{Cumulative Distribution Function of the Generalized Exponential Poisson Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Generalized Exponential Poisson distribution given parameters (\eqn{\beta}, \eqn{\kappa}, and \eqn{h}) computed by \code{\link{pargep}}. The cumulative distribution function is
\deqn{F(x) =  \left(\frac{1 - \exp[-h + h\exp(-\eta x)]}{1 - \exp(-h)}\right)^\kappa\mbox{,}}
where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x > 0}, \eqn{\eta = 1/\beta}, \eqn{\beta > 0} is a scale parameter, \eqn{\kappa > 0} is a shape parameter, and \eqn{h > 0} is another shape parameter.
}
\usage{
cdfgep(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{pargep}} or  \code{\link{vec2par}}.}
}
\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Barreto-Souza, W., and Cribari-Neto, F., 2009, A generalization of the exponential-Poisson distribution: Statistics and Probability, 79, pp. 2493--2500.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdfgep}}, \code{\link{quagep}}, \code{\link{lmomgep}}, \code{\link{pargep}}}
\examples{
gep <- list(para=c(2, 1.5, 3), type="gep")
cdfgep(0.48,gep)
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Generalized Exponential Poisson}
