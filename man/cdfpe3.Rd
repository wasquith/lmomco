\name{cdfpe3}
\alias{cdfpe3}
\title{Cumulative Distribution Function of the Pearson Type III Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Pearson Type III distribution given parameters (\eqn{\mu}, \eqn{\sigma}, and \eqn{\gamma}) computed by \code{\link{parpe3}}. These parameters are equal to the product moments: mean, standard deviation, and skew (see \code{\link{pmoms}}). The cumulative distribution function is
\deqn{F(x) = \frac{G\left(\alpha,\frac{Y}{\beta}\right)}{\Gamma(\alpha)} \mbox{,}}
for \eqn{\gamma \ne 0} and where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x}, \eqn{G} is defined below and is related to the incomplete gamma function of \R (\code{pgamma()}), \eqn{\Gamma} is the complete gamma function, \eqn{\xi} is a location parameter, \eqn{\beta} is a scale parameter, \eqn{\alpha} is a shape parameter, and \eqn{Y = x - \xi} if \eqn{\gamma > 0} and \eqn{Y = \xi - x} if \eqn{\gamma < 0}. These three \dQuote{new} parameters are related to the product moments by
\deqn{\alpha = 4/\gamma^2 \mbox{,}}
\deqn{\beta  = \frac{1}{2}\sigma |\gamma| \mbox{,}}
\deqn{\xi = \mu - 2\sigma/\gamma \mbox{.}}
Lastly,  the function \eqn{G(\alpha,x)} is
\deqn{G(\alpha,x) = \int_0^x t^{(a-1)} \exp(-t)\, \mathrm{d}t \mbox{.}}

If \eqn{\gamma = 0}, the distribution is symmetrical and simply is the normal distribution with mean and standard deviation of \eqn{\mu} and \eqn{\sigma}, respectively. Internally, the \eqn{\gamma = 0} condition is implemented by \code{pnorm()}. If \eqn{\gamma > 0}, the distribution is right-tail heavy, and \eqn{F(x)} is the returned nonexceedance probability. On the other hand if \eqn{\gamma < 0}, the distribution is left-tail heavy and \eqn{1-F(x)} is the actual nonexceedance probability that is returned.
}
\usage{
cdfpe3(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parpe3}} or  \code{\link{vec2par}}.}
}
\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal
of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments:
Version 3, IBM Research Report RC20525, T.J. Watson Research Center,
Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{pdfpe3}}, \code{\link{quape3}}, \code{\link{lmompe3}}, \code{\link{parpe3}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  cdfpe3(50,parpe3(lmr))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Pearson Type III}
