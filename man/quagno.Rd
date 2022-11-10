\name{quagno}
\alias{quagno}
\title{Quantile Function of the Generalized Normal Distribution}
\description{
This function computes the quantiles of the Generalized Normal (Log-Normal3) distribution given parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) computed by \code{\link{pargno}}.  The quantile function has no explicit form. The parameters have the following interpretations: \eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter, and \eqn{\kappa} is a shape parameter.
}
\usage{
quagno(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{pargno}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
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
\seealso{\code{\link{cdfgno}}, \code{\link{pdfgno}}, \code{\link{lmomgno}}, \code{\link{pargno}}, \code{\link{qualn3}}
}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  quagno(0.5,pargno(lmr))
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Generalized Normal}
