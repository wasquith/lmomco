\name{quagum}
\alias{quagum}
\title{Quantile Function of the Gumbel Distribution}
\description{
This function computes the quantiles of the Gumbel distribution given parameters (\eqn{\xi} and \eqn{\alpha})  computed by \code{\link{pargum}}. The quantile function  is
\deqn{x(F) = \xi - \alpha\log(-\log(F)) \mbox{,}}
where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F}, \eqn{\xi} is a location parameter, and \eqn{\alpha} is a scale parameter.
}
\usage{
quagum(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{pargum}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal
of the Royal Statistical Society, Series B, v. 52, p. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments:
Version 3, IBM Research Report RC20525, T.J. Watson Research Center,
Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfgum}}, \code{\link{pdfgum}}, \code{\link{lmomgum}}, \code{\link{pargum}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  quagum(0.5,pargum(lmr))
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Gumbel}
