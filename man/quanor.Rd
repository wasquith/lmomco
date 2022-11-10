\name{quanor}
\alias{quanor}
\title{Quantile Function of the Normal Distribution}
\description{
This function computes the quantiles of the Normal distribution given parameters (\eqn{\mu} and \eqn{\sigma})  computed by \code{\link{parnor}}. The quantile function  has no explicit form (see \code{cdfnor} and \code{qnorm}). The parameters have the following interpretations: \eqn{\mu} is the arithmetic mean and \eqn{\sigma} is the standard deviation. The \R function \code{qnorm} is used.
}
\usage{
quanor(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parnor}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the  quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfnor}}, \code{\link{pdfnor}}, \code{\link{lmomnor}}, \code{\link{parnor}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  quanor(0.5,parnor(lmr))
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Normal}
