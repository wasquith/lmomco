\name{quawak}
\alias{quawak}
\title{Quantile Function of the Wakeby Distribution}
\description{
This function computes the quantiles of the Wakeby distribution given
parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, and \eqn{\delta}) computed by \code{\link{parwak}}. The quantile function  is
\deqn{x(F) = \xi+\frac{\alpha}{\beta}(1-(1-F)^\beta)-
                 \frac{\gamma}{\delta}(1-(1-F))^{-\delta} \mbox{,}}
where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F}, \eqn{\xi} is a location parameter, \eqn{\alpha} and \eqn{\beta} are scale parameters, and \eqn{\gamma} and \eqn{\delta} are shape parameters. The five returned parameters from \code{parwak} in order are \eqn{\xi}, \eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, and \eqn{\delta}.
}
\usage{
quawak(f, wakpara, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{wakpara}{The parameters from \code{\link{parwak}} or \code{\link{vec2par}}.}
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
\seealso{\code{\link{cdfwak}}, \code{\link{pdfwak}}, \code{\link{lmomwak}}, \code{\link{parwak}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  quawak(0.5,parwak(lmr))
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Wakeby}
