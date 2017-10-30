\name{quast3}
\alias{quast3}
\title{Quantile Function of the 3-Parameter Student t Distribution}
\description{
This function computes the quantiles of the 3-parameter Student t distribution given parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\nu})  computed by \code{\link{parst3}}. There is no explicit solution for the quantile function for nonexceedance probability \code{F} but built-in \R functions can be used. For \eqn{\nu \ge 1000}, one can use \code{qnorm(F, mean=U,} \code{sd=A)} and for \code{U} = \eqn{\xi} and \code{A} = \eqn{\alpha} for \eqn{1.000001 \le \nu \le 1000}, one can use \code{U +} \code{A*qt(F, N)} for \code{N} = \eqn{\nu} and where the \R function \code{qnorm} is the Normal distribution and \R function \code{qt} is the 1-parameter Student t distribution.
}
\usage{
quast3(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parst3}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical on whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\author{ W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{cdfst3}}, \code{\link{pdfst3}}, \code{\link{lmomst3}}, \code{\link{parst3}} }
\examples{
lmr <- lmoms(c(123,34,4,654,37,78))
quast3(0.75,parst3(lmr))
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Student t (3-parameter)}
