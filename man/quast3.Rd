\name{quast3}
\alias{quast3}
\title{Quantile Function of the 3-Parameter Student t Distribution}
\description{
This function computes the quantiles of the 3-parameter Student t distribution given parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\nu}) computed by \code{\link{parst3}}. There is no explicit solution for the quantile function for nonexceedance probability \code{F} but built-in \R functions can be used. The implementation is \code{U} = \eqn{\xi} and \code{A} = \eqn{\alpha} for \eqn{1.001 \le \nu \le 10^5.5}, one can use \code{U + A*qt(F, N)} where \code{qt} is the 1-parameter Student t quantile function. The numerically accessible range of implementation here and consistency to the \eqn{\tau_4} and \eqn{\tau_6} is \eqn{10.001 \le \nu \le 10^5.5}. The limits for \eqn{\nu} stem from study of ability for theoretical integration of the quantile function to produce viable \eqn{\tau_4} and \eqn{\tau_6} (see \code{inst/doc/t4t6/studyST3.R}).
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
\author{W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{cdfst3}}, \code{\link{pdfst3}}, \code{\link{lmomst3}}, \code{\link{parst3}} }
\examples{
lmr <- lmoms(c(123, 34, 4, 654, 37, 78))
quast3(0.75, parst3(lmr))
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Student t (3-parameter)}
