\name{lmomexp}
\alias{lmomexp}
\title{L-moments of the Exponential Distribution}
\description{
This function estimates the L-moments of the Exponential distribution given the parameters (\eqn{\xi} and \eqn{\alpha}) from \code{\link{parexp}}. The L-moments in terms of the parameters are \eqn{\lambda_1 = \xi + \alpha}, \eqn{\lambda_2 = \alpha/2}, \eqn{\tau_3 = 1/3}, \eqn{\tau_4 = 1/6}, and \eqn{\tau_5 = 1/10}.
}
\usage{
lmomexp(para)
}
\arguments{
  \item{para}{The parameters of the distribution.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{0}.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL}.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL}.}
  \item{source}{An attribute identifying the computational
source of the L-moments: \dQuote{lmomexp}.}
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{\code{\link{parexp}}, \code{\link{cdfexp}}, \code{\link{pdfexp}}, \code{\link{quaexp}} }
\examples{
lmr <- lmoms(c(123,34,4,654,37,78))
lmr
lmomexp(parexp(lmr))
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Exponential}
