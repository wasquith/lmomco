\name{lmomnor}
\alias{lmomnor}
\title{L-moments of the Normal Distribution}
\description{
This function estimates the L-moments of the Normal distribution given the parameters (\eqn{\mu} and \eqn{\sigma}) from \code{\link{parnor}}. The L-moments in terms of the parameters are
\eqn{\lambda_1 = \mu},
\eqn{\lambda_2 = \sigma / \sqrt{pi}},
\eqn{\tau_3 = 0},
\eqn{\tau_4 = 0.122602}, and
\eqn{\tau_5 = 0}.
}
\usage{
lmomnor(para)
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
source of the L-moments: \dQuote{lmomnor}.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{parnor}}, \code{\link{cdfnor}}, \code{\link{pdfnor}}, \code{\link{quanor}} }
\examples{
lmr <- lmoms(c(123, 34, 4, 654, 37, 78))
lmr
lmomnor(parnor(lmr))
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Normal}
