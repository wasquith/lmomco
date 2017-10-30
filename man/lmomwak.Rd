\name{lmomwak}
\alias{lmomwak}
\title{L-moments of the Wakeby Distribution}
\description{
This function estimates the L-moments of the Wakeby distribution given the parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, and \eqn{\delta}) from \code{\link{parwak}}. The L-moments in terms of the parameters are complicated and solved numerically.
}
\usage{
lmomwak(wakpara)
}
\arguments{
  \item{wakpara}{The parameters of the distribution.}
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
source of the L-moments: \dQuote{lmomwak}.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{\code{\link{parwak}}, \code{\link{cdfwak}}, \code{\link{pdfwak}}, \code{\link{quawak}} }
\examples{
lmr <- lmoms(c(123,34,4,654,37,78))
lmr
lmomwak(parwak(lmr))
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Wakeby}
