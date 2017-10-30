\name{lmompe3}
\alias{lmompe3}
\title{L-moments of the Pearson Type III Distribution}
\description{
This function estimates the L-moments of the Pearson Type III distribution given the parameters (\eqn{\mu}, \eqn{\sigma}, and \eqn{\gamma}) from \code{\link{parpe3}} as the product moments: mean, standard deviation, and skew. The first three L-moments in terms of these parameters are complex and numerical methods are required. For simplier expression of the distribution functions (\code{\link{cdfpe3}}, \code{\link{pdfpe3}}, and \code{\link{quape3}}) the \dQuote{moment parameters} are expressed differently.

The Pearson Type III distribution is of considerable theoretical interest because the parameters, which are estimated via the L-moments, are in fact the product moments. Although, these values fitted by the method of L-moments will not be numerically equal to the sample product moments.  Further details are provided in the Examples section of the \code{\link{pmoms}} function documentation.
}
\usage{
lmompe3(para)
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
source of the L-moments: \dQuote{lmompe3}.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{\code{\link{parpe3}}, \code{\link{cdfpe3}}, \code{\link{pdfpe3}}, \code{\link{quape3}} }
\examples{
lmr <- lmoms(c(123,34,4,654,37,78))
lmr
lmompe3(parpe3(lmr))
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Pearson Type III}
