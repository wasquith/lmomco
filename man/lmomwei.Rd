\name{lmomwei}
\alias{lmomwei}
\title{L-moments of the Weibull Distribution}
\description{
This function estimates the L-moments of the Weibull distribution given the parameters
(\eqn{\zeta}, \eqn{\beta}, and \eqn{\delta}) from \code{\link{parwei}}. The Weibull distribution is a reverse Generalized Extreme Value distribution.  As result, the Generalized Extreme Value algorithms (\code{\link{lmomgev}}) are used for computation of the L-moments of the Weibull in this package (see \code{\link{parwei}}).
}
\usage{
lmomwei(para)
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
source of the L-moments: \dQuote{lmomwei}.}
}
\references{
Hosking, J.R.M. and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{\code{\link{parwei}}, \code{\link{cdfwei}}, \code{\link{pdfwei}}, \code{\link{quawei}} }
\examples{
lmr <- lmoms(c(123,34,4,654,37,78))
lmr
lmomwei(parwei(lmr))
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Weibull}
