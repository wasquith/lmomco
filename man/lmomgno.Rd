\name{lmomgno}
\alias{lmomgno}
\title{L-moments of the Generalized Normal Distribution}
\description{
This function estimates the L-moments of the Generalized Normal (Log-Normal3) distribution given the parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) from \code{\link{pargno}}. The L-moments in terms of the parameters are
\deqn{\lambda_1 = \xi + \frac{\alpha}{\kappa}(1-\mathrm{exp}(\kappa^2/2) \mbox{, and}}
\deqn{\lambda_2 = \frac{\alpha}{\kappa}(\mathrm{exp}(\kappa^2/2)(1-2\Phi(-\kappa/\sqrt{2})) \mbox{,}}
where \eqn{\Phi} is the cumulative distribution of the Standard Normal distribution. There are no simple expressions for \eqn{\tau_3}, \eqn{\tau_4}, and \eqn{\tau_5}. Logarthmic transformation of the data prior to fitting of the Generalized Normal distribution is not required. The distribution is algorithmically the same with subtle parameter modifications as the Log-Normal3 distribution (see \code{\link{lmomln3}}, \code{\link{parln3}}). If desired for user-level control of the lower bounds of a Log-Normal-like distribution is required, then see \code{\link{parln3}}.
}
\usage{
lmomgno(para)
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
source of the L-moments: \dQuote{lmomgno}.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{pargno}}, \code{\link{cdfgno}}, \code{\link{pdfgno}}, \code{\link{quagno}}, \code{\link{lmomln3}}
}
\examples{
lmr <- lmoms(c(123,34,4,654,37,78))
lmr
lmomgno(pargno(lmr))
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Generalized Normal}
