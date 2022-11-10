\name{lmomray}
\alias{lmomray}
\title{L-moments of the Rayleigh Distribution}
\description{
This function estimates the L-moments of the Rayleigh distribution given the parameters (\eqn{\xi} and \eqn{\alpha}) from \code{\link{parray}}. The L-moments in terms of the parameters are
\deqn{\lambda_1 = \xi + \alpha\sqrt{\pi/2} \mbox{,}}
\deqn{\lambda_2 = \frac{1}{2} \alpha(\sqrt{2} - 1)\sqrt{\pi}\mbox{,}}
\deqn{\tau_3 = \frac{1 - 3/\sqrt{2} + 2/\sqrt{3}}{1 - 1/\sqrt{2}} = 0.1140 \mbox{, and}}
\deqn{\tau_4 = \frac{1 - 6/\sqrt{2} +  10/\sqrt{3} - 5\sqrt{4}}{1 - 1/\sqrt{2}} = 0.1054 \mbox{.}}
}
\usage{
lmomray(para)
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
source of the L-moments: \dQuote{lmomray}.}
}
\references{
Hosking, J.R.M., 1986, The theory of probability weighted moments: Research Report RC12210, IBM Research Division, Yorkton Heights, N.Y.
}
\author{W.H. Asquith}
\seealso{\code{\link{parray}}, \code{\link{cdfray}}, \code{\link{pdfray}}, \code{\link{quaray}} }
\examples{
lmr <- lmoms(c(123,34,4,654,37,78))
lmr
lmomray(parray(lmr))
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Rayleigh}
