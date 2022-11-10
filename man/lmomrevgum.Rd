\name{lmomrevgum}
\alias{lmomrevgum}
\title{L-moments of the Reverse Gumbel Distribution}
\description{
This function estimates the L-moments of the Reverse Gumbel distribution given the parameters (\eqn{\xi} and \eqn{\alpha}) from \code{\link{parrevgum}}. The first two type-B L-moments in terms of the parameters are
\deqn{\lambda^B_1 = \xi - (0.5722\dots) \alpha  - \alpha\lbrace\mathrm{Ei}(-\log(1-\zeta))\rbrace\mbox{and}}
\deqn{\lambda^B_2 = \alpha\lbrace\log(2) + \mathrm{Ei}(-2\log(1-\zeta)) - \mathrm{Ei}(-\log(1-\zeta))\rbrace\mbox{,}}
where \eqn{\zeta} is the right-tail censoring fraction of the sample or the nonexceedance probability of the right-tail censoring threshold, and \eqn{\mathrm{Ei}(x)} is the exponential integral defined as
\deqn{ \mathrm{Ei}(X) = \int_X^{\infty} x^{-1}\mathrm{exp}(-x)\mathrm{d}x \mbox{,}}
where \eqn{\mathrm{Ei}(-\log(1-\zeta)) \rightarrow 0} as \eqn{\zeta \rightarrow 1} and \eqn{\mathrm{Ei}(-\log(1-\zeta))} can not be evaluated as \eqn{\zeta \rightarrow 0}.
}
\usage{
lmomrevgum(para)
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
  \item{zeta}{Number of samples observed (noncensored) divided by the total number of samples.}
  \item{source}{An attribute identifying the computational
source of the L-moments: \dQuote{lmomrevgum}.}
}
\references{
Hosking, J.R.M., 1995, The use of L-moments in the analysis of censored data, in Recent Advances in Life-Testing and Reliability, edited by N. Balakrishnan, chapter 29, CRC Press, Boca Raton, Fla., pp. 546--560.
}
\author{W.H. Asquith}
\seealso{\code{\link{parrevgum}}, \code{\link{cdfrevgum}}, \code{\link{pdfrevgum}}, \code{\link{quarevgum}} }
\examples{
lmr <- lmoms(c(123,34,4,654,37,78))
rev.para <- lmom2par(lmr,type='revgum')
lmomrevgum(rev.para)
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Gumbel (reversed)}
