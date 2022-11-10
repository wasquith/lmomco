\name{lmomlmrq}
\alias{lmomlmrq}
\title{L-moments of the Linear Mean Residual Quantile Function Distribution}
\description{
This function estimates the L-moments of the Linear Mean Residual Quantile Function distribution given the parameters (\eqn{\mu} and \eqn{\alpha}) from \code{\link{parlmrq}}. The first six L-moments in terms of the parameters are
\deqn{\lambda_1 = \mu \mbox{,}}
\deqn{\lambda_2 = (\alpha + 3\mu)/6 \mbox{,}}
\deqn{\lambda_3 = 0 \mbox{,}}
\deqn{\lambda_4 = (\alpha + \mu)/12 \mbox{,}}
\deqn{\lambda_5 = (\alpha + \mu)/20 \mbox{, and}}
\deqn{\lambda_6 = (\alpha + \mu)/30 \mbox{.}}

Because \eqn{\alpha + \mu > 0}, then \eqn{\tau_3 > 0}, so the distribution is positively skewed.  The coefficient of L-variation is in the interval \eqn{(1/3, 2/3)}.
}
\usage{
lmomlmrq(para)
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
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{lmomlmrq}.}
}
\references{
Midhu, N.N., Sankaran, P.G., and Nair, N.U., 2013, A class of distributions with linear mean residual quantile function and it's generalizations: Statistical Methodology, v. 15, pp. 1--24.
}
\author{W.H. Asquith}
\seealso{\code{\link{parlmrq}}, \code{\link{cdflmrq}}, \code{\link{pdflmrq}}, \code{\link{qualmrq}}
}
\examples{
lmr <- lmoms(c(3, 0.05, 1.6, 1.37, 0.57, 0.36, 2.2))
lmr
lmomlmrq(parlmrq(lmr))
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Linear Mean Residual Quantile Function}
