\name{lmomlap}
\alias{lmomlap}
\title{L-moments of the Laplace Distribution}
\description{
This function estimates the L-moments of the Laplace distribution given the parameters (\eqn{\xi} and \eqn{\alpha}) from \code{\link{parlap}}. The L-moments in terms of the parameters are 
\eqn{\lambda_1 = \xi}, 
\eqn{\lambda_2 = 3\alpha/4}, 
\eqn{\tau_3 = 0}, 
\eqn{\tau_4 = 17/22}, 
\eqn{\tau_5 = 0}, and 
\eqn{\tau_6 = 31/360}.

For \eqn{r} odd and \eqn{r \ge 3}, \eqn{\lambda_r = 0}, and for \eqn{r} even and \eqn{r \ge 4}, the L-moments using the hypergeometric function \eqn{{}_2F_1()} are
\deqn{\lambda_r = \frac{2\alpha}{r(r-1)}[1 - {}_2F_1(-r, r-1, 1, 1/2)]\mbox{,}}
where \eqn{{}_2F_1(a, b, c, z)} is defined as 
\deqn{{}_2F_1(a, b, c, z) = \sum_{n=0}^\infty \frac{(a)_n(b)_n}{(c)_n}\frac{z^n}{n!}\mbox{,}}
where \eqn{(x)_n} is the \emph{rising} Pochhammer symbol, which is defined by
\deqn{(x)_n =  1 \mbox{\ for\ } n = 0\mbox{, and}}
\deqn{(x)_n = x(x+1)\cdots(x+n-1) \mbox{\ for\ } n > 0\mbox{.}}
}
\usage{
lmomlap(para)
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
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{lmomlap}.}
}
\references{
Hosking, J.R.M., 1986, The theory of probability weighted moments: IBM Research Report RC12210, T.J. Watson Research Center, Yorktown Heights, New York.
}
\author{ W.H. Asquith}
\seealso{\code{\link{parlap}}, \code{\link{cdflap}}, \code{\link{pdflap}}, \code{\link{qualap}}  }
\examples{
lmr <- lmoms(c(123,34,4,654,37,78))
lmr
lmomlap(parlap(lmr))
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Laplace}
