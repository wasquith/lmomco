\name{lmomlap}
\alias{lmomlap}
\title{L-moments of the Laplace Distribution}
\description{
This function estimates the L-moments of the Laplace distribution given the parameters (\eqn{\xi} and \eqn{\alpha}) from \code{\link{parlap}}. The L-moments in terms of the parameters are
\eqn{\lambda_1 = \xi},
\eqn{\lambda_2 = 3\alpha/4},
\eqn{\tau_3 = 0},
\eqn{\tau_4 = 17/72},
\eqn{\tau_5 = 0}, and
\eqn{\tau_6 = 31/360}.

For \eqn{r} odd and \eqn{r \ge 3}, \eqn{\lambda_r = 0}, and for \eqn{r} even and \eqn{r \ge 4}, the L-moments using the hypergeometric function \eqn{{}_2F_1()} (refer to the \pkg{hypergeo} package) are
\deqn{\lambda_r = \frac{2\alpha}{r(r-1)}\bigl[1 - {}_2F_1(-r, r-1, 1, 1/2)\bigr]\mbox{,}}
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
\details{
Hosking (1986) presents the \eqn{\lambda_2}, \eqn{\tau_4}, and \eqn{\tau_6} as integer ratios relative to whatever the \eqn{\alpha} parameter is. We can cross check those readily with the \pkg{hypergeo} and \pkg{numbers} packages:
\preformatted{
  L2  <- (1-hypergeo::hypergeo(-2, 2-1, 1, 1/2)) * (2/(2*(2-1))) # 0.75+0i = 3/4
  r   <- c(4, 6)
  L46 <- (1-hypergeo::hypergeo(-r, r-1, 1, 1/2)) * (2/(r*(r-1)))
  print(L46 / L2, digits=16) # Tau4 and Tau6
  # 0.2361111111111111+0i 0.0861111111111111+0i
  Tau <- as.numeric( L46 / L2 )
  numbers::contfrac(Tau[1], tol=.0001 ) # no better than 17 /  72
  numbers::contfrac(Tau[2], tol=.00001) # no better than 31 / 360
  # These results are consistent with Hosking (1986) and no better accuracy is
  # possible on "modern" computers because the 0***11111...... and check by
  # continued fraction providing also the integer ratio to a tolerance.
}
}
\references{
Hosking, J.R.M., 1986, The theory of probability weighted moments: IBM Research Report RC12210, T.J. Watson Research Center, Yorktown Heights, New York.
}
\author{W.H. Asquith}
\seealso{\code{\link{parlap}}, \code{\link{cdflap}}, \code{\link{pdflap}}, \code{\link{qualap}}  }
\examples{
lmr <- lmoms(c(123,34,4,654,37,78)); para <- parlap(lmr, minvar=FALSE)
print( str(          lmr  )  )
print( str( lmomlap( para ) ))

\dontrun{
  para <- vec2par(c(0, 4), type="lap")
  q <- rlmomco(3000, para); lmr <- lmoms(q)
  parlap(lmr, minvar=TRUE )$para
  parlap(lmr, minvar=FALSE)$para # }
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Laplace}
\keyword{Package: hypergeo}
\keyword{Package: numbers}
