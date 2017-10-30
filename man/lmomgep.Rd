\name{lmomgep}
\alias{lmomgep}
\title{L-moments of the Generalized Exponential Poisson Distribution}
\description{
This function estimates the L-moments of the Generalized Exponential Poisson (GEP) distribution given the parameters (\eqn{\beta}, \eqn{\kappa}, and \eqn{h}) from \code{\link{pargep}}. The L-moments in terms of the parameters are best expressed in terms of the expectations of order statistic maxima \eqn{\mathrm{E}[X_{n:n}]} for the distribution. The fundamental relation is
\deqn{\lambda_r = \sum_{k=1}^r (-1)^{r-k}k^{-1}{r-1 \choose k-1}{r+k-2 \choose k-1}\mathrm{E}[X_{k:k}]\mbox{.}}
The L-moments do not seem to have been studied for the GEP. The challenge is the solution to \eqn{\mathrm{E}[X_{n:n}]} through an expression by Barreto-Souza and Cribari-Neto (2009) that is
\deqn{\mathrm{E}[X_{n:n}] = \frac{\beta\,h\,\Gamma(\kappa+1)\,\Gamma(n\kappa + 1)}{n\,\Gamma(n)\,(1 - \exp(-h))^{n\kappa}}\sum_{j=0}^{\infty} \frac{(-1)^j\exp(-h(j+1))}{\Gamma(n\kappa - j)\,\Gamma(j+1)}\;F^{12}_{22}(h(j+1))\mbox{,}}
where \eqn{F^{12}_{22}(h(j+1))} is the Barnes Extended Hypergeometric function with arguments reflecting those needed for the GEP (see comments under \code{\link{BEhypergeo}}).
}
\usage{
lmomgep(para, byqua=TRUE)
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{byqua}{A logical triggering the \code{\link{theoLmoms.max.ostat}} instead of using the mathematics of Barreto-Souza and Cribari-Neto (2009) (see Details).}
}
\details{
The mathematics (not of L-moments but \eqn{\mathrm{E}[X_{n:n}]}) shown by Barreto-Souza and Cribari-Neto (2009) are correct but are apparently subject to considerable numerical issues even with substantial use of logarithms and exponentiation in favor of multiplication and division in the above formula for \eqn{\mathrm{E}[X_{n:n}]}. Testing indicates that numerical performance is better if the non-\eqn{j}-dependent terms in the infinite sum remain \emph{inside} it.  Testing also indicates that the edges of performance can be readily hit with large \eqn{\kappa} and less so with large \eqn{h}. It actually seems superior to not use the above equation for L-moment computation based on \eqn{\mathrm{E}[X_{n:n}]} but instead rely on expectations of maxima order statistics (\code{\link{expect.max.ostat}}) from numerical integration of the quantile function (\code{\link{quagep}}) as is implementated in \code{\link{theoLmoms.max.ostat}}. This is the reason that the \code{byqua} argument is available and set to the shown default. Because the GEP is experimental, this function provides two approaches for \eqn{\lambda_r} computation for research purposes.
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
source of the L-moments: \dQuote{lmomgep}.}
}
\references{
Barreto-Souza, W., and Cribari-Neto, F., 2009, A generalization of the exponential-Poisson distribution: Statistics and Probability, 79, pp. 2493--2500.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pargep}}, \code{\link{cdfgep}}, \code{\link{pdfgep}}, \code{\link{quagep}} }
\examples{
\dontrun{
gep <- vec2par(c(2, 1.5, 3), type="gep")
lmrA <- lmomgep(gep, byqua=TRUE);   print(lmrA)
lmrB <- lmomgep(gep, byqua=FALSE);  print(lmrB)

# Because the L-moments of the Generalized Exponential Poisson are computed
# strictly from the expectations of the order statistic extrema, lets us evaluate
# by theoretical integration of the quantile function and simulation:
set.seed(10); gep <- vec2par(c(2, 1.5, 3), type="gep")
lmr  <- lmomgep(gep, byqua=FALSE)
E33a <- (lmr$lambdas[3] + 3*lmr$lambdas[2] + 2*lmr$lambdas[1])/2  # 2.130797
E33b <- expect.max.ostat(3, para=gep, qua=quagep)                 # 2.137250
E33c <- mean(replicate(20000, max(quagep(runif(3), gep))))        # 2.140226
# See how the E[X_{3:3}] by the formula shown in this documentation results in
# a value that is about 0.007 too small. Now this might now seem large but it
# is a difference.  Try gep <- list(para=c(2, 1.5, 13), type="gep") or
#  gep <- list(para=c(2, .08, 21), type="gep"), which fails on byqua=TRUE}
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Generalized Exponential Poisson}
