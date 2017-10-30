\name{theoLmoms}
\alias{theoLmoms}
\title{ The Theoretical L-moments and L-moment Ratios using Integration of the Quantile Function}
\description{
Compute the theoretrical L-moments for a vector. A theoretrical L-moment in integral form is
\deqn{ \lambda_r = \frac{1}{r}
                               \sum^{r-1}_{k=0}{(-1)^k {r-1 \choose k}
        \frac{r!\:I_r}{(r-k-1)!\,k!}
       } \mbox{,}}
in which
\deqn{ I_r = \int^1_0 x(F) \times F^{r-k-1}(1-F)^{k}\,\mathrm{d}F \mbox{,}}
where \eqn{x(F)} is the quantile function of the random variable \eqn{X} for nonexceedance probability \eqn{F}, and \eqn{r} represents the order of the L-moments. This function actually dispatches to \code{\link{theoTLmoms}} with \code{trim=0} argument.
}
\usage{
theoLmoms(para, nmom=5, verbose=FALSE, minF=0, maxF=1)
}
\arguments{
  \item{para}{A distribution parameter object such as from \code{\link{vec2par}}.}
  \item{nmom}{The number of moments to compute. Default is 5.}
  \item{verbose}{Toggle verbose output. Because the \R function \code{integrate} is used to perform the numerical integration, it might be useful to see selected messages regarding the numerical integration.}
  \item{minF}{The end point of nonexceedance probability in which to perform the integration. Try setting to non-zero (but very small)  if the integral is divergent.}
  \item{maxF}{The end point of nonexceedance probability in which to perform the integration. Try setting to non-unity (but still very close [perhaps \code{1 - minF}]) if the integral is divergent.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the TL-moments. First element is \eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is \eqn{\tau_2}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which will equal zero (the ordinary L-moments).}
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{theoTLmoms}.}
}
\note{
 The actual function used is \code{theoTLmoms(para,nmom=nmom,trim=0,verbose=verbose)}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.
}
\author{ W.H. Asquith}
\seealso{\code{\link{theoTLmoms}}}
\examples{
para <- vec2par(c(0,1),type='nor') # standard normal
TL00 <- theoLmoms(para) # compute ordinary L-moments
}
\keyword{L-moment (theoretical)}

