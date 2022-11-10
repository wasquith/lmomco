\name{BEhypergeo}
\alias{BEhypergeo}
\alias{BarnesExtendedHypergeometric}
\title{Barnes Extended Hypergeometric Function}
\description{
This function computes the Barnes Extended Hypergeometric function, which in \pkg{lmomco} is useful in applications involving expectations of order statistics for the Generalized Exponential Poisson (GEP) distribution (see \code{\link{lmomgep}}). The function is
\deqn{
F_{p,q}(\bm{\mathrm{n}},\bm{\mathrm{d}}; \lambda) = \sum_{k=0}^\infty \frac{\lambda^k}{\Gamma(k+1)}\frac{\Pi_{i=1}^{p} \Gamma(n_i + k)\Gamma^{-1}{(n_i)}}{\Pi_{i=1}^{q} \Gamma(d_i + k)\Gamma^{-1}{(d_i)}}\mbox{,} 
}
where \eqn{\bm{\mathrm{n}} = [n_1, n_2, \ldots, n_p]} for \eqn{p} operands and \eqn{\bm{\mathrm{d}} = [d_1, d_2, \ldots, d_q]} for \eqn{q} operands, and \eqn{\lambda > 0} is a parameter.
}
\usage{
BEhypergeo(p,q, N,D, lambda, eps=1E-12, maxit=500)
}
\arguments{
  \item{p}{An integer value.}
  \item{q}{An integer value.}
  \item{N}{A scalar or vector associated with the \eqn{p} summation (see Details).}
  \item{D}{A scalar or vector associated with the \eqn{q} summation (see Details).}
  \item{lambda}{A real value \eqn{\lambda > 0}.}
  \item{eps}{The relative convergence error on which to break an infinite loop.}
  \item{maxit}{The maximum number of interations before a mandatory break on the loop, and a warning will be issued.}
}
\details{
 For the GEP both \eqn{\bm{\mathrm{n}}} and \eqn{\bm{\mathrm{d}}} are vectors of the same value, such as \eqn{\bm{\mathrm{n}} = [1, \ldots, 1]} and \eqn{\bm{\mathrm{d}} = [2, \ldots, 2]}. This implementation is built around this need by the GEP and if the length of either vector is not equal to the operand then the first value of the vector is repeated the operand times. For example for \eqn{\bm{\mathrm{n}}}, if \code{n = 1}, then \code{n = rep(n[1], length(p))} and so on for \eqn{\bm{\mathrm{d}}}. Given that \code{n} and \code{d} are vectorized for the GEP, then a shorthand is used for the GEP mathematics shown herein:
\deqn{F^{12}_{22}(h(j+1)) \equiv F_{2,2}([1,\ldots,1], [2,\ldots,2]; h(j+1))\mbox{,}} 
for the \eqn{h} parameter of the distribution.

Lastly, for \pkg{lmomco} and the GEP the arguments only involve \eqn{p = q = 2} and \eqn{N = 1}, \eqn{D = 2}, so the function is uniquely a function of the \eqn{h} parameter of the distribution:
\preformatted{
  H <- 10^seq(-10,10, by=0.01)
  F22 <- sapply(1:length(H), function(i) BEhypergeo(2,2,1,1, H[i])$value
  plot(log10(H), log10(F22), type="l")
}
For this example, the solution increasingly wobbles towards large \eqn{h}, which is further explored by
\preformatted{
  plot(log10(H[1:(length(H)-1)]), diff(log10(F22)), type="l", xlim=c(0,7))
  plot(log10(H[H > 75 & H < 140]), c(NA,diff(log10(F22[H > 75 & H < 140]))),
       type="b"); lines(c(2.11,2.11), c(0,10))
}
It can be provisionally concluded that the solution to \eqn{F^{12}_{22}(\cdot)} begins to be suddenly questionable because of numerical difficulties beyond \eqn{\log(h) = 2.11}. Therefore, it is given that \eqn{h < 128} might be an operational numerical upper limit.
}
\value{
An \R \code{list} is returned.

\item{value}{The value for the function.}
\item{its}{The number of iterations \eqn{j}.}
\item{error}{The error of convergence.} 
}
\author{W.H. Asquith}
\references{
Kus, C., 2007, A new lifetime distribution: Computational Statistics and Data Analysis, v. 51, pp. 4497--4509.
}
\seealso{\code{\link{lmomgep}}}
\examples{
BEhypergeo(2,2,1,2,1.5)
}
\keyword{special function}
