\name{pp.f}
\alias{pp.f}
\title{ Quantile Function of the Ranks of Plotting Positions }
\description{
There are two major forms (outside of the general plotting-position formula \code{\link{pp}}) for estimation of the \eqn{p_r}th probability of the \eqn{r}th order statistic for a sample of size \eqn{n}: the mean is \eqn{pp'_r = r/(n+1)} (Weibull plotting position) and the Beta quantile function is \eqn{pp_r(F) = IIB(F, r, n+1-r)}, where \eqn{F} represents the nonexceedance probability of the plotting position.  \eqn{IIB} is the \dQuote{inverse of the incomplete beta function} or the quantile function of the Beta distribution as provided in \R by \code{qbeta(f, a, b)}.  If \eqn{F=0.5}, then the median is returned but that is conveniently implemented in \code{\link{pp.median}}.   Readers might consult Gilchrist (2011, chapter 12) and Karian and Dudewicz (2011, p. 510).
}
\usage{
pp.f(f, x)
}
\arguments{
  \item{f}{A nonexceedance probability.}
  \item{x}{A vector of data. The ranks and the length of the vector are computed within the function.}
}
\value{
  An \R \code{vector} is returned.
}
\note{
The function uses the \R function \code{rank}, which has specific settings to handle tied data. For implementation here, the \code{ties.method="first"} method to \code{rank} is used.
}
\references{
Gilchrist, W.G., 2000, Statistical modelling with quantile functions: Chapman and Hall/CRC, Boca Raton.

Karian, Z.A., and Dudewicz, E.J., 2011, Handbook of fitting statistical distributions with R: Boca Raton, FL, CRC Press.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{pp}}, \code{\link{pp.median}} }
\examples{
X <- sort(rexp(10))
PPlo <- pp.f(0.25, X)
PPhi <- pp.f(0.75, X)
plot(c(PPlo,NA,PPhi), c(X,NA,X))
points(pp(X), X) # Weibull i/(n+1)
}
\keyword{plotting position}
\keyword{rankit}
