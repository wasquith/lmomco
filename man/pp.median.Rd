\name{pp.median}
\alias{pp.median}
\title{ Quantile Function of the Ranks of Plotting Positions }
\description{
The median of a plotting position. The median is \eqn{pp^\star_r = IIB(0.5, r, n+1-r)}.  \eqn{IIB} is the \dQuote{inverse of the incomplete beta function} or the quantile function of the Beta distribution as provided in \R by \code{qbeta(f, a, b)}.  Readers might consult Gilchrist (2011, chapter 12) and Karian and Dudewicz (2011, p. 510). The \eqn{pp'_r} are known in some fields as \dQuote{mean rankit} and \eqn{pp^\star_r} as \dQuote{median rankit.}
}
\usage{
pp.median(x)
}
\arguments{
  \item{x}{A real value vector. The ranks and the length of the vector are computed within the function.}
}
\value{
  An \R \code{vector} is returned.
}
\note{
The function internally calls \code{\link{pp.f}} (see Note in for that function).
}
\references{
Gilchrist, W.G., 2000, Statistical modelling with quantile functions: Chapman and Hall/CRC, Boca Raton.

Karian, Z.A., and Dudewicz, E.J., 2011, Handbook of fitting statistical distributions with R: Boca Raton, FL, CRC Press.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{pp}}, \code{\link{pp.f}} }
\examples{
\dontrun{
X <- rexp(10)*rexp(10)
means  <- pp(X, sort=FALSE)
median <- pp.median(X)
supposed.median <- pp(X, a=0.3175, sort=FALSE)
lmr <- lmoms(X)
par <- parwak(lmr)
F <- nonexceeds()
plot(F, qlmomco(F,par), type="l", log="y")
points(means, X)
points(median, X, col=2)
points(supposed.median, X, pch=16, col=2, cex=0.5)
# The plot shows that the median and supposed.median by the plotting-position
# formula are effectively equivalent. Thus, the partial application it seems
# that a=0.3175 would be good enough in lieu of the complexity of the
# quantile function of the Beta distribution.
}
}
\keyword{plotting position}
\keyword{rankit}
