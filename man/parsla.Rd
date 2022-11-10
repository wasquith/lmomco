\name{parsla}
\alias{parsla}
\title{Estimate the Parameters of the Slash Distribution}
\description{
This function estimates the parameters of the Slash distribution from the trimmed L-moments (TL-moments) having trim level 1. The relations between distribution parameters and TL-moments are shown under \code{\link{lmomsla}}.
}
\usage{
parsla(lmom, ...)
}
\arguments{
  \item{lmom}{A TL-moment object from \code{\link{TLmoms}} with \code{trim=1}.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{sla}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parsla}.}
}
\references{
Rogers, W.H., and Tukey, J.W., 1972, Understanding some long-tailed symmetrical distributions: Statistica Neerlandica, v. 26, no. 3, pp. 211--226.
}
\author{W.H. Asquith}
\seealso{\code{\link{TLmoms}}, \code{\link{lmomsla}}, \code{\link{cdfsla}}, \code{\link{pdfsla}}, \code{\link{quasla}}
}
\examples{
\dontrun{
par1 <- vec2par(c(-100, 30), type="sla")
X   <- rlmomco(500, par1)
lmr <- TLmoms(X, trim=1)
par2 <- parsla(lmr)
F <- seq(0.001,.999, by=0.001)
plot(qnorm(pp(X)), sort(X), pch=21, col=8,
     xlab="STANDARD NORMAL VARIATE",
     ylab="QUANTILE")
lines(qnorm(F), quasla(F, par1), lwd=3)
lines(qnorm(F), quasla(F, par2), col=2)
}
}
\keyword{distribution (parameters)}
\keyword{Distribution: Slash}
