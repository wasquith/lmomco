\name{pdftexp}
\alias{pdftexp}
\title{Probability Density Function of the Truncated Exponential Distribution}
\description{
This function computes the probability density of the Truncated Exponential distribution given parameters (\eqn{\psi} and \eqn{\alpha})  computed by \code{\link{partexp}}. The parameter \eqn{\psi} is the right truncation, and \eqn{\alpha} is a scale parameter. The probability density function, letting \eqn{\beta = 1/\alpha} to match nomenclature of Vogel and others (2008), is
\deqn{f(x) = \frac{\beta\,\exp(-\beta{t})}{1 - \mathrm{exp}(-\beta\psi)}\mbox{,}}
where \eqn{x(x)} is the probability density for the quantile \eqn{0 \le x \le \psi} and \eqn{\psi > 0} and \eqn{\alpha > 0}. This distribution represents a nonstationary Poisson process.

The distribution is restricted to a narrow range of L-CV (\eqn{\tau_2 = \lambda_2/\lambda_1}). If \eqn{\tau_2 = 1/3}, the process represented is a stationary Poisson for which the probability density function is simply the uniform distribution and \eqn{f(x) = 1/\psi}. If \eqn{\tau_2 = 1/2}, then the distribution is represented as the usual exponential distribution with a location parameter of zero and a scale parameter \eqn{1/\beta}. Both of these limiting conditions are supported.
}
\usage{
pdftexp(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{partexp}} or \code{\link{vec2par}}.}
}

\value{
  Probability density (\eqn{F}) for \eqn{x}.
}
\references{
Vogel, R.M., Hosking, J.R.M., Elphick, C.S., Roberts, D.L., and Reed, J.M., 2008, Goodness of fit of probability distributions for sightings as species approach extinction: Bulletin of Mathematical Biology, DOI 10.1007/s11538-008-9377-3, 19 p.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdftexp}}, \code{\link{quatexp}}, \code{\link{lmomtexp}}, \code{\link{partexp}}}
\examples{
lmr <- vec2lmom(c(40,0.38), lscale=FALSE)
pdftexp(0.5,partexp(lmr))
\dontrun{
F <- seq(0,1,by=0.001)
A <- partexp(vec2lmom(c(100, 1/2), lscale=FALSE))
x <- quatexp(F, A)
plot(x, pdftexp(x, A), pch=16, type='l')
by <- 0.01; lcvs <- c(1/3, seq(1/3+by, 1/2-by, by=by), 1/2)
reds <- (lcvs - 1/3)/max(lcvs - 1/3)
for(lcv in lcvs) {
    A <- partexp(vec2lmom(c(100, lcv), lscale=FALSE))
    x <- quatexp(F, A)
    lines(x, pdftexp(x, A),
          pch=16, col=rgb(reds[lcvs == lcv],0,0))
}
}
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Exponential (trimmed)}
