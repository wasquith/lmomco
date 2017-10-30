\name{quatexp}
\alias{quatexp}
\title{Quantile Function of the Truncated Exponential Distribution}
\description{
This function computes the quantiles of the Truncated Exponential distribution given parameters (\eqn{\psi} and \eqn{\alpha})  computed by \code{\link{partexp}}. The parameter \eqn{\psi} is the right truncation, and \eqn{\alpha} is a scale parameter. The quantile function, letting \eqn{\beta = 1/\alpha} to match nomenclature of Vogel and others (2008), is
\deqn{x(F) = -\frac{1}{\beta}\log(1-F[1-\mathrm{exp}(-\beta\psi)])\mbox{,}}
where \eqn{x(F)} is the quantile \eqn{0 \le x \le \psi} for nonexceedance probability \eqn{F} and \eqn{\psi > 0} and \eqn{\alpha > 0}. This distribution represents a nonstationary Poisson process.

The distribution is restricted to a narrow range of L-CV (\eqn{\tau_2 = \lambda_2/\lambda_1}). If \eqn{\tau_2 = 1/3}, the process represented is a stationary Poisson for which the quantile function is simply the uniform distribution and \eqn{x(F) = \psi\,F}. If \eqn{\tau_2 = 1/2}, then the distribution is represented as the usual exponential distribution with a location parameter of zero and a scale parameter \eqn{1/\beta}. Both of these limiting conditions are supported.
}
\usage{
quatexp(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{partexp}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Vogel, R.M., Hosking, J.R.M., Elphick, C.S., Roberts, D.L., and Reed, J.M., 2008, Goodness of fit of probability distributions for sightings as species approach extinction: Bulletin of Mathematical Biology, DOI 10.1007/s11538-008-9377-3, 19 p.
}
\author{ W.H. Asquith}
\seealso{\code{\link{cdftexp}}, \code{\link{pdftexp}}, \code{\link{lmomtexp}}, \code{\link{partexp}}}
\examples{
lmr <- vec2lmom(c(40,0.38), lscale=FALSE)
quatexp(0.5,partexp(lmr))
\dontrun{
F <- seq(0,1,by=0.001)
A <- partexp(vec2lmom(c(100, 1/2), lscale=FALSE))
plot(qnorm(F), quatexp(F, A), pch=16, type='l')
by <- 0.01; lcvs <- c(1/3, seq(1/3+by, 1/2-by, by=by), 1/2)
reds <- (lcvs - 1/3)/max(lcvs - 1/3)
for(lcv in lcvs) {
    A <- partexp(vec2lmom(c(100, lcv), lscale=FALSE))
    lines(qnorm(F), quatexp(F, A), col=rgb(reds[lcvs == lcv],0,0))
}
}
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Exponential (trimmed)}
