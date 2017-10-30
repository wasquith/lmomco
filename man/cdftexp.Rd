\name{cdftexp}
\alias{cdftexp}
\title{Cumulative Distribution Function of the Truncated Exponential Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Truncated Exponential distribution given parameters (\eqn{\psi} and \eqn{\alpha}) computed by \code{\link{partexp}}. The parameter \eqn{\psi} is the right truncation of the distribution and \eqn{\alpha} is a scale parameter. The cumulative distribution function, letting \eqn{\beta = 1/\alpha} to match nomenclature of Vogel and others (2008), is
\deqn{F(x) = \frac{1-\mathrm{exp}(-\beta{t})}{1-\mathrm{exp}(-\beta\psi)}\mbox{,}}
where \eqn{F(x)} is the nonexceedance probability for the quantile \eqn{0 \le x \le \psi} and \eqn{\psi > 0} and \eqn{\alpha > 0}. This distribution represents a nonstationary Poisson process.

The distribution is restricted to a narrow range of L-CV (\eqn{\tau_2 = \lambda_2/\lambda_1}). If \eqn{\tau_2 = 1/3}, the process represented is a stationary Poisson for which the cumulative distribution function is simply the uniform distribution and \eqn{F(x) = x/\psi}. If \eqn{\tau_2 = 1/2}, then the distribution is represented as the usual exponential distribution with a location parameter of zero and a rate parameter \eqn{\beta} (scale parameter \eqn{\alpha = 1/\beta}). These two limiting conditions are supported.
}
\usage{
cdftexp(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{partexp}} or  \code{\link{vec2par}}.}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Vogel, R.M., Hosking, J.R.M., Elphick, C.S., Roberts, D.L., and Reed, J.M., 2008, Goodness of fit of probability distributions for sightings as species approach extinction: Bulletin of Mathematical Biology, DOI 10.1007/s11538-008-9377-3, 19 p.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdftexp}}, \code{\link{quatexp}}, \code{\link{lmomtexp}}, \code{\link{partexp}}}
\examples{
cdftexp(50,partexp(vec2lmom(c(40,0.38), lscale=FALSE)))
\dontrun{
F <- seq(0,1,by=0.001)
A <- partexp(vec2lmom(c(100, 1/2), lscale=FALSE))
x <- quatexp(F, A)
plot(x, cdftexp(x, A), pch=16, type='l')
by <- 0.01; lcvs <- c(1/3, seq(1/3+by, 1/2-by, by=by), 1/2)
reds <- (lcvs - 1/3)/max(lcvs - 1/3)
for(lcv in lcvs) {
    A <- partexp(vec2lmom(c(100, lcv), lscale=FALSE))
    x <- quatexp(F, A)
    lines(x, cdftexp(x, A), pch=16, col=rgb(reds[lcvs == lcv],0,0))
}

  # Vogel and others (2008) example sighting times for the bird
  # Eskimo Curlew, inspection shows that these are fairly uniform.
  # There is a sighting about every year to two.
  T <- c(1946, 1947, 1948, 1950, 1955, 1956, 1959, 1960, 1961,
         1962, 1963, 1964, 1968, 1970, 1972, 1973, 1974, 1976,
         1977, 1980, 1981, 1982, 1982, 1983, 1985)
  R <- 1945 # beginning of record
  S <- T - R
  lmr <- lmoms(S)
  PARcurlew <- partexp(lmr)
  # read the warning message and then force the texp to the
  # stationary process model (min(tau_2) = 1/3).
  lmr$ratios[2] <- 1/3
  lmr$lambdas[2] <- lmr$lambdas[1]*lmr$ratios[2]
  PARcurlew <- partexp(lmr)
  Xmax <- quatexp(1, PARcurlew)
  X <- seq(0,Xmax, by=.1)
  plot(X, cdftexp(X,PARcurlew), type="l")
  # or use the MVUE estimator
  TE <- max(S)*((length(S)+1)/length(S)) # Time of Extinction
  lines(X, punif(X, min=0, max=TE), col=2)}
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Exponential (trimmed)}
