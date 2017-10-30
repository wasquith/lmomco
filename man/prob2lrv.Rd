\name{prob2lrv}
\alias{prob2lrv}
\title{Convert a Vector of Annual Nonexceedance Probabilities to Logistic Reduced Variates}
\description{
This function converts a vector of annual nonexceedance probabilities \eqn{F} to logistic reduced variates (LRV, \eqn{lrv}) 
\deqn{lrv = 1/(\exp(-lrv) + 1)\mbox{,}}
where \eqn{0 \le F \le 1}. The logistic distribution, which is generalized by the Generalized Logistic (\code{\link{quaglo}}) with \eqn{\kappa = 0}, will plot as a straightline when the horizontal axis is LRV transformed.
}
\usage{
prob2lrv(f)
}
\arguments{
  \item{f}{A vector of annual nonexceedance probabilities.}
}
\value{
   A vector of logistic reduced variates.
}
\author{ W.H. Asquith}
\references{
Bradford, R.B., 2002, Volume-duration growth curves for flood estimation in permeable catchments: Hydrology and Earth System Sciences, v. 6, no. 5, pp. 939--947.
}
\seealso{\code{\link{lrv2prob}}, \code{\link{prob2T}}
}
\examples{
F <- nonexceeds()
lrv <- prob2lrv(F)
\dontrun{
X <- rlmomco(10040, vec2par(c(0,1,0), type="glo"))
plot(prob2lrv(pp(X, a=0.4)), sort(X)); abline(0,1)
}
}
\keyword{utility (probability)}
\concept{Logistic reduced variate}
\concept{Logistic Reduced Variate}
\concept{logistic reduced variate}
