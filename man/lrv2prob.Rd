\name{lrv2prob}
\alias{lrv2prob}
\title{Convert a Vector of Logistic Reduced Variates to Annual Nonexceedance Probabilities}
\description{
This function converts a vector of logistic reduced variates (\eqn{lrv}) to annual nonexceedance probabilities \eqn{F}
\deqn{F = -\log((1-lrv)/lrv)\mbox{,}}
where \eqn{0 \le F \le 1}.
}
\usage{
lrv2prob(lrv)
}
\arguments{
  \item{lrv}{A vector of logistic reduced variates.}
}
\value{
   A vector of annual nonexceedance probabilities.
}
\author{ W.H. Asquith}
\references{
Bradford, R.B., 2002, Volume-duration growth curves for flood estimation in permeable catchments: Hydrology and Earth System Sciences, v. 6, no. 5, pp. 939--947.
}
\seealso{\code{\link{prob2lrv}}, \code{\link{prob2T}}
}
\examples{
T <- c(1, 2, 5, 10, 25, 50, 100, 250, 500); lrv <- prob2grv(T2prob(T))
F <- lrv2prob(lrv)
}
\keyword{utility (probability)}
\concept{Logistic reduced variate}
\concept{Logistic Reduced Variate}
\concept{logistic reduced variate}
