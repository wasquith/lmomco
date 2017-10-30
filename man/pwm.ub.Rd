\name{pwm.ub}
\alias{pwm.ub}
\title{ Unbiased Sample Probability-Weighted Moments }
\description{
Unbiased sample probability-weighted moments (PWMs) are computed from a sample. The \eqn{\beta_r}'s are computed using
\deqn{\beta_r = n^{-1} {n-1 \choose r}^{-1} \sum^n_{j=1} {j-1 \choose r} x_{j:n}\mbox{.}}
}
\usage{
pwm.ub(x, nmom=5, sort=TRUE)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{nmom}{Number of PWMs to return (\eqn{r =} \code{nmom - 1}).}
  \item{sort}{Do the data need sorting? The computations require sorted data. This option is provided to optimize processing speed if presorted data already exists.}
}
\value{
  An \R \code{list} is returned.

  \item{betas}{The PWMs. Note that convention is the have a \eqn{\beta_0}, but this is placed in the first index \code{i=1} of the \code{betas} vector.}
  \item{source}{Source of the PWMs: \dQuote{pwm.ub}.}
}
\note{
Through a user inquiry, it came to the author's attention in May 2014 that some unrelated studies using PWMs in the earth-system sciences have published erroneous sample PWMs formula. Because \pkg{lmomco} is intended to be an authoritative source, here are some computations to further prove correctness with provenance: 
\preformatted{
"pwm.handbookhydrology" <- function(x, nmom=5) {
   x <- sort(x, decreasing = TRUE); n <- length(x); betas <- rep(NA, nmom)
   for(r in 0:(nmom-1)) {
      tmp <- sum(sapply(1:(n-r),
          function(j) { choose(n - j, r) * x[j] / choose(n - 1, r) }))
      betas[(r+1)] <- tmp/n
   }
   return(betas)
}
}
and a demonstration with alternative algebra in Stedinger and others (1993)
\preformatted{
set.seed(1)
glo <- vec2par(c(123,1123,-.5), type="glo"); X <- rlmomco(100, glo)
lmom2pwm(lmoms(X, nmom=5))$betas # unbiased L-moments flipped to PWMs
[1]  998.7932 1134.0658 1046.4906  955.8872  879.3349
pwm.ub(X, nmom=5)$betas  # Hosking and Wallis (1997) and Asquith (2011)
[1]  998.7932 1134.0658 1046.4906  955.8872  879.3349
pwm.handbookhydrology(X) # ** alert reverse sort, opposite usually seen**
[1]  998.7932 1134.0658 1046.4906  955.8872  879.3349
}
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Greenwood, J.A., Landwehr, J.M., Matalas, N.C., and Wallis, J.R., 1979, Probability weighted moments---Definition and relation to parameters of several distributions expressable in inverse form: Water Resources Research, v. 15, pp. 1,049--1,054.

Stedinger, J.R., Vogel, R.M., Foufoula-Georgiou, E., 1993, Frequency analysis of extreme events: \emph{in} Handbook of Hydrology, ed. Maidment, D.R., McGraw-Hill, Section 18.6 Partial duration series, mixtures, and censored data, pp. 18.37--18.39.

Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{pwm.pp}}, \code{\link{pwm.gev}}, \code{\link{pwm2lmom}} }
\examples{
pwm <- pwm.ub(rnorm(20))
}
\keyword{probability-weighted moment (sample)}
