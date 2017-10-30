\name{pwm}
\alias{pwm}
\title{ Unbiased Sample Probability-Weighted Moments }
\description{
Unbiased sample probability-weighted moments (PWMs) are computed from a sample. The \eqn{\beta_r}'s are computed using
\deqn{\beta_r = n^{-1}\sum^n_{j=1} {j-1 \choose r} x_{j:n}\mbox{.}}
}
\usage{
pwm(x, nmom=5, sort=TRUE)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{nmom}{Number of PWMs to return (\eqn{r = } \code{nmom - 1}).}
  \item{sort}{Do the data need sorting? The computations require sorted data. This option is provided to optimize processing speed if presorted data already exists.}
}
\value{
  An \R \code{list} is returned.

  \item{betas}{The PWMs. Note that convention is the have a \eqn{\beta_0}, but this is placed in the first index \code{i=1} of the \code{betas} vector.}
  \item{source}{Source of the PWMs: \dQuote{pwm}.}
}
\references{
Greenwood, J.A., Landwehr, J.M., Matalas, N.C., and Wallis, J.R., 1979, Probability weighted moments---Definition and relation to parameters of several distributions expressable in inverse form: Water Resources Research, v. 15, pp. 1,049--1,054.

Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal
of the Royal Statistical Society, Series B, v. 52, pp. 105--124.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{lmoms}}, \code{\link{pwm2lmom}}, \code{\link{pwm}} }
\examples{
# Data listed in Hosking (1995, table 29.2, p. 551)
H <- c(3,4,5,6,6,7,8,8,9,9,9,10,10,11,11,11,13,13,13,13,13,
       17,19,19,25,29,33,42,42,51.9999,52,52,52)
# 51.9999 was really 52, but a real non censored data point.
z <-  pwmRC(H,52,checkbetas=TRUE)
str(z)
# Hosking(1995) reports that A-type L-moments for this sample are
# lamA1=15.7 and lamAL-CV=.389, and lamAL-skew=.393
pwm2lmom(z$Abetas)
# WHA gets 15.666, 0.3959, and 0.4030

# See p. 553 of Hosking (1995)
# Data listed in Hosking (1995, table 29.3, p. 553)
D <- c(-2.982, -2.849, -2.546, -2.350, -1.983, -1.492, -1.443,
       -1.394, -1.386, -1.269, -1.195, -1.174, -0.854, -0.620,
       -0.576, -0.548, -0.247, -0.195, -0.056, -0.013,  0.006,
        0.033,  0.037,  0.046,  0.084,  0.221, 0.245, 0.296)
D <- c(D,rep(.2960001,40-28)) # 28 values, but Hosking mentions
                              # 40 values in total
z <-  pwmRC(D,.2960001)
# Hosking reports B-type L-moments for this sample are
# lamB1 = -.516 and lamB2 = 0.523
pwm2lmom(z$Bbetas)
# WHA gets -.5162 and 0.5218
}
\keyword{probability-weighted moment (sample)}
