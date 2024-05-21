\name{cdfgdd}
\alias{cdfgdd}
\title{Cumulative Distribution Function of the Gamma Difference Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Gamma Difference distribution (Klar, 2015) given parameters (\eqn{\alpha_1 > 0}, \eqn{\beta_1 > 0}, \eqn{\alpha_2 > 0}, \eqn{\beta_2 > 0}) computed by \code{\link{pargdd}}. The cumulative distribution function is complex and numerical integration is used.

\deqn{
F(x) = \frac{\beta_2^{\alpha_2}}{\Gamma(\alpha_1)\Gamma(\alpha_2)} \int_{\mathrm{max}\{0, -t\}}^\infty \!\!\!\!\!\!\!\!x^{\alpha_2 - 1} e^{-\beta_2x}\gamma\bigl(\alpha_1, \beta_1(x+t)\bigr)\,\mathrm{d}x\mbox{,}}

where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x \in (-\infty, \infty)}, \eqn{\Gamma(y)} is the complete gamma function, and \eqn{\gamma(a, y)} is the lower incomplete gamma function

\deqn{\gamma(a, y) = \int_0^y t^{a-1}e^{-t}\,\mathrm{d}t\mbox{.}}

The so-called Gamma Difference distribution is the distribution for the difference of two Gamma random variables \eqn{X_1 \sim \Gamma(\alpha_1, \beta_1)} and \eqn{X_1 \sim \Gamma(\alpha_2, \beta_2)}; \eqn{X = X_1 - X_2} is a Gamma Difference random variable. The distribution has other names in the literature.
}
\usage{
cdfgdd(x, para, paracheck=TRUE, silent=TRUE, ...)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{pargdd}} or  \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity.}
  \item{silent}{The argument of \code{silent} for the \code{try()} operation wrapped on \code{integrate()}.}
  \item{...}{Additional argument to pass.}
}
\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\author{W.H. Asquith}
\references{
Klar, B., 2015, A note on gamma difference distributions: Journal of Statistical Computation and Simulation v. 85, no. 18, pp. 1--8, \doi{10.1080/00949655.2014.996566}.
}
\seealso{\code{\link{pdfgdd}}, \code{\link{quagdd}}, \code{\link{lmomgdd}}, \code{\link{pargdd}}}
\examples{
\dontrun{
x <- seq(-5, 7, by=0.01)
para <- list(para=c(3,   1, 1, 1), type="gdd")
plot(x, cdfgdd(x, para), type="l", xlim=c(-5,7), ylim=c(0, 1),
     xlab="x", ylab="distribution function of gamma difference distribution")
para <- list(para=c(2,   1, 1, 1), type="gdd")
lines(x, cdfgdd(x, para), lty=2)
para <- list(para=c(1,   1, 1, 1), type="gdd")
lines(x, cdfgdd(x, para), lty=3)
para <- list(para=c(0.5, 1, 1, 1), type="gdd")
lines(x, cdfgdd(x, para), lty=4) # }
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Gamma Difference}
