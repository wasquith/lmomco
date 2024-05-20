\name{pdfgdd}
\alias{pdfgdd}
\title{Probability Density Function of the Gamma Difference Distribution}
\description{
This function computes the probability density of the Gamma Difference distribution (Klar, 2015) given parameters (\eqn{\alpha_1 > 0}, \eqn{\beta_1 > 0}, \eqn{\alpha_2 > 0}, \eqn{\beta_2 > 0}) computed by \code{\link{pargdd}}.

\deqn{f(x, x > 0) = c e^{+\beta_2x}\int_{+x}^\infty  z^{\alpha_1-1}
                      (z-x)^{\alpha_2 - 1} e^{-(\beta_1+\beta_2)z}\, \mathrm{d}z\mbox{,}}

and

\deqn{f(x, x < 0) = c e^{-\beta_1x}\int_{-x}^\infty z^{\alpha_2-1}
                      (z+x)^{\alpha_1 - 1} e^{-(\beta_1+\beta_2)z}\, \mathrm{d}z\mbox{,}}

where \eqn{c} is defined as

\deqn{c = \frac{\beta_1^{\alpha_1} \beta_2^{\alpha_2}}{\Gamma(\alpha_1) \Gamma(\alpha_2)}\mbox{,}}

where \eqn{\Gamma(y)} is the complete gamma function.
}
\usage{
pdfgdd(x, para, paracheck=TRUE, silent=TRUE, ...)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{pargdd}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity.}
  \item{silent}{The argument of \code{silent} for the \code{try()} operation wrapped on \code{integrate()}.}
  \item{...}{Additional argument to pass.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Klar, B., 2015, A note on gamma difference distributions: Journal of Statistical Computation and Simulation v. 85, no. 18, pp. 1--8, \doi{10.1080/00949655.2014.996566}.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfgdd}}, \code{\link{quagdd}}, \code{\link{lmomgdd}}, \code{\link{pargdd}}}
\examples{
x <- seq(-8, 8, by=0.01) # the operations on x are to center
para <- list(para=c(3,   1, 1, 1), type="gdd")
plot(x-(3  /1 - 1/1), pdfgdd(x, para), type="l",
     xlim=c(-6,6), ylim=c(0, 0.7),
     xlab="x", ylab="density of gamma difference distribution")
para <- list(para=c(2,   1, 1, 1), type="gdd")
lines(x-(2  /1 - 1/1), pdfgdd(x, para), lty=2)
para <- list(para=c(1,   1, 1, 1), type="gdd")
lines(x-(1  /1 - 1/1), pdfgdd(x, para), lty=3)
para <- list(para=c(0.5, 1, 1, 1), type="gdd")
lines(x-(0.5/1 - 1/1), pdfgdd(x, para), lty=4)
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Gamma Difference}
