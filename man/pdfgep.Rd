\name{pdfgep}
\alias{pdfgep}
\title{Probability Density Function of the Generalized Exponential Poisson Distribution}
\description{
This function computes the probability density of the Generalized Exponential Poisson distribution given parameters (\eqn{\beta}, \eqn{\kappa}, and \eqn{h})  computed by \code{\link{pargep}}. The probability density function is
\deqn{f(x) = \frac{\kappa h \eta}{[1 - \exp(-h)]^\kappa}{1 - \exp[-h + h\exp(-\eta x)}\times\exp[-h - \eta x + h\exp(-\eta x)]\mbox{,}}
where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x > 0}, \eqn{\eta = 1/\beta}, \eqn{\beta > 0} is a scale parameter, \eqn{\kappa > 0} is a shape parameter, and \eqn{h > 0} is another shape parameter.
}
\usage{
pdfgep(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{pargep}} or \code{\link{vec2par}}.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Barreto-Souza, W., and Cribari-Neto, F., 2009, A generalization of the exponential-Poisson distribution: Statistics and Probability, 79, pp. 2493--2500.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdfgep}}, \code{\link{quagep}}, \code{\link{lmomgep}}, \code{\link{pargep}}}
\examples{
pdfgep(0.5, vec2par(c(10,2.9,1.5), type="gep"))
\dontrun{
x <- seq(0,3, by=0.01); ylim <- c(0,1.5)
plot(NA,NA, xlim=range(x), ylim=ylim, xlab="x", ylab="f(x)")
mtext("Barreto-Souza and Cribari-Neto (2009, fig. 1)")
K <- c(0.1, 1, 5, 10)
for(i in 1:length(K)) {
   gep <- vec2par(c(2,K[i],1), type="gep"); lines(x, pdfgep(x, gep), lty=i)
}
}
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Generalized Exponential Poisson}
