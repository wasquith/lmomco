\name{quagep}
\alias{quagep}
\title{Quantile Function of the Generalized Exponential Poisson Distribution}
\description{
This function computes the quantiles of the Generalized Exponential Poisson distribution given parameters (\eqn{\beta}, \eqn{\kappa}, and \eqn{h}) of the distribution computed by \code{\link{pargep}}. The quantile function of the distribution is
\deqn{x(F) = \eta^{-1} \log[1 + h^{-1}\log(1 - F^{1/\kappa}[1 - \exp(-h)])]\mbox{,}}
where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x > 0}, \eqn{\eta = 1/\beta}, \eqn{\beta > 0} is a scale parameter, \eqn{\kappa > 0} is a shape parameter, and \eqn{h > 0} is another shape parameter.
}
\usage{
quagep(f, para, paracheck=TRUE)
}
\details{If \code{f = 1} or is so close to unity that \code{NaN} in the computations of the quantile function, then the function enters into an infinite loop for which an \dQuote{order of magnitude decrement} on the value of \cr \code{.Machine$double.eps} is made until a numeric hit is encountered.  Let \eqn{\eta} be this machine value, then \eqn{F = 1 - \eta^{1/j}} where \eqn{j} is the iteration in the infinite loop. Eventually \eqn{F} becomes small enough that a finite value will result. This result is an estimate of the maximum numerical value the function can produce on the current running platform.  This feature assists in the numerical integration of the quantile function for L-moment estimation (see \code{\link{expect.max.ostat}}).  The \code{\link{expect.max.ostat}} was zealous on reporting errors related to lack of finite integration. However with the \dQuote{order magnitude decrementing,} then the errors in \code{\link{expect.max.ostat}} become fewer and are either
\preformatted{
Error in integrate(fnb, lower, upper, subdivisions = 200L) : 
  extremely bad integrand behaviour
}
or
\preformatted{
Error in integrate(fnb, lower, upper, subdivisions = 200L) : 
  maximum number of subdivisions reached
}
and are shown here to aid in research into Generalized Exponential Power implementation.
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{pargep}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Barreto-Souza, W., and Cribari-Neto, F., 2009, A generalization of the exponential-Poisson distribution: Statistics and Probability, 79, pp. 2493--2500.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfgep}}, \code{\link{pdfgep}}, \code{\link{lmomgep}}, \code{\link{pargep}}}
\examples{
gep <- list(para=c(2, 1.5, 3), type="gep")
quagep(0.5, gep)
\dontrun{
  pdf("gep.pdf")
  F <- nonexceeds(f01=TRUE)
  K <- seq(-1,2,by=.2); H <- seq(-1,2,by=.2)
  K <- 10^(K); H <- 10^(H)
  for(i in 1:length(K)) {
    for(j in 1:length(H)) {
      gep <- vec2par(c(2,K[i],H[j]), type="gep")
      message("(K,H): ",K[i]," ",H[j])
      plot(F, quagep(F, gep), lty=i, col=j, type="l", ylim=c(0,4),
           xlab="NONEXCEEDANCE PROBABILITY", ylab="X(F)")
      mtext(paste("(K,H): ",K[i]," ",H[j]))
    }
  }
  dev.off()
}
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Generalized Exponential Poisson}
