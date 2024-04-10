\name{quaaep4}
\alias{quaaep4}
\title{Quantile Function of the 4-Parameter Asymmetric Exponential Power Distribution}
\description{
This function computes the quantiles of the 4-parameter Asymmetric Exponential Power distribution given parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\kappa}, and \eqn{h}) of the distribution computed by \code{\link{paraep4}}. The quantile function of the distribution given the cumulative distribution function \eqn{F(x)} for \eqn{F < F(\xi)} is
\deqn{x(F) = \xi - \alpha\kappa\biggl[\gamma^{(-1)}\bigl((1+\kappa^2)F/\kappa^2,\; 1/h\bigr)\biggr]^{1/h}\mbox{,}}
and for \eqn{F \ge F(\xi)} is
\deqn{x(F) = \xi + \frac{\alpha}{\kappa}\biggl[\gamma^{(-1)}\bigl((1+\kappa^2)(1-F),\; 1/h\bigr)\biggr]^{1/h} \mbox{,}}
where \eqn{x(F)} is the quantile \eqn{x} for nonexceedance probability \eqn{F},
\eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter,
\eqn{\kappa} is a shape parameter, \eqn{h} is another shape parameter, \eqn{\gamma^{(-1)}(Z, shape)} is the inverse of the upper tail of the incomplete gamma function. The range of the distribution is \eqn{-\infty < x < \infty}. The inverse upper tail of the incomplete gamma function is \code{qgamma(Z, shape, lower.tail=FALSE)} in \R. The mathematical definition of the upper tail of the incomplete gamma function shown in documentation for \code{\link{cdfaep4}}. If the \eqn{\tau_3} of the distribution is zero (symmetrical), then the distribution is known as the Exponential Power (see \code{\link{lmrdia46}}).
}
\usage{
quaaep4(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{paraep4}} or similar.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\author{W.H. Asquith}
\references{
Asquith, W.H., 2014, Parameter estimation for the 4-parameter asymmetric exponential power distribution by the method of L-moments using R: Computational Statistics and Data Analysis, v. 71, pp. 955--970.

Delicado, P., and Goria, M.N., 2008, A small sample comparison of maximum likelihood,
moments and L-moments methods for the asymmetric exponential power distribution:
Computational Statistics and Data Analysis, v. 52, no. 3, pp. 1661--1673.
}
\seealso{\code{\link{cdfaep4}}, \code{\link{pdfaep4}}, \code{\link{lmomaep4}}, \code{\link{paraep4}}}
\examples{
para <- vec2par(c(0,1, 0.5, 2), type="aep4");
IQR <- quaaep4(0.75,para) - quaaep4(0.25,para);
cat("Interquartile Range=",IQR,"\n")

\dontrun{
F <- c(0.00001, 0.0001, 0.001, seq(0.01, 0.99, by=0.01),
       0.999, 0.9999, 0.99999);
delx <- 0.1;
x <- seq(-10,10, by=delx);
K <- .67

PAR <- list(para=c(0,1, K, 0.5), type="aep4");
plot(x,cdfaep4(x, PAR), type="n",
     ylab="NONEXCEEDANCE PROBABILITY",
     ylim=c(0,1), xlim=c(-20,20));
lines(x,cdfaep4(x,PAR), lwd=3);
lines(quaaep4(F, PAR), F, col=4);

PAR <- list(para=c(0,1, K, 1), type="aep4");
lines(x,cdfaep4(x, PAR), lty=2, lwd=3);
lines(quaaep4(F, PAR), F, col=4, lty=2);

PAR <- list(para=c(0,1, K, 2), type="aep4");
lines(x,cdfaep4(x, PAR), lty=3, lwd=3);
lines(quaaep4(F, PAR), F, col=4, lty=3);

PAR <- list(para=c(0,1, K, 4), type="aep4");
lines(x,cdfaep4(x, PAR), lty=4, lwd=3);
lines(quaaep4(F, PAR), F, col=4, lty=4);
}
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Asymmetric Exponential Power}
\keyword{Distribution: Exponential Power}
