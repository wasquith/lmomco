\encoding{utf8}
\name{quapdq3}
\alias{quapdq3}
\title{Quantile Function of the Polynomial Density-Quantile3 Distribution}
\description{
This function computes the quantiles of the Polynomial Density-Quantile3 distribution (PDQ3) given parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) computed by \code{\link{parpdq3}}. The quantile function is

\deqn{x(F) = \xi + \alpha \biggl[\log\biggl(\frac{F}{1-F}\biggr) + \kappa \log\bigg(\frac{[1-\kappa(2F-1)]^2}{4F(1-F)}\biggr)\biggr]\mbox{,}}

where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F},
\eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter,
and \eqn{\kappa} is a shape parameter. The range of the distribution is \eqn{-\infty < x < \infty}. This formulation of logistic distribution generalization is unique in the literature.
}
\usage{
quapdq3(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parpdq3}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\details{
The PDQ3 was proposed by Hosking (2007) with the core justification of maximizing entropy and that \dQuote{maximizing entropy subject to a set of constraints can be regarded as deriving a distribution that is consistent with the information specified in the constraints while making minimal assumptions about the form of the distribution other than those embodied in the constraints.} The PDQ3 is that family constrained to the \eqn{\lambda_1}, \eqn{\lambda_2}, and \eqn{\tau_3} values of the L-moments. (See also the Polynomial Density-Quantile4 function for constraint on \eqn{\lambda_1}, \eqn{\lambda_2}, and \eqn{\tau_4} values of the L-moments, \code{\link{quapdq4}}.)

The PDQ3 has maximum entropy conditional on having specified values for the L-moments of \eqn{\lambda_1}, \eqn{\lambda_2}, and \eqn{\lambda_3 = \tau_3\lambda_2}.  The tails of the PDQ3 are exponentially decreasing and the distribution could be useful in distributional analysis with data showing similar tail characteristics. The attainable L-kurtosis range is \eqn{\tau_4 = (5\tau_3/\kappa) - 1}.
}
\references{
Hosking, J.R.M., 2007, Distributions with maximum entropy subject to constraints on their L-moments or expected order statistics: Journal of Statistical Planning and Inference, v. 137, no. 9, pp. 2870--2891, \doi{10.1016/j.jspi.2006.10.010}.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfpdq3}}, \code{\link{pdfpdq3}}, \code{\link{lmompdq3}}, \code{\link{parpdq3}}}
\examples{
lmr <- lmoms(c(123, 34, 4, 654, 37, 78))
quapdq3(0.5, parpdq3(lmr)) # [1] 51.22802

\dontrun{
  FF <- seq(0.002475, 1 - 0.002475, by=0.001)
  para <- list(para=c(0.6933, 1.5495, 0.5488), type="pdq3")
  plot(log(FF/(1-FF)), quapdq3(FF, para), type="l", col=grey(0.8), lwd=4,
       xlab="Logistic variate, log(f/(1-f))", ylab="Quantile, Q(f)")
  lines(log(FF/(1-FF)), log(qf(FF, df1=7, df2=1)), lty=2)
  legend("topleft", c("log F(7,1) distribution with same L-moments",
                      "PDQ3 distribution with same L-moments as the log F(7,1)"),
         lwd=c(1, 4), lty=c(2, 1), col=c(1, grey(0.8)), cex=0.8)
  mtext("Mimic Hosking (2007, fig. 2 [right])") # }
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Polynomial Density-Quantile3}
