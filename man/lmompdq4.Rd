\encoding{utf8}
\name{lmompdq4}
\alias{lmompdq4}
\title{L-moments of the Polynomial Density-Quantile4 Distribution}
\description{
This function estimates the L-moments of the Polynomial Density-Quantile4 distribution given the parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) from \code{\link{parpdq4}}. The L-moments in terms of the parameters are
\deqn{\lambda_1 = \xi\mbox{,}}
\deqn{\lambda_2 = \frac{\alpha}{\kappa} (1-\kappa^2)\, \mathrm{atanh}(\kappa)\mathrm{\ for\ } \kappa > 0\mbox{,}}
\deqn{\lambda_2 = \frac{\alpha}{\kappa} (1+\kappa^2)\, \mathrm{atan}(\kappa)\mathrm{\ for\ } \kappa < 0\mbox{,}}
\deqn{\tau_3 = 0 \mbox{, and}}
\deqn{\tau_4 = -\frac{1}{4} + \frac{5}{4\kappa}\biggl(\frac{1}{\kappa} - \frac{1}{\mathrm{atanh}(\kappa)} \biggr) \mathrm{\ for\ } \kappa > 0\mbox{,}}
\deqn{\tau_4 = -\frac{1}{4} - \frac{5}{4\kappa}\biggl(\frac{1}{\kappa} - \frac{1}{\mathrm{atan}(\kappa)} \biggr) \mathrm{\ for\ } \kappa < 0\mbox{,}}
}
\usage{
lmompdq4(para, paracheck=TRUE)
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{paracheck}{A logical switch as to whether the validity of the parameters should be checked. Default is \code{paracheck=TRUE}.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{0}.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL}.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL}.}
  \item{source}{An attribute identifying the computational
source of the L-moments: \dQuote{lmompdq4}.}
}
\references{
Hosking, J.R.M., 2007, Distributions with maximum entropy subject to constraints on their L-moments or expected order statistics: Journal of Statistical Planning and Inference, v. 137, no. 9, pp. 2,870--2,891, \doi{10.1016/j.jspi.2006.10.010}.
}
\author{W.H. Asquith}
\seealso{\code{\link{parpdq4}}, \code{\link{cdfpdq4}}, \code{\link{pdfpdq4}}, \code{\link{quapdq4}} }
\examples{
\dontrun{
para <- list(para=c(20, 1, -0.5), type="pdq4")
lmoms(quapdq4(runif(100000), para))$lambdas
lmompdq4(para)$lambdas #}

\dontrun{
para <- list(para=c(20, 1, +0.5), type="pdq4")
lmoms(quapdq4(runif(100000), para))$lambdas
lmompdq4(para)$lambdas #}

\dontrun{
K1 <- seq(-5, 0, by=0.001)
K2 <- seq( 0, 1, by=0.001)
suppressWarnings(mono_decrease_part1 <- -(1/4) + (5/(4*K1)) * (1/K1 - 1/atanh(K1)))
                 mono_increase_part2 <- -(1/4) - (5/(4*K1)) * (1/K1 - 1/atan( K1))
                 mono_increase_part1 <- -(1/4) + (5/(4*K2)) * (1/K2 - 1/atanh(K2))
                 mono_decrease_part2 <- -(1/4) - (5/(4*K2)) * (1/K2 - 1/atan( K2))

plot( 0, 0, type="n", xlim=range(c(K1,K2)), ylim=c(-0.25, 1),
     xlab="Kappa shape parameter PDQ4 distribution", ylab="L-kurtosis (Tau4)")
lines(K1, mono_decrease_part1, col=4, lwd=0.3)
lines(K2, mono_increase_part1, col=4, lwd=3)
lines(K2, mono_decrease_part2, col=2, lwd=0.3)
lines(K1, mono_increase_part2, col=2, lwd=3)

abline(h= 1/6, lty=2, lwd=0.6)
abline(h=-1/4, lty=2, lwd=0.6)
text(-5, -1/4, "Tau4 lower bounds", pos=4, cex=0.8)
abline(v=0,    lty=2, lwd=0.6)
abline(v=1,    lty=1, lwd=0.9)
points(-0.7029, 0.1226, pch=15, col="darkgreen")

legend("topleft", c("Monotonic increasing for kappa < 0",
                    "Monotonic increasing for kappa > 0",
                    "Monotonic decreasing for kappa > 0 (not used)",
                    "Monotonic decreasing for kappa < 0 (not used)",
                    "Normal distribution"), cex=0.8,
        pch=c(NA, NA, NA, NA, 15), lwd=c(3,3, 0.3, 0.3, NA),
        pt.cex=c(NA, NA, NA, NA, 1), col=c(2, 4, 2, 4, "darkgreen")) # }
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Polynomial Density-Quantile4}
